#
# Client for retrieving conserved-domain data for SEED proteins.
#
# Three backends, selected in this order by new():
#
#   rpsblast  local rps-blast against a CDD profile database. Preferred:
#             ~0.06s/protein warm, no network, no rate limit. Enabled by
#             setting $FIG_Config::cdd_rps_db.
#   ncbi      NCBI Batch CD-Search (bwrpsb.cgi). Submit/poll; a cold lookup
#             has been measured at anywhere from 20s to over 200s, so this
#             path is asynchronous and may report results as still pending.
#   jsonrpc   the original ConservedDomainSearch JSON-RPC service, if
#             $FIG_Config::ConservedDomainSearchURL is set.
#
# Results are cached per protein md5 under the owning genome's directory, so
# repeat views cost nothing and a re-annotation that changes a sequence
# naturally misses the cache.
#
# Configuration (all optional except cdd_rps_db for the local backend):
#
#   $FIG_Config::cdd_rps_db        path prefix of the rps database, e.g.
#                                  "/scratch/olson/cdd/Cdd" (the directory
#                                  holding Cdd.pal and its volumes).
#   $FIG_Config::cdd_rps_prog      rpsblast executable; searched on PATH if unset.
#   $FIG_Config::cdd_rps_threads   -num_threads for rpsblast (default 4).
#   $FIG_Config::cdd_evalue        E-value cutoff (default 0.01).
#   $FIG_Config::cdd_wait_budget   seconds to wait on NCBI before returning
#                                  pending (default 10).
#   $FIG_Config::ConservedDomainSearchURL   legacy JSON-RPC endpoint.
#

package ConservedDomainSearch;

use strict;
use warnings;

use FIG_Config;
use Data::Dumper;
use JSON::XS;
use LWP::UserAgent;
use Digest::MD5 'md5_hex';
use File::Path 'make_path';
use File::Temp;
use SeedUtils;
use base 'Class::Accessor';

__PACKAGE__->mk_accessors(qw(fig url ua backend rps_db rps_prog rps_threads evalue));

our $idx = 1;

our $NCBI_URL = "https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi";

#
# A hit row has eleven columns. The first ten are exactly the NCBI Batch
# CD-Search "hits" table with its leading Query column stripped, which is the
# layout create_cdd_features() and domains_of() already index into. The
# eleventh is our own addition: the domain's long description, which NCBI
# serves from a separate pssmid_lookup call but which rpsblast hands us for
# free in the subject title.
#
use constant {
    COL_TYPE        => 0,
    COL_PSSMID      => 1,
    COL_FROM        => 2,
    COL_TO          => 3,
    COL_EVALUE      => 4,
    COL_BITSCORE    => 5,
    COL_ACCESSION   => 6,
    COL_SHORTNAME   => 7,
    COL_INCOMPLETE  => 8,
    COL_SUPERFAMILY => 9,
    COL_DESCRIPTION => 10,
};

sub new
{
    my($class, $fig) = @_;

    my $self = {
        fig         => $fig,
        url         => $FIG_Config::ConservedDomainSearchURL,
        ua          => LWP::UserAgent->new(),
        rps_db      => $FIG_Config::cdd_rps_db,
        rps_prog    => $FIG_Config::cdd_rps_prog,
        rps_threads => (defined($FIG_Config::cdd_rps_threads) ? $FIG_Config::cdd_rps_threads : 4),
        evalue      => (defined($FIG_Config::cdd_evalue) ? $FIG_Config::cdd_evalue : 0.01),
        _mem        => {},
    };

    #
    # A 3600 second timeout (what this module used to carry) means a wedged
    # service pins a web worker for an hour. Nothing here should take a minute.
    #
    $self->{ua}->timeout(60);
    $self->{ua}->agent("SEED-ConservedDomainSearch/1.0");

    bless $self, $class;

    if ($self->{rps_db})
    {
        my $prog = $self->{rps_prog} || _find_rps_prog();
        if ($prog && -x $prog && _rps_db_present($self->{rps_db}))
        {
            $self->{rps_prog} = $prog;
            $self->{backend}  = 'rpsblast';
        }
        else
        {
            warn "ConservedDomainSearch: cdd_rps_db is set to '$self->{rps_db}' but " .
                 ($prog ? "the database is not readable" : "no rpsblast executable was found") .
                 "; falling back\n";
        }
    }

    $self->{backend} ||= $self->{url} ? 'jsonrpc' : 'ncbi';

    return $self;
}

sub _find_rps_prog
{
    for my $name (qw(rpsblast+ rpsblast))
    {
        for my $dir (split(/:/, $ENV{PATH} || ''))
        {
            my $p = "$dir/$name";
            return $p if -x $p;
        }
    }
    return undef;
}

#
# A multi-volume database has <db>.pal; a single-volume one has <db>.rps.
#
sub _rps_db_present
{
    my($db) = @_;
    return (-f "$db.pal" || -f "$db.rps") ? 1 : 0;
}

=head3 create_cdd_features

    my @features = $cdd->create_cdd_features($fid, $options);

Look up the given fid and create a set of quasi-features whose locations are
mapped from protein coordinates back onto the contig. Each returned element is
C<[$sfid, $type, $annotation, $location, $translation]>.

=cut

sub create_cdd_features
{
    my($self, $fid, $options) = @_;

    my $cdd = $self->lookup($fid, $options);

    $cdd = $cdd->{$fid};
    return () unless $cdd;

    my $loc = $self->fig->feature_location($fid);
    return () unless $loc;

    #
    # We are going to simplify this code by assuming contiguous locations.
    # It's just an approximation anyway.
    #

    my($contig, $left, $right, $strand) = SeedUtils::boundaries_of($loc);

    my $translation = $self->fig->get_translation($fid);

    my @out;
    my $subid = 1;

    #
    # Accession, short name and description all travel in the hit row, so the
    # pssmid_lookup round trip the original did here is gone.
    #
    my %info;
    for my $what (qw(domain_hits site_annotations structural_motifs))
    {
        for my $ent (@{$cdd->{$what} || []})
        {
            my $pssmid = $what eq 'domain_hits'       ? $ent->[1]
                       : $what eq 'structural_motifs' ? $ent->[3]
                       :                                $ent->[5];
            next unless defined($pssmid);
            next if $info{$pssmid};
            my $acc   = $ent->[COL_ACCESSION];
            my $short = $ent->[COL_SHORTNAME];
            my $desc  = $ent->[COL_DESCRIPTION];
            $info{$pssmid} = [$acc, $short, $desc]
                if defined($acc) && $acc ne '';
        }
    }

    for my $what (qw(domain_hits site_annotations structural_motifs))
    {
        my $list = $cdd->{$what};
        my $pssmid;
        for my $ent (@$list)
        {
            my @locs;
            my $anno;

            if ($what eq 'domain_hits')
            {
                @locs = ([$ent->[2], $ent->[3]]);
                $pssmid = $ent->[1];
                $anno = $ent->[7];
            }
            elsif ($what eq 'structural_motifs')
            {
                @locs = ([$ent->[1], $ent->[2]]);
                $anno = $ent->[0];
                $pssmid = $ent->[3];
            }
            else
            {
                $anno = $ent->[1];
                $pssmid = $ent->[5];
                my @x = split(/,/, $ent->[2]);
                for my $i (0..$#x)
                {
                    my($v) = $x[$i] =~ /(\d+)/;
                    push(@locs, [$v, $v, $i+1]);
                }
            }

            #
            # Now do the coordinate mapping & create features for each.
            #

            for my $loc (@locs)
            {
                my($pstart, $pend, $xidx) = @$loc;
                next unless defined($pstart) && defined($pend);
                my($lbeg, $lend);
                if ($strand eq '+')
                {
                    #
                    # Residue N occupies bases (N-1)*3 .. N*3-1 counting from
                    # $left. The original omitted the -1 and so pushed every
                    # plus-strand domain three bases downstream, which at
                    # residue 1 produced lbeg > lend -- an inverted location.
                    #
                    $lbeg = $left + ($pstart - 1) * 3;
                    $lend = $left + $pend * 3 - 1;
                }
                else
                {
                    $lbeg = $right - ($pstart - 1) * 3;
                    $lend = $right - $pend * 3 + 1;
                }

                my $floc = join("_", $contig, $lbeg, $lend);
                my $trans = substr($translation, $pstart - 1, $pend - $pstart + 1);
                my $sfid = "$fid.$pstart-$pend.$what.$subid";
                my $info = $info{$pssmid};
                if ($info)
                {
                    $sfid = $info->[0] . "-$fid-$subid";
                    $anno = $info->[1];
                    $anno .= ": " . $info->[2] if defined($info->[2]) && $info->[2] ne '';
                    $sfid .= "-$xidx" if (defined($xidx));
                }
                my $type = $what;
                $type =~ s/s$//;
                #
                # The accession rides along as a sixth element so callers can
                # group identical domains without having to pick the id apart.
                #
                push(@out, [$sfid, $type, $anno, $floc, $trans,
                            ($info ? $info->[0] : undef)]);
                $subid++;
            }
        }
    }
    return @out;
}

=head3 lookup

    my $res = $cdd->lookup($fid, $options);

Return domain data for a single feature, as C<< { $fid => { domain_hits => [...],
site_annotations => [...], structural_motifs => [...] } } >>.

=cut

sub lookup
{
    my($self, $fid, $options) = @_;

    my $res = $self->lookup_fids([$fid], $options);
    return $res;
}

=head3 lookup_fids

    my $res = $cdd->lookup_fids(\@fids, $options);

Bulk form of L</lookup>. One backend call covers every uncached sequence.

=cut

sub lookup_fids
{
    my($self, $fids, $options) = @_;

    $options = {} unless ref($options) eq 'HASH';

    my($seq_of, $md5_of) = $self->_sequences_for($fids);

    my $by_md5 = $self->_hits_for_md5s($fids, $seq_of, $md5_of, $options);

    my %out;
    for my $fid (@$fids)
    {
        my $md5 = $md5_of->{$fid};
        my $rows = (defined($md5) && $by_md5->{$md5}) ? $by_md5->{$md5} : [];
        $out{$fid} = {
            domain_hits      => _apply_mode($rows, $options->{data_mode}),
            site_annotations => [],
            structural_motifs=> [],
        };
    }
    return \%out;
}

=head3 prefetch

    my $status = $cdd->prefetch(\@fids, $options);

Warm the cache for a list of features in one backend call. Returns
C<< { ready => $n, pending => \@md5s, backend => $name } >>. A non-empty
C<pending> list means the NCBI job has been submitted and recorded but had not
finished within the wait budget; call again later to collect it.

=cut

sub prefetch
{
    my($self, $fids, $options) = @_;

    $options = {} unless ref($options) eq 'HASH';

    my($seq_of, $md5_of) = $self->_sequences_for($fids);
    my $by_md5 = $self->_hits_for_md5s($fids, $seq_of, $md5_of, $options);

    my %pending = map { $_ => 1 } @{$self->{_pending} || []};

    return {
        ready   => scalar(keys %$by_md5),
        pending => [sort keys %pending],
        backend => $self->{backend},
    };
}

#
# Resolve translations and md5s for a list of fids. We md5 the sequence
# ourselves rather than calling md5_of_peg, which falls back to a per-peg
# get_translation loop for pegs that are not already in its index.
#
sub _sequences_for
{
    my($self, $fids) = @_;

    my $fig = $self->fig;
    my(%seq, %md5);

    my $bulk;
    if ($fig->can('get_translation_bulk'))
    {
        $bulk = eval { $fig->get_translation_bulk($fids) };
    }

    for my $fid (@$fids)
    {
        my $s = $bulk ? $bulk->{$fid} : undef;
        $s = $fig->get_translation($fid) unless defined($s) && $s ne '';
        next unless defined($s) && $s ne '';
        $s =~ s/\s+//g;
        $seq{$fid} = $s;
        $md5{$fid} = md5_hex(uc($s));
    }
    return (\%seq, \%md5);
}

#
# The core: resolve a set of md5s to hit rows, consulting the in-process memo,
# then the on-disk cache, then the backend.
#
sub _hits_for_md5s
{
    my($self, $fids, $seq_of, $md5_of, $options) = @_;

    my %want;           # md5 => seq
    my %genome_of_md5;  # md5 => a genome that contains it (for cache placement)

    for my $fid (@$fids)
    {
        my $md5 = $md5_of->{$fid};
        next unless defined($md5);
        $want{$md5} = $seq_of->{$fid};
        $genome_of_md5{$md5} ||= FIG::genome_of($fid);
    }

    my %have;
    my @missing;

    for my $md5 (keys %want)
    {
        if (my $m = $self->{_mem}->{$md5})
        {
            $have{$md5} = $m;
            next;
        }
        my $rows = $self->_read_cache($genome_of_md5{$md5}, $md5);
        if (defined($rows))
        {
            $self->{_mem}->{$md5} = $rows;
            $have{$md5} = $rows;
            next;
        }
        push(@missing, $md5);
    }

    $self->{_pending} = [];

    if (@missing && !$options->{cached_only})
    {
        my $fetched;
        if ($self->{backend} eq 'rpsblast')
        {
            $fetched = $self->_run_rpsblast(\@missing, \%want, $options);
        }
        elsif ($self->{backend} eq 'jsonrpc')
        {
            $fetched = $self->_run_jsonrpc(\@missing, \%want, $options);
        }
        else
        {
            $fetched = $self->_run_ncbi(\@missing, \%want, $options);
        }

        for my $md5 (keys %$fetched)
        {
            my $rows = $fetched->{$md5};
            $self->{_mem}->{$md5} = $rows;
            $have{$md5} = $rows;
            $self->_write_cache($genome_of_md5{$md5}, $md5, $rows);
        }

        #
        # Anything we asked for and did not get back is still in flight.
        #
        my @still = grep { !exists $fetched->{$_} } @missing;
        $self->{_pending} = \@still;
    }

    return \%have;
}

#--------------------------------------------------------------------------
# Cache
#--------------------------------------------------------------------------

#
# organism_directory returns undef for a genome that is not present, and is
# overridden by FIGV/FIGM, so go through the method rather than building the
# path. When there is nowhere to put the cache we simply do not cache.
#
sub cache_dir
{
    my($self, $genome) = @_;
    return undef unless defined($genome) && $genome ne '';
    my $d = eval { $self->fig->organism_directory($genome) };
    return undef unless defined($d) && $d ne '' && -d $d;
    return "$d/CDD";
}

=head3 is_cached

    my $bool = $cdd->is_cached($genome, $md5);

True when this protein already has a cache entry, including the negative
entry written for a protein that genuinely has no domains. Lets a bulk
caller skip work without going through a lookup.

=cut

sub is_cached
{
    my($self, $genome, $md5) = @_;
    my $f = $self->_cache_file($genome, $md5);
    return (defined($f) && -f $f) ? 1 : 0;
}

sub _cache_file
{
    my($self, $genome, $md5) = @_;
    my $dir = $self->cache_dir($genome);
    return undef unless $dir;
    return "$dir/" . substr($md5, 0, 2) . "/$md5";
}

#
# Returns undef on a miss, or an arrayref of hit rows (possibly empty, for a
# protein we have looked up and which genuinely has no domains).
#
sub _read_cache
{
    my($self, $genome, $md5) = @_;

    my $f = $self->_cache_file($genome, $md5);
    return undef unless defined($f) && -f $f;

    open(my $fh, "<", $f) or return undef;
    my @rows;
    my $marked;
    while (<$fh>)
    {
        chomp;
        next if $_ eq '';
        if (/^#/)
        {
            $marked = 1 if /^#none/;
            next;
        }
        push(@rows, [split(/\t/, $_, -1)]);
    }
    close($fh);

    #
    # An empty file is ambiguous -- it could be a truncated write -- so a
    # domain-free protein is recorded with an explicit #none marker and an
    # otherwise-empty file is treated as a miss.
    #
    return undef if !@rows && !$marked;
    return \@rows;
}

sub _write_cache
{
    my($self, $genome, $md5, $rows) = @_;

    my $f = $self->_cache_file($genome, $md5);
    return unless defined($f);

    my $dir = $f;
    $dir =~ s,/[^/]+$,,;
    if (! -d $dir)
    {
        eval { make_path($dir) };
        return unless -d $dir;
    }

    #
    # Write-then-rename: NFS append is not atomic, and a half-written cache
    # file would be indistinguishable from a real result.
    #
    my $tmp = "$f~$$";
    open(my $fh, ">", $tmp) or return;
    if (@$rows)
    {
        print $fh join("\t", map { defined($_) ? $_ : '' } @$_), "\n" for @$rows;
    }
    else
    {
        print $fh "#none\n";
    }
    unless (close($fh))
    {
        unlink($tmp);
        return;
    }
    rename($tmp, $f) or unlink($tmp);
}

#--------------------------------------------------------------------------
# Local rpsblast backend
#--------------------------------------------------------------------------

sub _run_rpsblast
{
    my($self, $md5s, $seq_of, $options) = @_;

    my $tmp = File::Temp->new(TEMPLATE => "cddqueryXXXXXX", TMPDIR => 1, SUFFIX => ".fa");
    my $n = 0;
    for my $md5 (@$md5s)
    {
        my $s = $seq_of->{$md5};
        next unless defined($s) && $s ne '';
        print $tmp ">$md5\n$s\n";
        $n++;
    }
    close($tmp);
    return {} unless $n;

    my @cmd = ($self->{rps_prog},
               '-query',       "$tmp",
               '-db',          $self->{rps_db},
               '-evalue',      $self->{evalue},
               '-num_threads', $self->{rps_threads},
               '-outfmt',      '6 qseqid sseqid qstart qend evalue bitscore stitle');

    my %raw;
    my $pid = open(my $fh, "-|", @cmd);
    if (!$pid)
    {
        warn "ConservedDomainSearch: cannot run $self->{rps_prog}: $!\n";
        return {};
    }
    while (<$fh>)
    {
        chomp;
        my($qid, $sid, $qstart, $qend, $evalue, $bits, $stitle) = split(/\t/, $_, 7);
        next unless defined($stitle);

        my($pssmid) = $sid =~ /CDD\|(\d+)/;
        $pssmid = $sid unless defined($pssmid);

        #
        # rpsblast's subject title is "<accession>, <short name>, <description>".
        # Descriptions contain commas, so the split is limited to three fields.
        #
        my($acc, $short, $desc) = split(/, /, $stitle, 3);
        $acc   = '' unless defined($acc);
        $short = '' unless defined($short);
        $desc  = '' unless defined($desc);
        $desc =~ s/\s+$//;

        push(@{$raw{$qid}}, ['', $pssmid, $qstart, $qend, $evalue, $bits,
                             $acc, $short, '-', '', $desc]);
    }
    close($fh);

    #
    # Every md5 we submitted gets an entry, so a protein with no domains is
    # cached as a definite negative rather than re-queried forever.
    #
    my %out;
    for my $md5 (@$md5s)
    {
        my $s = $seq_of->{$md5};
        $out{$md5} = _classify_hits($raw{$md5} || [], defined($s) ? length($s) : 0);
    }
    return \%out;
}

#
# NCBI's "rep" data mode returns one representative hit per region rather than
# every profile that matched. rpsblast has no equivalent, and without it a
# single protein yields a dozen stacked hits that render as an unreadable
# smear. Approximate it by walking hits best-score-first and keeping one only
# if it does not substantially overlap a better one already kept.
#
# Classification into Specific/Non-specific proper needs NCBI's
# bitscore_specific.txt thresholds, which we do not ship; "best hit covering
# this region" is the useful approximation and is what gets drawn.
#
sub _classify_hits
{
    my($hits, $qlen) = @_;

    my @sorted = sort { $b->[COL_BITSCORE] <=> $a->[COL_BITSCORE] } @$hits;

    #
    # CDD mixes single-domain models (cd, pfam, smart) with whole-protein
    # cluster models (PRK, COG, TIGR, PLN). The latter match end to end and,
    # being the highest-scoring hit, would mask the architecture underneath
    # them -- E. coli thrA's best hit is PRK09436 across all 820 residues,
    # which as a drawn feature says nothing the gene arrow does not. Prefer
    # hits that cover part of the protein, and fall back to the full-length
    # ones only when a protein really is a single domain.
    #
    my @candidates = @sorted;
    if ($qlen && $qlen > 0)
    {
        my @partial = grep {
            my $len = abs($_->[COL_TO] - $_->[COL_FROM]) + 1;
            $len < 0.9 * $qlen;
        } @sorted;
        @candidates = @partial if @partial;
    }

    $_->[COL_TYPE] = 'Non-specific' for @sorted;

    my @kept;
    for my $h (@candidates)
    {
        my($s, $e) = ($h->[COL_FROM], $h->[COL_TO]);
        ($s, $e) = ($e, $s) if $s > $e;
        my $len = $e - $s + 1;
        next unless $len > 0;

        my $redundant = 0;
        for my $k (@kept)
        {
            my($ks, $ke) = ($k->[COL_FROM], $k->[COL_TO]);
            ($ks, $ke) = ($ke, $ks) if $ks > $ke;
            my $os = $s > $ks ? $s : $ks;
            my $oe = $e < $ke ? $e : $ke;
            my $ov = $oe - $os + 1;
            if ($ov > 0 && $ov > 0.5 * $len)
            {
                $redundant = 1;
                last;
            }
        }
        next if $redundant;
        $h->[COL_TYPE] = 'Specific';
        push(@kept, $h);
    }

    #
    # Every hit is returned, and cached: the rep/full choice is a property of
    # the caller, not of the data, so filtering before the cache write would
    # make data_mode => 'full' permanently unable to see what rep discarded.
    #
    return \@sorted;
}

sub _apply_mode
{
    my($rows, $mode) = @_;
    return $rows if ($mode || 'rep') eq 'full';
    return [grep { lc($_->[COL_TYPE]) eq 'specific' } @$rows];
}

#--------------------------------------------------------------------------
# NCBI Batch CD-Search backend
#--------------------------------------------------------------------------

#
# A pending registry lets a submission survive the request that made it: the
# cdsid is recorded before we ever poll, so a timeout or a crash cannot strand
# an NCBI job. It lives in FIG_Config::var rather than with the genome because
# one cdsid spans many md5s across many genomes and is transient -- NCBI keeps
# a job for about two days.
#
sub _pending_dir
{
    my $d = ($FIG_Config::var || '') . "/cdd_pending";
    return undef unless $FIG_Config::var;
    if (! -d $d)
    {
        eval { make_path($d) };
        return undef unless -d $d;
    }
    return $d;
}

sub _read_pending
{
    my($self, $md5) = @_;
    my $d = _pending_dir() or return undef;
    open(my $fh, "<", "$d/$md5") or return undef;
    my $l = <$fh>;
    close($fh);
    return undef unless defined($l);
    chomp $l;
    my($cdsid, $when) = split(/\t/, $l);
    return $cdsid;
}

sub _write_pending
{
    my($self, $md5s, $cdsid) = @_;
    my $d = _pending_dir() or return;
    my $now = time;
    for my $md5 (@$md5s)
    {
        my $tmp = "$d/$md5~$$";
        open(my $fh, ">", $tmp) or next;
        print $fh "$cdsid\t$now\n";
        close($fh) or do { unlink($tmp); next };
        rename($tmp, "$d/$md5") or unlink($tmp);
    }
}

sub _clear_pending
{
    my($self, $md5s) = @_;
    my $d = _pending_dir() or return;
    unlink("$d/$_") for @$md5s;
}

sub _run_ncbi
{
    my($self, $md5s, $seq_of, $options) = @_;

    my %out;
    my @todo = @$md5s;

    #
    # Anything already submitted by an earlier request gets one cheap poll
    # before we consider submitting it again.
    #
    my %by_cdsid;
    for my $md5 (@todo)
    {
        my $cdsid = $self->_read_pending($md5) or next;
        push(@{$by_cdsid{$cdsid}}, $md5);
    }
    for my $cdsid (keys %by_cdsid)
    {
        my $res = $self->_ncbi_poll($cdsid);
        next unless $res->{done};
        my $hits = $res->{hits};
        for my $md5 (@{$by_cdsid{$cdsid}})
        {
            my $s = $seq_of->{$md5};
            $out{$md5} = _classify_hits($hits->{$md5} || [], defined($s) ? length($s) : 0);
        }
        $self->_clear_pending($by_cdsid{$cdsid});
    }

    @todo = grep { !exists $out{$_} && !$self->_read_pending($_) } @todo;
    return \%out unless @todo;

    my $fasta = '';
    for my $md5 (@todo)
    {
        my $s = $seq_of->{$md5};
        next unless defined($s) && $s ne '';
        $fasta .= ">$md5\n$s\n";
    }
    return \%out unless $fasta ne '';

    my $res = $self->ua->post($NCBI_URL, [
        queries  => $fasta,
        db       => 'cdd',
        smode    => 'auto',
        useid1   => 'true',
        evalue   => $self->{evalue},
        maxhit   => 500,
        dmode    => 'rep',
        tdata    => 'hits',
    ]);
    if (!$res->is_success)
    {
        warn "ConservedDomainSearch: NCBI submit failed: " . $res->status_line . "\n";
        return \%out;
    }
    my($cdsid) = $res->content =~ /^#cdsid\s+(\S+)/m;
    if (!$cdsid)
    {
        warn "ConservedDomainSearch: no cdsid in NCBI response\n";
        return \%out;
    }

    #
    # Record before polling: if we die or time out from here on, the job is
    # still recoverable by the next request.
    #
    $self->_write_pending(\@todo, $cdsid);

    my $budget = defined($FIG_Config::cdd_wait_budget) ? $FIG_Config::cdd_wait_budget : 10;
    my $deadline = time + $budget;
    while (time < $deadline)
    {
        sleep 2;
        my $p = $self->_ncbi_poll($cdsid);
        next unless $p->{done};
        for my $md5 (@todo)
        {
            my $s = $seq_of->{$md5};
            $out{$md5} = _classify_hits($p->{hits}->{$md5} || [], defined($s) ? length($s) : 0);
        }
        $self->_clear_pending(\@todo);
        last;
    }

    return \%out;
}

#
# Returns { done => bool, hits => { md5 => [rows] } }.
#
sub _ncbi_poll
{
    my($self, $cdsid) = @_;

    my $res = $self->ua->post($NCBI_URL, [ cdsid => $cdsid, tdata => 'hits',
                                           dmode => 'rep', cddefl => 'true' ]);
    return { done => 0 } unless $res->is_success;

    my $body = $res->content;

    #
    # A completed response carries two #status lines (0, then "success"); a
    # running one carries a single "#status 3".
    #
    my($status) = $body =~ /^#status\s+(\S+)/m;
    return { done => 0 } if !defined($status) || $status eq '3';
    if ($status ne '0' && $status !~ /^success$/i)
    {
        warn "ConservedDomainSearch: NCBI job $cdsid failed with status $status\n";
        return { done => 1, hits => {} };
    }

    my %hits;
    for my $line (split(/\n/, $body))
    {
        next if $line =~ /^\s*$/;
        next if $line =~ /^#/;
        next if $line =~ /^Query\b/;
        my @f = split(/\t/, $line, -1);
        next unless @f >= 10;

        #
        # Column 0 is "Q#<n> - ><our fasta id>"; recover the md5 from it and
        # drop the column so the rest lines up with the documented layout.
        #
        my $q = shift @f;
        my($md5) = $q =~ /^Q#\d+\s*-\s*>?\s*(\S+)/;
        next unless defined($md5);

        push(@f, '') while @f < 11;      # no description column from NCBI
        push(@{$hits{$md5}}, \@f);
    }
    return { done => 1, hits => \%hits };
}

#--------------------------------------------------------------------------
# Legacy JSON-RPC backend
#--------------------------------------------------------------------------

sub _run_jsonrpc
{
    my($self, $md5s, $seq_of, $options) = @_;

    my @seqs = map { [$_, $_, $seq_of->{$_}] } grep { defined($seq_of->{$_}) } @$md5s;
    return {} unless @seqs;

    my $res = $self->lookup_seqs(\@seqs, $options);
    return {} unless ref($res) eq 'HASH';

    my %out;
    for my $md5 (@$md5s)
    {
        my $ent = $res->{$md5} or next;
        $out{$md5} = $ent->{domain_hits} || [];
    }
    return \%out;
}

sub lookup_seqs
{
    my($self, $seqs, $options) = @_;

    $options = {} unless ref($options) eq 'HASH';

    my $req = {
        id => $idx++,
        method => 'ConservedDomainSearch.cdd_lookup',
        params => [$seqs, $options],
    };
    my $res = $self->ua->post($self->url, Content => encode_json($req));
    if (!$res->is_success)
    {
        die "Failure invoking cdd_lookup: " . $res->status_line . "\n" .  $res->content;
    }
    my $data = decode_json($res->content);
    return $data->{result}->[0];
}

=head3 domains_of

    my $fidHash = $cdd->domains_of(\@fids);

Compute the conserved domains for a list of features. For each feature, this
method will return a list of the accessions of the specific conserved domains
found.

=cut

sub domains_of {
    my ($self, $fids) = @_;

    my %opts = (data_mode => 'rep');
    my $results = $self->lookup_fids($fids, \%opts);

    my %retVal;
    for my $fid (@$fids) {
        my @doms;
        my $hits = $results->{$fid}{domain_hits};
        for my $hit (@$hits) {
            #
            # The service emits this lowercase; the original compared against
            # 'Specific' and so never matched anything.
            #
            if (lc($hit->[COL_TYPE]) eq 'specific') {
                push @doms, $hit->[COL_ACCESSION];
            }
        }
        $retVal{$fid} = \@doms;
    }
    return \%retVal;
}

1;
