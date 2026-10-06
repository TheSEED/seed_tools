#
# Warm the conserved-domain cache for one or more genomes, so that the
# compare-regions CDD track draws from cache instead of waiting on a backend.
#
# Against the local rpsblast backend a bacterial genome takes a few minutes;
# against NCBI it is measured in hours and will mostly leave work pending, so
# this is really a tool for the local backend.
#
# Usage:
#   cdd_prefill [options] genome-id [genome-id ...]
#   cdd_prefill [options] -a
#
# Options:
#   -a            every genome in the SEED, smallest first
#   -b N          proteins per backend call (default 200; NCBI's own limit
#                 is 1000 sequences per submission)
#   -f            re-fetch and overwrite entries that are already cached
#   -n            say what would be done, fetch nothing
#   -v            report each batch as it completes
#

use strict;
use warnings;

use FIG;
use FIG_Config;
use ConservedDomainSearch;
use Getopt::Long;
use Digest::MD5 'md5_hex';

my($all, $batch_size, $force, $dry_run, $verbose, $help) = (0, 200, 0, 0, 0, 0);

GetOptions("a"   => \$all,
           "b=i" => \$batch_size,
           "f"   => \$force,
           "n"   => \$dry_run,
           "v"   => \$verbose,
           "h|help" => \$help) or die usage();
die usage() if $help;

my $fig = FIG->new;

my @genomes = @ARGV;
if ($all)
{
    #
    # Smallest first: a run that is cut short has then finished the most
    # genomes rather than the fewest.
    #
    @genomes = sort { $fig->genome_szdna($a) <=> $fig->genome_szdna($b) } $fig->genomes('complete');
}
die usage() unless @genomes;

my $cds = ConservedDomainSearch->new($fig);
print "backend: $cds->{backend}\n";
if ($cds->{backend} ne 'rpsblast')
{
    print STDERR "warning: the $cds->{backend} backend is slow enough that a full genome\n" .
                 "         will mostly leave work pending. Set \$FIG_Config::cdd_rps_db.\n";
}

my $t_start = time;
my($g_done, $total_skipped, $total_pending) = (0, 0, 0);

for my $genome (@genomes)
{
    my @pegs = $fig->pegs_of($genome);
    unless (@pegs)
    {
        print "$genome: no pegs, skipping\n";
        next;
    }

    #
    # Filter to what is actually missing before doing any work, so a
    # re-run over a warm genome costs a directory walk rather than a
    # backend call.
    #
    my @todo;
    if ($force)
    {
        @todo = @pegs;
    }
    else
    {
        for my $peg (@pegs)
        {
            my $seq = $fig->get_translation($peg);
            next unless defined($seq) && $seq ne '';
            $seq =~ s/\s+//g;
            my $md5 = md5_hex(uc($seq));
            push(@todo, $peg) unless $cds->is_cached($genome, $md5);
        }
    }

    my $have = @pegs - @todo;
    $total_skipped += $have;

    printf "%-14s %5d pegs, %5d cached, %5d to fetch\n",
           $genome, scalar(@pegs), $have, scalar(@todo);

    next if $dry_run || !@todo;

    my $t0 = time;
    my $pending = 0;
    while (@todo)
    {
        my @batch = splice(@todo, 0, $batch_size);
        my $st = eval { $cds->prefetch(\@batch, { data_mode => 'rep' }) };
        if ($@)
        {
            warn "$genome: batch failed: $@";
            next;
        }
        $pending += scalar @{$st->{pending} || []};
        printf "  batch of %d: ready=%d pending=%d\n",
               scalar(@batch), $st->{ready}, scalar(@{$st->{pending} || []})
            if $verbose;

        #
        # The memo exists to make one page render cheap; across a whole
        # genome it would grow without bound.
        #
        $cds->{_mem} = {};
    }

    $total_pending += $pending;
    $g_done++;
    printf "%-14s done in %ds%s\n", $genome, time - $t0,
           ($pending ? ", $pending still pending" : "");
}

printf "\n%d genome%s, %d already cached, %d pending, %ds total\n",
       $g_done, ($g_done == 1 ? "" : "s"), $total_skipped, $total_pending, time - $t_start;

sub usage
{
    return <<'END';
Usage: cdd_prefill [-b N] [-f] [-n] [-v] genome-id [genome-id ...]
       cdd_prefill [-b N] [-f] [-n] [-v] -a

  -a   all complete genomes, smallest first
  -b N proteins per backend call (default 200)
  -f   re-fetch entries that are already cached
  -n   dry run: report what is missing, fetch nothing
  -v   report each batch
END
}
