use Data::Dumper;

use SeedSearch;

use File::Slurp;
use LWP::UserAgent;
use strict;
use Encode;
use FIG_Config;
use FIG;
my $fig = new FIG;

@ARGV == 1 or die "Usage: sphinx_index_genome function|attrib|subsystem\n";

my $which = shift;

my @fields;
my %attr_types;
if ($which eq 'function')
{
    @fields = qw(annotation fid genome);
    $attr_types{$_} = "attr=\"string\"" for @fields;
}
elsif ($which eq 'genome')
{
    @fields = qw(genome name taxonomy contigs);
    $attr_types{$_} = "attr=\"string\"" for qw(genome name);
}
elsif ($which eq 'attrib')
{
    @fields = qw(genome alias subsystem annotation fid);
    $attr_types{$_} = 'attr="string"' foreach qw(annotation);
}
elsif ($which eq 'subsystem')
{
    @fields = qw(subsystem curator version classification role);
    $attr_types{$_} = "attr=\"string\"" for @fields;
}
elsif ($which eq 'dlit')
{
    @fields = qw(title abstract protein genome fid pmid pmcid);
    $attr_types{$_} = "attr=\"string\"" foreach (qw(title fid pmid pmcid));
}
else
{
    die "Unknown type $which\n";
}

print <<END;
<?xml version="1.0" encoding="utf-8"?>
<sphinx:docset>
<sphinx:schema>
END
print "<sphinx:field name=\"$_\" $attr_types{$_}/>\n" for @fields;
print <<END;
</sphinx:schema>
END

my %tmap = (peg => 1, rna => 2);

my @genomes;
if (my $glist = $ENV{SPHINX_INDEX_ONLY})
{
    @genomes = split(/,/, $glist);
}
else
{
    @genomes = $fig->genomes();
}

if ($which eq 'subsystem')
{
    &index_subsystems;
    print "</sphinx:docset>\n";
    exit;
}
elsif ($which eq 'genome')
{
    &index_genomes(\@genomes);
    print "</sphinx:docset>\n";
    exit;
}
elsif ($which eq 'dlit')
{
    &index_dlits();
    print "</sphinx:docset>\n";
    exit;
}

#
# Ingest the subsystem index.
#

my $next_id = 1;
    
for my $genome (@genomes)
{
	# print STDERR "$genome\n";
    my $gs = escape($fig->genus_species($genome));
    if ($which eq 'function')
    {
	my @fids = $fig->all_features($genome);
	my $fns = $fig->function_of_bulk(\@fids);

	for my $fid (keys %$fns)
	{
	    my $fn = escape($fns->{$fid});
	    my $docid = SeedSearch::fid_to_docid($fid);
	    if (!defined($docid))
	    {
		    #print STDERR "Skipping $fid docid not defined\n";
		next;
	    }
	    print <<END;
<sphinx:document id="$docid">
<genome>$gs</genome>
<annotation>$fn</annotation>
</sphinx:document>
END
	}
    }
    else
    {
	my %ss_info;
	
	my $sth = $fig->db_handle->{_dbh}->prepare(qq(SELECT i.protein, i.subsystem
						      FROM subsystem_index i LEFT JOIN aux_roles a ON i.role = a.role
						      WHERE i.protein LIKE 'fig|$genome.peg.%' AND a.subsystem IS NULL),
					       { mysql_use_result => 1});
	$sth->execute();
	while (my $ent = $sth->fetchrow_arrayref())
	{
	    my($prot, $ss) = @$ent;
	    $ss_info{$prot}->{$ss} = 1;
	}
	
	my $gs = escape($fig->genus_species($genome));
	
	my $all_data = $fig->all_features_detailed_fast($genome);
	
	my $ext_aliases_l = $fig->db_handle->SQL(qq(SELECT id, alias
						    FROM ext_alias
						    WHERE id like 'fig|${genome}.%'));
	my %ext_aliases;
	map { $ext_aliases{$_->[0]}->{$_->[1]}++ } @$ext_aliases_l;
	
	for my $feature (@$all_data)
	{
	    my($fid, $loc, $aliases, $type, $b, $e, $func, $who) = @$feature;
	    
	    # my @ss = $fig->peg_to_subsystems($fid, 1, 1);
	    my @ss = keys %{$ss_info{$fid}};
	    @ss = map { defined($_) ? encode_utf8($_) : () } @ss;
	    my $ss = escape(join("\n", map { s/_/ /g; $_ } @ss));
	    
	    $func = defined($func) ? escape($func) : "";
	    my $efid = escape($fid);
	    
	    my %aliases = map { $_ => 1 } split(",", $aliases);
	    map { $aliases{$_} = 1 } keys %{$ext_aliases{$fid}};
	    my @aliases = keys %aliases;
	    my $alias_txt = "";
	    if (@aliases)
	    {
		$alias_txt = escape(join("\n",
					 map { s/&/&amp;/g;
					       s/</&lt;/g;
					       s/>/&gt;/g;
					       $_ } @aliases));
	    }
	    my $docid = SeedSearch::fid_to_docid($fid);
	    if (!$docid)
	    {
		    #warn "Skipping $fid docid not defined\n";
		next;
	    }
	    print <<END;
<sphinx:document id="$docid">
<genome>$genome $gs</genome>
<fid>$efid</fid>
<annotation>$func</annotation>
<alias>$alias_txt</alias>
<subsystem>$ss</subsystem>
</sphinx:document>
END
	}
    }
}
print "</sphinx:docset>\n";

sub index_genomes
{
    my($genomes) = @_;
    my $i = 1;
    for my $g (@$genomes)
    {
	my $gs = escape($fig->genus_species($g));
	my $tax =escape($fig->taxonomy_of($g));

	my $contigs = escape(join(" ", $fig->all_contigs($g)));

	print "<sphinx:document id=\"$i\">\n";

	print "<genome>$g</genome>\n";
	print "<name>$gs</name>\n";
	print "<contigs>$contigs</contigs>\n";
	print "<taxonomy>$tax</taxonomy>\n";
	print "</sphinx:document>\n";
	$i++;
    }
}

sub index_dlits
{

    my @dlits = $fig->get_attributes(undef, 'evidence_code', 'dlit%');

    my $ua = LWP::UserAgent->new;
    my $url_fmt = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=%s&retmode=text&rettype=abstract';
    # @fields = qw(title abstract protein fid);
    my $cache = "$FIG_Config::var/pm_abstract_cache";
    -d $cache || mkdir $cache;
    my $id = 1;
    my %by_pmid;
    for my $dlit (@dlits)
    {
	my($pmid) = $dlit->[2] =~ /dlit\((\d+)\)/;
	next unless $pmid && $pmid > 10;
	push(@{$by_pmid{$pmid}}, $dlit);
    }
    for my $pmid (keys %by_pmid)
    {
	last if ($ENV{LIMIT_DLIT_INDEX} && $id > $ENV{LIMIT_DLIT_INDEX});

	my $abs;
	my $cache_file = "$cache/$pmid";
	if (-s $cache_file)
	{
	    $abs = read_file($cache_file);
	}
	else
	{
	    my $url = sprintf($url_fmt, $pmid);
	    my $resp = $ua->get($url);
	    if ($resp->is_success)
	    {
		$abs = $resp->content;
		open(S, ">", $cache_file);
		print S $abs;
		close(S);
	    }
	}
	my $eabs = escape($abs);
	my $details = Dlits::get_pubmed_document_details($fig, $pmid);
	next unless ref($details);

	my %pegs;
	my %md5;

	for my $ent (@{$by_pmid{$pmid}})
	{
	    if ($ent->[0] =~ /Protein:(\S+)/)
	    {
		my $md5 = $1;
		my @pegs = $fig->pegs_with_md5($md5);
		$md5{$md5} = 1;
		$pegs{$_} = 1 foreach @pegs;
	    }
	    elsif ($ent->[0] =~ /^fig\|/)
	    {
		$pegs{$ent->[0]} = 1;
	    }
	}
	    
	my($pmcid) = $abs =~ /PMCID:\s+(\S+)/mg;

	my @pegs = sort { &FIG::by_fig_id($a, $b) } keys %pegs;

	my %genomes = map { &FIG::genome_of($_) } @pegs;
	my @genomes = sort { &FIG::by_genome_id($a, $b) } keys %genomes;
	my @genome_names = map { $fig->genus_species($_) } @genomes;

	my $genomes = escape(join(" ", @genomes));
	my $genome_names = escape(join(" ", @genome_names));
	my $fids = escape(join(" ", @pegs));
	my $md5 = escape(join(" ", keys %md5));

	print "<sphinx:document id=\"$id\">\n";
	print "<title>" . escape($details->{title}) . "</title>\n";
	print "<abstract>$eabs</abstract>\n";
	print "<protein>$md5</protein>\n";
	print "<fid>$fids</fid>\n";
	print "<genome>$genomes $genome_names</genome>\n";
	print "<pmid>$pmid</pmid>\n";
	print "<pmcid>$pmcid</pmcid>\n";
	print "</sphinx:document>\n";
	$id++;
    }
}

sub index_subsystems
{
    my $i = 1;
    my @ss = $fig->all_subsystems_detailed();

    #
    # Pull subsystem roles from database.
    #
    my %ss_to_role;
    my $res = $fig->db_handle->SQL("SELECT subsystem, role FROM subsystem_nonaux_role");
    for my $ent (@$res)
    {
	my($ss, $role) = @$ent;
	$ss_to_role{$ss} .= " " . $role;
    }
    undef $res;

    for my $ent (@ss)
    {
	print "<sphinx:document id=\"$i\">\n";
	my $roles = $ss_to_role{$ent->{subsystem}};
	$ent->{subsystem} =~ s/_/ /g;
	for my $f (grep { $_ ne 'role' } @fields)
	{
	    print "<$f>" . escape($ent->{$f}) . "</$f>\n";
	}
	print "<role>" . escape($roles) . "</role>\n";
	print "</sphinx:document>\n";
	$i++;
    }
}

sub escape
{
    my($s) = @_;
    return "" unless defined($s);
    $s = encode_utf8($s);
    $s =~ s/&/&amp;/g;
    $s =~ s/</&lt;/g;
    $s =~ s/>/&gt;/g;
    return $s;
}
