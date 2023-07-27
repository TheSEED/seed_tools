use List::Util 'first';
use Data::Dumper;
use strict;
use SeedUtils;
use FIGV;
use FIG;
use Getopt::Long;

#
# Write a SEED organism as a tab-delimited file.
#

my $org_dir;
my $out_file;
my $url;
my $rc = GetOptions("url=s" => \$url,
		    "orgdir=s" => \$org_dir,
		    "outfile=s" => \$out_file);

($rc && @ARGV == 1) or die "Usage: seed2txt [--orgdir genome-dir] [--outfile file] [--url url] genome-id \n";

my $genome = shift;

my $out_fh;
if ($out_file ne '')
{
    open($out_fh, ">", $out_file) or die "Cannot write $out_file: $!";
}
else
{
    open($out_fh, ">&", \*STDOUT);
}

my $fig;
if ($org_dir)
{
    $fig = FIGV->new($org_dir);
}
else
{
    $fig = FIG->new();
}

my $all_features = $fig->all_features_detailed_fast($genome);

push(@$_, boundaries_of($_->[1])) foreach @$all_features;
my $dir = $fig->organism_directory($genome);
my %figfam;
if (open(F, "<", "$dir/found"))
{
    while (<F>)
    {
	chomp;
	my($id, $ff) = split(/\t/);
	$figfam{$id} = $ff;
    }
    close(F);
}

my %ev_code;
if (open(F, "<", "$dir/evidence.codes"))
{
    while (<F>)
    {
	chomp;
	my($id, $ev) = split(/\t/);
	push(@{$ev_code{$id}}, $ev);
    }
    close(F);
}

if ($url)
{
    print $out_fh join("\t", qw(contig_id feature_id link type location start stop strand
				function aliases figfam evidence_codes nucleotide_sequence aa_sequence)), "\n";
}
else
{
    print $out_fh join("\t", qw(contig_id feature_id type location start stop strand
				function aliases figfam evidence_codes nucleotide_sequence aa_sequence)), "\n";
}

for my $f (sort { $a->[9] cmp $b->[9] or $a->[10] <=> $b->[10] } @$all_features)
{
    my($fid, $loc, $aliases, $type, $first, $last, $func, $who, undef, $contig, $left, $right, $strand) = @$f;

    my $trans = '';
    $trans = $fig->get_translation($fid) if $type eq 'peg';
    my $dna = $fig->dna_seq($genome, $loc);
    my $ev = $ev_code{$fid};
    $ev = ref($ev) ? join(" ", @$ev) : "";

    my @link = ();
    if ($url)
    {
	my $link = $url;
	$link =~ s/PEG/$fid/g;
	$link = "<a href='$link'>$link</a>";
	@link = ($link);
    }
    print $out_fh join("\t", $contig, $fid, @link, $type, $loc,
		       ($strand eq '+' ? ($first, $last) : ($last, $first)), $strand,
		       $func, $aliases,
		       $figfam{$fid},
		       $ev,
		       $dna, $trans), "\n";
}
