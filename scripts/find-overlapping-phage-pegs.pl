use Data::Dumper;
use strict;
use FIG;

my $fig = FIG->new;

for my $g ($fig->genomes)
{
    my $tax = $fig->taxonomy_of($g);
    next unless $tax =~ /^Viruses/;
    my $feats = $fig->all_features_detailed_fast($g);

    my %by_stop;
    for my $feat (@$feats)
    {
	my($fid, $loc, $alias, $type, $left, $right, $fn) = @$feat;
	next unless $type eq 'peg';
	if ($fn =~ /pp1ab|pp1a|pp1b/i)
	{
	    print "$fid\t$fn\n";
	}
	my($ctg, $start, $stop) = $loc =~ /^(.*?)_(\d+)_(\d+)$/;
	push(@{$by_stop{$ctg, $stop}}, [$fid, $loc, $fn, $right - $left]);
    }
    my @mult = grep { @{$by_stop{$_}} > 1 } keys %by_stop;
    for my $k (@mult)
    {
	my @cand = sort { $b->[3] <=> $a->[3] } @{$by_stop{$k}};
	my $longest = shift @cand;
	#
	# we print the ones we want to ignore
	#
	print "$_->[0]\t$_->[2]\n" foreach @cand;
    }
}
