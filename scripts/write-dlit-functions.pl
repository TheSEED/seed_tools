#
# Write the set of SEED functions that are covered by dlits
#

use strict;
use FIG;
use SeedUtils;
use Data::Dumper;

my $fig = new FIG;

my %md5;

for my $dlit ($fig->all_dlits)
{
    $md5{$dlit->[1]} = 1;
}

my @md5s = sort keys %md5;

my $batch_size = 100;

my $rdbH = $fig->db_handle;

my %funcs;
while (@md5s)
{
    my @batch = splice(@md5s, 0, $batch_size);
    my $q = join(", ", map { "?" } @batch);

    my $res = $rdbH->SQL("SELECT assigned_function
 		FROM protein_sequence_MD5 m JOIN assigned_functions f ON m.id = f.prot
		WHERE m.md5 in ($q)", undef, @batch);
    $funcs{SeedUtils::strip_func($_->[0])} = 1 foreach @$res;
}

for my $func (sort { $a cmp $b } keys %funcs)
{
    print "$func\n";
}
