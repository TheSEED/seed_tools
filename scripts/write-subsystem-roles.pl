use FIG;
use strict;
my $fig = new FIG;

@ARGV == 1 or die "Usage: $0 role-file\n";
my $out = shift;

my $list = $fig->subsystem_roles();
open(D, ">", $out) or die "Cannot write $out: $!";
print D "$_\n" foreach sort { $a cmp $b } keys %$list;
close(D);
