#
# Finish the sims update.
#

use strict;
use Getopt::Long::Descriptive;
use IPC::Run qw(run);
use File::Path qw(make_path);
use File::Basename;

my($opt, $usage) = describe_options("%c %o sims-dir",
				    ["keep=i", "Keep this many matches", { default => 200 }],
				    ["help|h", "Show this help message"]);

print($usage->text), exit 0 if $opt->help;
die($usage->text) if @ARGV != 1;

my $sims_dir = shift;

#
# If there are previous sims, we need to flip and update.
#

my($tag) = $sims_dir =~ /(\d+)$/;
make_path("$sims_dir/Sims.$tag");


if (-d "$sims_dir/prev_sims" && <$sims_dir/prev_sims/*>)
{
    my $ok = run(["flip_sims", "$sims_dir/sims.split", "$sims_dir/sims.flipped"]);
    $ok or die "flip_sims failed\n";

    $ok = run(["update_sims2", 
	       "$sims_dir/nr-len.btree", 
	       "$sims_dir/peg.synonyms.reduce_sims_index.btree", 
	       $opt->keep,
	       "$sims_dir/prev_sims",
	       "$sims_dir/Sims.$tag",
	       "$sims_dir/sims.flipped",
	       "$sims_dir/ids.deleted"]);
    $ok or die "update_sims2 failed\n";
}

for my $split (<$sims_dir/sims.split/*>)
{
    my $f = basename($split);
    symlink($split, "$sims_dir/Sims.$tag/$f") or die "symlink $split $sims_dir/Sims.$tag/$f failed: $!";
}
    # my $ok = run(["cp", "-v", <$sims_dir/sims.split/*>, "$sims_dir/Sims.$tag"]);
    # $ok or die "Error copying sims into place\n";
