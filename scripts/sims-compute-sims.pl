#
# Compute sims with diamond
#

use strict;
use Getopt::Long::Descriptive;
use File::Path qw(make_path);
use IPC::Run 'run';

my($opt, $usage) = describe_options("%c %o sims-dir",
				    ["threads=i", "Number of threads to use", { default => 4 }],
				    ["keep=i", "Keep this many matches", { default => 200 }],
				    ["diamond-sensitivity=s", "Diamond sensitivity flag to use", { default => "sensitive" }],
				    ["diamond-block-size=s", "Diamond block size", { default => 5 }],
				    ["tmpdir|t=s", "Diamond temp dir"],
				    ["help|h", "Show this help message"]);

print($usage->text), exit 0 if $opt->help;
die($usage->text) if @ARGV != 1;

my $sims_dir = shift;

my($tag) = $sims_dir =~ /(\d+)$/;
make_path("$sims_dir/sims.split", "$sims_dir/sims.raw");

#
# We run diamond on each of the files in the seqs directory, generating a gz file in the sims directory.
# Then these are gzcatted into the reduce/split pipeline.
#

my $seq_dir = "$sims_dir/seqs";
my $raw_dir = "$sims_dir/sims.raw";

opendir(D, $seq_dir) or die "Cannot opendir $seq_dir: $!";
my @seq_files = sort grep { -f "$seq_dir/$_" } readdir(D);
closedir(D);

my @out_files;

for my $seq_file (@seq_files)
{
    my $sim_file = $seq_file;
    $sim_file =~ s/seq/sim/;
    my $sim_path = "$raw_dir/$sim_file";

    my @cmd = ("diamond", "blastp", "-v",
		  "--block-size", $opt->diamond_block_size,
		  "-c1",
		  "-e", "1e-4",
		  "--masking", 0,
		  "--compress", 1,
		  "--threads", $opt->threads,
		  "--" . $opt->diamond_sensitivity,
		  "-k", $opt->keep,
		  ($opt->tmpdir ? ("--tmpdir", $opt->tmpdir) : ()),
		  "--outfmt", qw(6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen),
		  "-d", "$sims_dir/nr",
		  "-q", "$seq_dir/$seq_file",
		  "-o", $sim_path);
print "@cmd\n";
    my $rc = system(@cmd);

    $rc == 0 or die "Error running @cmd\n";

    push(@out_files, "$sim_path.gz");
}

my $ok = run(["zcat", @out_files],
	  "|", 
	  ["reduce_sims", "$sims_dir/nr-len.btree", $opt->keep],
	  "|",
	  ["split_sims", "$sims_dir/sims.split", "sims.$tag"]);
$ok or die "Reduce/split pipeline failed\n";
