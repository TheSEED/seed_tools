#
# Copyright (c) 2003-2006 University of Chicago and Fellowship
# for Interpretations of Genomes. All Rights Reserved.
#
# This file is part of the SEED Toolkit.
# 
# The SEED Toolkit is free software. You can redistribute
# it and/or modify it under the terms of the SEED Toolkit
# Public License. 
#
# You should have received a copy of the SEED Toolkit Public License
# along with this program; if not write to the University of Chicago
# at info@ci.uchicago.edu or the Fellowship for Interpretation of
# Genomes at veronika@thefig.info or download a copy from
# http://www.theseed.org/LICENSE.TXT.
#

#
# Rename the genome from whatever it is in the given SEED skeleton dir to
# the genome ID passed on the command line. Used to regenerate the genome dir
# after registering the taxon id with the clearinghouse.
#
# Allow the explicit identification of old-id for the purposes of RAST upload
# where we renumber out of a parse_genbank output that just has a tax id in the
# identifiers.
#

use strict;
use File::Copy;
use Getopt::Long;


my $usage = "Usage: renumber_seed_dir [--old-id ID] [--exists-ok] old-dir new-genome-id new-dir";

my($exists_ok, $old_id);
my $rc = GetOptions("exists-ok" => \$exists_ok,
		    "old-id=s" => \$old_id);

($rc && @ARGV == 3) or die $usage;

my ($old_dir, $new_id, $new_dir) = @ARGV;

(-d $old_dir) or die "Old dir $old_dir not found\n";
$old_dir =~ s/\/$//o;

if (!$old_id)
{
    ($old_id) = ($old_dir =~ m/(\d+\.\d+)/o);
}

my $old_id_regexp = qr/$old_id/;

if (-d $new_dir && !$exists_ok)
{
    die "New directory $new_dir must not exist\n";
}

$new_id =~ /^\d+\.\d+/ or die "Invalid genome id $new_id\n";

#
# Make directories first.
#
(-d $new_dir) or mkdir($new_dir) or die "Cannot mkdir $new_dir $!";

my @feature_types;
if (opendir(D, "$old_dir/Features"))
{
    mkdir "$new_dir/Features" or die "mkdir $new_dir/Features failed: $!";
    
    @feature_types = grep { $_ !~ /^\./ and -d "$old_dir/Features/$_" } readdir(D);
    closedir(D);
    
    for my $ft (@feature_types)
    {
	mkdir("$new_dir/Features/$ft") or die "Cannot mkdir $new_dir/Features/$ft: $!";
    }
}

#
# Copy the plain files over.
#

opendir(D, $old_dir) or die "Cannot open dir $old_dir: $!";

my @top_files = grep { $_ !~ /^\./ and -f "$old_dir/$_" } readdir(D);
closedir(D);

for my $file (@top_files)
{
    #
    # original code only did this for some; is there any reason we should not replace all?
    #
    copy_and_replace("$old_dir/$file", "$new_dir/$file", $old_id_regexp, $new_id)
	or die "copy $file failed: $!";
}

#
# Handle the features.
#

for my $ft (@feature_types)
{
    my $ofd = "$old_dir/Features/$ft";
    my $nfd = "$new_dir/Features/$ft";

    copy_and_replace("$ofd/fasta", "$nfd/fasta", $old_id_regexp, $new_id);
    copy_and_replace("$ofd/tbl", "$nfd/tbl", $old_id_regexp, $new_id);
}

sub copy_and_replace
{
    my ($old, $new, $old_id_regexp, $new_id) = @_;

    open(OLD, "<$old") or die "Cannot open $old: $!";
    open(NEW, ">$new") or die "Cannot open $new: $!";

    while (<OLD>)
    {
	s/fig\|$old_id_regexp\.(\w+)/fig|$new_id.$1/g;
	print NEW $_;
    }
    close(OLD);
    close(NEW);
}
