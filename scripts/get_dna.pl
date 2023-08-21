# -*- perl -*-
########################################################################
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
########################################################################

use FIG;
use Carp;
use Data::Dumper;

# usage: get_dna Contigs Tbl

(($contigs = shift @ARGV) && (-s $contigs))
    || die "\n   usage: get_dna Contigs Tbl\n\n";

open(CONTIGS,"<$contigs") || die "Could not read-open contigs file $contigs";

if (@ARGV && (-s $ARGV[0]) && ($tbl = shift @ARGV))
{
    open($TBL, "<$tbl") || die "Could not read-open tbl file $tbl";
}
else
{
    $TBL = \*STDIN;
}

$trouble = 0;
while (($id, $seqP) = &FIG::read_fasta_record(\*CONTIGS))
{
    if (not defined($contigs{$id})) {
	$contigs{$id} = $$seqP;
    }
    else {
	$trouble = 1;
	warn "Duplicate contig ID: $id\n";
    }
}
close(CONTIGS);

if ($trouble) {
    die "\nAborting due to duplicate contig IDs\n\n";
}


while (defined($entry = <$TBL>)) {
    if ($entry =~ /^(\S+)\t(\S+)/) {
	($fid, $loc) = ($1, $2);
	
	if ($seq = lc(&FIG::extract_seq(\%contigs, $loc))) {
	    &FIG::display_id_and_seq($fid, \$seq);
	}
	else {
	    warn "Could not extract sequence for entry: $entry";
	}
    }
    else {
	die "Malformed tbl entry: $entry";
    }
}
close($TBL);
