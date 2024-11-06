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

# Usage: count_bases < something.fasta or
#        cat something.fasta | count_bases

$dictP = base_dict_init();

$/ = "\n>";
while (defined($_ = <STDIN>))
{
    chomp;
    if ($_ =~ /^>?(\S+)[^\n]*\n(.*)/s)
    {
	$sid  =  $1;
	$seq =  $2;
	$seq =~ s/\n//gs;
	$seq =~ s/ //gs;
	$seq =~ s/[\-\.\~]//gs;
	$seq = lc $seq;
	$seq =~ s/u/t/g;

	if ($seq && $sid)
	{
	    print "$sid\t", chars_in_string(\$seq,$dictP), "\n";
	}
    }
}

sub chars_in_string {

# Counts the number certain characters in a string. Only
# the characters included in the incoming dictionary are 
# counted. So this can work for amino acids, bases or .. 

# In: Pointer to string
#     Pointer to dictionary hash
# Out: Integer

    local ($stringP,$dictP) = @_;
    my $i = 0;

    foreach $ch (keys %$dictP) {
	$i += $$stringP =~ s/$ch/$ch/g;
    }

    return $i;
}

sub base_dict_init {

# Returns a dictionary that says whether a character is a 
# valid base symbol. 

# In:  Nothing
# Out: Pointer to hash

    local %dict;

    $dict{"A"} = $dict{"a"} = 1;
    $dict{"G"} = $dict{"g"} = 1;
    $dict{"C"} = $dict{"c"} = 1;
    $dict{"U"} = $dict{"u"} = 1;
    $dict{"T"} = $dict{"t"} = 1;
    $dict{"R"} = $dict{"r"} = 1;
    $dict{"Y"} = $dict{"y"} = 1;
    $dict{"W"} = $dict{"w"} = 1;
    $dict{"S"} = $dict{"s"} = 1;
    $dict{"M"} = $dict{"m"} = 1;
    $dict{"K"} = $dict{"k"} = 1;
    $dict{"H"} = $dict{"h"} = 1;
    $dict{"D"} = $dict{"d"} = 1;
    $dict{"V"} = $dict{"v"} = 1;
    $dict{"B"} = $dict{"b"} = 1;
    $dict{"N"} = $dict{"n"} = 1;

    return \%dict;
}

sub complement {
#
# Returns the complement of a sequence with preservation of case 
# complemented ambiguity codes. 
#
    local ($seq) = @_;
    my (%dict); undef %dict; 
    my ($cseq) = "";

    $dict{"A"} = "U"; $dict{"a"} = "u";
    $dict{"G"} = "C"; $dict{"g"} = "c";
    $dict{"C"} = "G"; $dict{"c"} = "g";
    $dict{"U"} = "A"; $dict{"u"} = "a";
    $dict{"T"} = "A"; $dict{"t"} = "a";
    $dict{"R"} = "Y"; $dict{"r"} = "y";
    $dict{"Y"} = "R"; $dict{"y"} = "r";
    $dict{"W"} = "W"; $dict{"w"} = "w";
    $dict{"S"} = "S"; $dict{"s"} = "s";
    $dict{"M"} = "K"; $dict{"m"} = "k";
    $dict{"K"} = "M"; $dict{"k"} = "m";
    $dict{"H"} = "D"; $dict{"h"} = "d";
    $dict{"D"} = "H"; $dict{"d"} = "h";
    $dict{"V"} = "B"; $dict{"v"} = "b";
    $dict{"B"} = "V"; $dict{"b"} = "v";
    $dict{"N"} = "N"; $dict{"n"} = "n";

    foreach (reverse split(//,$seq)) {
        $cseq .= $dict{$_};
    }

    return $cseq;
}
