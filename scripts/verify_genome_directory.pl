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

use DB_File;
use File::stat;
use File::Basename;
# $fig = new FIG;

use Digest::MD5 qw(md5_base64);

$usage = "verify_genome_directory [-no_fatal_stops] [-code=genetic_code_number] Dir";

my $trouble        =  0;
my $no_fatal_stops =  0;
my $code_number    = 11;
while (@ARGV && ($ARGV[0] =~ m/^-/)) {
    if    ($ARGV[0] =~ m/-no_fatal_stops/o) {
	$no_fatal_stops = 1;
    }
    elsif ($ARGV[0] =~ m/-code=(\d+)/)  {
	$code_number = $1;
    }
    else {
	$trouble = 1;
	warn "FATAL: Invalid argument $ARGV[0]\n";
    }
    shift @ARGV;
}
die "\n\n   usage: $usage\n\n" if $trouble;

($dir = shift @ARGV)
    || die "\n   usage: $usage\n\n";

$dir =~ s/\/$//o;
if (!-d $dir) {
    print STDERR "Subdirectory $dir does not exist --- trying $FIG_Config::organisms/$dir\n" if $ENV{SEED_VERBOSE};
    $dir = "$FIG_Config::organisms/$dir";
    if (!-d $dir) {
	die "Organism directory $dir does not exist --- aborting";
    }
}



my %is_special;
if (-s qq($dir/special_pegs)) {
    %is_special = map { m/^(\S+)\s+([^\n]+)/o ? ($1 => $2) : () } &FIG::file_read(qq($dir/special_pegs));
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#...Set up Genetic Code tables...
#-----------------------------------------------------------------------
if (-s qq($dir/GENETIC_CODE)) {
    my $genetic_code_metadata = &FIG::file_head(qq($dir/GENETIC_CODE), 1);
    if ($genetic_code_metadata =~ m/^(\d+)/o) {
	$code_number = $1;
	print STDERR "Using genetic code $code_number\n" if $ENV{VERBOSE};
    }
    else {
	die "Could not handle contents of $dir/GENETIC_CODE: $_";
    }
}
$trans_tblP = &FIG::genetic_code($code_number);

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#...Add ambiguous codons that _could_ be STARTs...
#-----------------------------------------------------------------------
my @ambigs = qw(a c g t u m r w s y k b d h v n x);

my $startP = {};
foreach my $x (@ambigs) {
    foreach my $y (@ambigs) {
	foreach my $z (@ambigs) {
	    my $codon = $x.$y.$z;
	    if (grep { $_ =~ m/^[agt]tg/ } (&_expand_ambigs($codon))) {
		$startP->{uc($codon)} = qq(M);
	    }
	}
    }
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#...Also find ambiguous codons that could be STOPs...
#-----------------------------------------------------------------------
my $alt_STOPs = &FIG::genetic_code($code_number);
foreach my $x (@ambigs) {
    foreach my $y (@ambigs) {
	foreach my $z (@ambigs) {
	    my $codon = $x.$y.$z;
	    if (grep { &FIG::translate($_, $trans_tblP) eq qq(*) } (&_expand_ambigs($codon))) {
		$alt_STOPs->{uc($codon)} = qq(*);
	    }
	}
    }
}



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#...Extract Taxon-ID from directory path...
#-----------------------------------------------------------------------
my $genome_id = basename $dir;
if ($genome_id =~ /^(\d+\.\d+)/) {
    $genome_id = $1;
}
else {
    die "ERROR: path $dir does not end in a valid genome id\n";
}

if (! (-s "$dir/GENOME"))   {
    print STDOUT "ERROR: For $dir, missing GENOME file\n";
    $bad{"no_genome"} = 1;
} else {
    $genome = &FIG::file_head(qq($dir/GENOME), 1);
    chomp $genome; 
}

if (! (-s "$dir/PROJECT"))  { 
    print STDOUT "ERROR: For $dir, missing PROJECT file\n";
    $bad{"no_project"} = 1;
}

if (! (-s "$dir/TAXONOMY")) { 
    print STDOUT "ERROR: For $dir, missing TAXONOMY file\n";
    $taxonomy = $domain = "Unknown";
    $bad{"no_taxonomy"} = 1;
} else { 
    $taxonomy = &FIG::file_head(qq($dir/TAXONOMY), 1);
    chomp $taxonomy;
    $taxonomy =~ m/^([^;]+)/;
    $domain   =  $1;
}


if (!-d "$dir/Features") {
    die "ERROR: For $dir, features directory $dir/Features does not exist";
}
else {
    if (! opendir(TMP,$dir)) {
	++$bad{no_org_dir};
	print STDOUT "ERROR: could not open $dir";
    }
    else {
	@contigs = grep { $_ !~ /(\~|\.nhr|\.nin|\.nnd|\.nni|\.nsd|\.nsi|\.nsq|\.btree|\.properties)$/o } grep { $_ =~ /^contig/ } readdir(TMP);
	closedir(TMP) || warn "Could not close directory $dir";
	
	if (@contigs == 0) {
	    ++$bad{contigs_missing_files};
	    print STDOUT "ERROR: $dir contains no contigs files\n";
	}
	else {
	    $size = 0;

	    #
	    # Sum up the size of all the contigs to see if we need to tie
	    # %subs to a db file.
	    #

	    my $contig_total = 0;
	    foreach $file (@contigs) {
		my $s = stat("$dir/$file");
		if ($s) {
		    $contig_total += $s->size;
		}
	    }

	    # the DB_File below changes an 11' run into 19', so i've upped the
            # threshold to 20GB to turn this off for now until a better
	    # fix is done (read tbls first and hash it rather than contigs)

	    if ($contig_total > 20_000_000_000) {
		$tied_file = "$FIG_Config::temp/verify_db.$$";
		tie %seq_of, 'DB_File', $tied_file or warn "Could not tie to $tied_file: $!\n";
	    }
	    
	    foreach $file (@contigs) {
		if (not open(TMP, "<$dir/$file")) {
		    ++$bad{contigs_unopenable_file};
		    print STDOUT "ERROR: Could not read-open $dir/$file";
		}
		else {
		    $last_id = "";
		    while (($id, $seqP, $rest) = &FIG::read_fasta_record(\*TMP)) {
			if (not $id) {
			    print STDOUT "ERROR: For $dir, file $file, record $., null contig ID after $last_id";
			}
			
			if (not defined($len_of_contig{$id})) {
			    $seq_of{$id} = $$seqP;
			    $size += $len_of_contig{$id} = length($$seqP);
			}
			else {
			    ++$bad{contigs_duplicates};
			    print STDOUT "ERROR: For $dir, file $file, record $., duplicate contig ID $id ($rest)\n"
				, substr($$seqP, 0, 50), "\n";
			}
			
			$checksum = md5_base64($$seqP);
			if (defined($used{$checksum})) {
			    print STDERR "WARNING: For $dir, file $file, record $., "
				, "$id has the same checksum as file $used{$checksum} --- it may be a duplicate\n";
			}
			else {
			    $used{$checksum} = "$file, $id";
			}
			
			$tmp =  $$seqP;
			$tmp =~ s/[ACGTRYMKSWHBVDN]//ig;
			if ($tmp) {
			    ++$bad{contigs_bad_chars};
			    print STDOUT "ERROR: For $dir, file $file, record $., "
				, "contig $id contains invalid chars:\n$tmp\n\n";
			}
			
			$last_id = $id;
		    }
		}
	    }
	    
	    if ((-s "$dir/COMPLETE") && ($size < 300000)) {
		print STDERR "WARNING: For $dir, contigs for $dir contain only $size bp of data\n";
	    }
	    
	    if (!-d "$dir/Features/peg")
	    {
		++$bad{peg_no_directory};
		print STDOUT "ERROR: For $dir, no peg directory $dir/Features/peg\n";
	    }
	    else
	    {
		if (!-s "$dir/Features/peg/tbl") {
		    ++$bad{peg_no_tbl};
		    print STDOUT "ERROR: For $dir, no peg tbl file $dir/Features/peg/tbl\n";
		}
		else {
		    if (not open(TMP, "<$dir/Features/peg/tbl")) {
			++$bad{peg_unopenable_tbl};
			print STDOUT "ERROR: For $dir, could not read-open $dir/Features/peg/tbl";
		    }
		    else {
			while (defined($entry = <TMP>)) {
			    chomp $entry;
			    ($id, $loc) = split(/\t/, $entry);
			    
			    if ($id =~ /fig\|(\d+\.\d+)\.([^\.]+)\.\d+/) {
				my($org, $type) = ($1, $2);
				
				if ($org ne $genome_id) {
				    ++$bad{invalid_peg_genome};
				    print STDOUT "ERROR: PEG id $id does not match taxon-ID $genome_id in $dir/Features/peg/tbl\n";
				}
				
				if ($type ne 'peg') {
				    ++$bad{invalid_peg_type};
				    print STDOUT "ERROR: PEG id $id has wrong type\n";
				}
			    }
			    
			    if (defined($loc_of_peg{$id})) {
				++$bad{peg_tbl_dup_id};
				print STDOUT "ERROR: For $dir, $id, duplicate PEG ID in: $entry\n"; 
			    }
			    else {
				$loc_of_peg{$id} = $loc;
				
				@exons = split(/,/, $loc);
				if ($exons[-1] =~ m/^(\S+)_(\d+)_(\d+)$/) {
				    
				    ($contig, $beg, $end) = ($1, $2, $3);
				    
				    unless (defined($is_special{$id}) && ($is_special{$id} =~ m/pseudogene/io)) {
					if ($beg < $end) {
					    $codon = uc &FIG::extract_seq(\%seq_of, join(qq(_), ($contig, ($end-2), $end)));
					    if ((&FIG::translate($codon, $alt_STOPs) ne qq(*))
						) {
						if (($end + 3) <= $len_of_contig{$contig}) {
						    $codon = uc &FIG::extract_seq(\%seq_of, join(qq(_), ($contig, ($end+1), ($end+3))));
						    if ((&FIG::translate($codon, $alt_STOPs) ne qq(*)) &&
							(not &possibly_truncated($contig, $beg, $end))
							) {
							print STDERR "PEG STOP missing for $id at $loc\n";
							++$warnings{"STOP-$codon"};
							
							if ($no_fatal_stops) {
							    ++$warnings{peg_tbl_stop_missing};
							}
							else {
							    ++$bad{peg_tbl_stop_missing};
							}
						    }
						    else {
							### ++$warnings{$codon};
							++$warnings{peg_tbl_stop_not_included};
							print STDERR "PEG STOP not included for $id at $loc\n";
						    }
						}
					    }
					    else {
						### ++$warnings{$codon};
					    }
					}
					else {
					    $codon = uc &FIG::extract_seq(\%seq_of, join(qq(_), ($contig, ($end+2), $end)));
					    if ((&FIG::translate($codon, $alt_STOPs) ne qq(*)) &&
						(not &possibly_truncated($contig, $beg, $end))
						) {
						
						if (($end - 3) >= 1) {
						    
						    $codon = uc &FIG::extract_seq(\%seq_of, join(qq(_), ($contig, ($end-1), ($end-3))));
						    if (&FIG::translate($codon, $alt_STOPs) ne qq(*)) {
							print STDERR "PEG STOP missing for $id at $loc\n";
							++$warnings{"STOP-$codon"};
							
							if ($no_fatal_stops) {
							    ++$warnings{peg_tbl_stop_missing};
							}
							else {
							    ++$bad{peg_tbl_stop_missing};
							}
						    }
						    else {
							### ++$warnings{$codon};
							++$warnings{peg_tbl_stop_not_included};
							print STDERR "PEG STOP not included for $id at $loc\n";
						    }
						}
					    }
					    else {
						### ++$warnings{$codon};
					    }
					}
				    }
				}
				else {
				    ++$bad{peg_tbl_unparsable_terminal_exon};
				}
				
				foreach $exon (@exons) {
				    $exon =~ m/^(\S+)_(\d+)_(\d+)$/;
				    ($contig, $beg, $end) = ($1, $2, $3);
				    if ($taxonomy =~ m/Eukaryot/) {
					if (not defined($len_of_contig{$contig})) {
					    print STDERR "WARNING: For $dir, $id, $loc refers to undefined contig $contig\n";
					}
					else {
					    if (&FIG::max($beg,$end) > $len_of_contig{$contig}) {
						print STDERR "WARNING: For $dir, $id, "
						    , "$exon in $loc exceeds contig length $len_of_contig{$contig}\n";
					    }
					}
				    }
				    else {
					if (not defined($len_of_contig{$contig})) {
					    ++$bad{peg_undef_contig_ref};
					    print STDOUT "ERROR: For $dir, $id, $loc refers to undefined contig $contig\n";
					}
					else {
					    if (&FIG::max($beg,$end) > $len_of_contig{$contig}) {
						++$bad{peg_out_of_bounds};
						print STDOUT "ERROR: For $dir, $id, "
						    , "$exon in $loc exceeds contig length $len_of_contig{$contig}\n";
					    }
					}
				    }
				}
			    }
			}
			close(TMP) || warn "could not close $dir/Features/peg/tbl";
		    }
		    
		    if (not open(TMP, "<$dir/Features/peg/fasta")) {
			++$bad{peg_unopenable_fasta};
			print STDOUT "ERROR: For $dir, could not read-open $dir/Features/peg/fasta";
		    }
		    else {
			while (($id, $seqP) = &FIG::read_fasta_record(\*TMP)) {
			    if (defined($len_of_peg{$id})) {
				++$bad{peg_fasta_dup_id};
				print STDOUT "ERROR: For $dir, $id duplicated in PEG fasta file\n";
			    }
			    else {
				$len_of_peg{$id} = length($$seqP);
				if (not defined($loc_of_peg{$id})) {
				    ++$bad{peg_fasta_no_tbl};
				    print STDOUT "ERROR: For $dir, $id is in PEG fasta, but not in tbl\n";
				}
			    }
			    
			    if ($$seqP =~ s/\*$//) { ++$bad{peg_terminal_stop}; }
			    
			    if ($$seqP =~ s/\*//g) {
				++$bad{peg_embedded_stop};
				print STDOUT "ERROR: PEG $id contains embedded STOPs\n";
			    }
			    
			    if ($$seqP =~ m/[^ABCDEFGHIJKLMNPQRSTUVWXxYZ]/) {
				$_ =  $$seqP;
				$_ =~ s/[ABCDEFGHIJKLMNPQRSTUVWXxYZ]//go;
				++$bad{peg_bad_chars}; 
				print STDOUT "ERROR: PEG $id contains invalid characters: $_\n";
			    }
			    
			    if (($num_acgt = ($$seqP =~ tr/acgtACGT//)) && ($num_acgt > 0.9*$len_of_peg{$id})) {
				++$bad{peg_untranslated};
				print STDOUT "ERROR: PEG $id appears to be DNA, not protein\n";
			    }
			}
			close(TMP) || warn "could not close $dir/Features/peg/fasta";
			print STDERR "Read ", (scalar keys %len_of_peg), " pegs from $dir/Features/peg/fasta\n" 
			    if $ENV{SEED_VERBOSE};
		    }
		    
		    foreach $id (sort keys %loc_of_peg) {
			if (not defined($len_of_peg{$id})) {
			    ++$bad{peg_tbl_no_fasta};
			    print STDOUT "ERROR: For $dir, $id is in peg tbl, but not in fasta\n";
			}
		    }
		}
		
		if (!-d "$dir/Features/rna") {
		    print STDERR "WARNING: For $dir, no rna directory\n";
		}
		else {
		    if (!-s "$dir/Features/rna/tbl") {
			++$bad{rna_no_tbl};
			print STDOUT "ERROR: For $dir, no rna tbl file $dir/Features/rna/tbl\n";
		    }
		    else {
			if (not open(TMP, "<$dir/Features/rna/tbl")) {
			    ++$bad{rna_unopenable_tbl};
			    print STDOUT "ERROR: For $dir, could not read-open $dir/Features/rna/tbl";
			}
			else {
			    while (defined($entry = <TMP>)) {
				chomp $entry;
				($id, $loc) = split(/\t/, $entry);
				
				if ($id =~ /fig\|(\d+\.\d+)\.([^\.]+)\.\d+/) {
				    my($org, $type) = ($1, $2);
				    
				    if ($org ne $genome_id) {
					++$bad{invalid_rna_genome};
					print STDOUT "ERROR: RNA id $id does not match taxon-ID $genome_id in $dir/Features/rna/tbl\n";
				    }
				    
				    if ($type ne 'rna') {
					++$bad{invalid_rna_type};
					print STDOUT "ERROR: RNA id $id has wrong type\n";
				    }
				}
				
				if (defined($loc_of_rna{$id})) {
				    ++$bad{rna_tbl_dup_id};
				    print STDOUT "ERROR: For $dir, $id, duplicate RNA ID in: $entry\n"; 
				}
				else {
				    $loc_of_rna{$id} = $loc;
				    foreach $exon (split(/,/, $loc)) {
					$exon =~ m/^(\S+)_(\d+)_(\d+)$/;
					($contig, $beg, $end) = ($1, $2, $3);
					if (not defined($len_of_contig{$contig})) {
					    ++$bad{rna_undef_contig_rna};
					    print STDOUT "ERROR: For $dir, $id, $loc refers to undefined contig $contig\n";
					}
					else {
					    if (&FIG::max($beg,$end) > $len_of_contig{$contig}) {
						++$bad{rna_out_of_bounds};
						print STDOUT "ERROR: For $dir, $id, $exon in $loc exceeds contig length $len_of_contig{$contig}\n";
					    }
					}
				    }
				}
			    }
			    close(TMP) || warn "could not close $dir/Features/rna/tbl";
				
			    if (not open(TMP, "<$dir/Features/rna/fasta")) {
				++$bad{rna_unopenable_fasta};
				print STDOUT "ERROR: For $dir, could not read-open $dir/Features/rna/fasta";
			    }
			    else {
				while (($id, $seqP) = &FIG::read_fasta_record(\*TMP)) {
				    if (defined($len_of_rna{$id})) {
					++$bad{peg_fasta_dup_id};
					print STDOUT "ERROR: For $dir, $id duplicated in RNA fasta file\n";
				    }
				    else {
					$len_of_rna{$id} = length($$seqP);
					($x = $$seqP) =~ s/[acgtn]//gi;
					($y = $$seqP) =~ s/[acgtumrwsykbdhvn]//gi;
					if ($y || (length($x) > (0.1) * length($$seqP))) {
					    ++$warnings{rna_fasta_bad_chars};
					    print STDERR "WARNING: RNA $id does not look like an RNA sequence\n";
					}
				    }
				    
				    if (not defined($loc_of_rna{$id})) {
					++$bad{rna_fasta_no_tbl};
					print STDOUT "ERROR: For $dir, $id is in rna fasta, but not in tbl\n";
				    }
				}
				close(TMP) || warn "could not close $dir/Features/rna/fasta";
				print STDERR "Read ", (scalar keys %len_of_rna), " rnas from $dir/Features/rna/fasta\n" 
				    if $ENV{SEED_VERBOSE};
				
				foreach $id (sort keys %loc_of_rna) {
				    if (not defined($len_of_rna{$id})) {
					++$bad{rna_tbl_no_fasta};
					print STDOUT "ERROR: For $dir, $id is in rna tbl, but not in fasta\n";
				    }
				}
			    }
			}
		    }
		}
	    }
	}

	if (open(FF,"<$dir/assigned_functions")) {
	    while ($entry = <FF>) {
		chomp $entry;
		my($id,$fn) = split /\t/, $entry;
		
		if ($id =~ /fig\|(\d+\.\d+)\.(peg|rna)\.\d+/) {
		    my($org, $type) = ($1, $2);

		    if ($org ne $genome_id) {
			++$bad{invalid_peg_in_assignments};
			print STDOUT "ERROR: PEG id $id does not match taxon-ID $genome_id in $dir/assigned_functions\n";
		    }
		    elsif ($type eq "peg" and !$loc_of_peg{$id}) {
			++$bad{invalid_peg_in_assignments};
			print STDOUT "ERROR: PEG id $id in $dir/assigned_functions was not found in features\n";
		    }
		    elsif ($type eq "rna" and !$loc_of_rna{$id}) {
			++$bad{invalid_rna_in_assignments};
			print STDOUT "ERROR: RNA id $id in $dir/assigned_functions was not found in features\n";
		    }
		}
		elsif ($id =~ /fig\|(\d+\.\d+)\.(\w+)\.\d+/) {
		    my($org, $type) = ($1, $2);
		    if ($org ne $genome_id) {
		   	++$bad{invalid_featureg_in_assignments};
			print STDOUT "ERROR: Feature id $id does not match taxon-ID $genome_id in $dir/assigned_functions\n";
		    }
		}
		else {
		    ++$bad{invalid_assignments};
		    print STDOUT "ERROR: Invalid assignment in line $.: $entry\n";
		}
	    }
	}
	
	if (-f "$dir/assigned.functions") {
	    ++$bad{misnamed_assignments_file};
	    print STDOUT "ERROR: assigned.functions file should be named assigned_functions\n";
	}
    }
}

if (keys %warnings) {
    print STDERR "$dir has warnings (", join(", ", map { "$_=$warnings{$_}" } sort keys %warnings), ")\t$domain\t$genome\n";
}

if (keys %bad) {
    print STDOUT "$dir is corrupt (", join(", ", map { "$_=$bad{$_}" } sort keys %bad), ")\t$domain\t$genome\n";
}
else {
    print STDERR "$dir is O.K.\n";
}
# print STDERR "$0 has completed\n";

exit(scalar keys %bad);


sub possibly_truncated {
    my($contig, $beg, $end) = @_;

    return 1 if (($beg < $end) && ($end < 300));
    return 1 if (($beg > $end) && ($beg < 300));
    
    return 1 if (($beg < $end) && ($beg > ($len_of_contig{$contig} - 300)));
    return 1 if (($beg > $end) && ($end > ($len_of_contig{$contig} - 300)));
    
    return 0;
}

sub _expand_ambigs {
    my (@stack) = @_;
    print STDERR "Expanding ", join(", ", @stack), "\n"
	if ($ENV{VERBOSE} && (@stack > 1));
    
    my @out;
    while (@stack > 0) {
	# m = (a|c)
	if ($stack[0] =~ m/^([^m]*)m(.*)$/) {
	    shift( @stack );
	    unshift( @stack, ( "$1a$2", "$1c$2") );
	}
	
	# r = (a|g)
	if ($stack[0] =~ m/^([^r]*)r(.*)$/) {
	    shift( @stack );
	    unshift( @stack, ("$1a$2", "$1g$2") );
	}
	
	# w = (a|t)
	if ($stack[0] =~ m/^([^w]*)w(.*)$/) {
	    shift( @stack );
	    unshift( @stack, ("$1a$2", "$1t$2") );
	}
	
	# s = (c|g)
	if ($stack[0] =~ m/^([^s]*)s(.*)$/) {
	    shift( @stack );
	    unshift( @stack, ("$1c$2", "$1g$2") );
	}
	
	# y = (c|t)
	if ($stack[0] =~ m/^([^y]*)y(.*)$/) {
	    shift( @stack );
	    unshift( @stack, ("$1c$2", "$1t$2") );
	}
	
	# k = (g|t)
	if ($stack[0] =~ m/^([^k]*)k(.*)$/) {
	    shift( @stack );
	    unshift( @stack, ("$1g$2", "$1t$2") );
	}
	
	# v = (a|c|g)
	if ($stack[0] =~ m/^([^v]*)v(.*)$/) {
	    shift( @stack );
	    unshift( @stack, ("$1a$2", "$1c$2", "$1g$2") );
	}
	
	# h = (a|c|t)
	if ($stack[0] =~ m/^([^h]*)h(.*)$/) {
	    shift( @stack );
	    unshift( @stack, ("$1a$2", "$1c$2", "$1t$2") );
	}
	
	# d = (a|g|t)
	if ($stack[0] =~ m/^([^d]*)d(.*)$/) {
	    shift( @stack );
	    unshift( @stack, ("$1a$2", "$1g$2", "$1t$2") );
	}
	
	# b = (c|g|t)
	if ($stack[0] =~ m/^([^b]*)b(.*)$/) {
	    shift( @stack );
	    unshift( @stack, ("$1c$2", "$1g$2", "$1t$2") );
	}
	
	# (n|x) = (a|c|g|t)
	if ($stack[0] =~ m/^([^xn]*)[xn](.*)$/) {
	    shift( @stack );
	    unshift( @stack, ("$1a$2", "$1c$2", "$1g$2", "$1t$2") );
	}
	
	while ( (@stack > 0) && ($stack[0] !~ m/[mrwsykvhdbxn]/)) {
	    push( @out, shift(@stack) );
	}
	
	last if (@stack == 0);
    }
    
    return @out;
}
