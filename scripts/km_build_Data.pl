########################################################################
########################################################################
use strict;
use Data::Dumper;
use SeedUtils;
use Getopt::Long;
use IPC::Run;
use IO::Handle;

my $usage = "usage: km_build_Data -d DataDir -k 8\n";
my $dataD;
my $k = 8;
my $rast_dirs;
my $sort_mem = "4G";

my $use_pub_seed = 0;     #### default is set to core_seed ###
my $rc  = GetOptions('d=s' => \$dataD,
		     'p'   => \$use_pub_seed,
		     'r=s' => \$rast_dirs,
		     'm=s' => \$sort_mem,
		     'k=i' => \$k);

my $min_reps_required = $ENV{KM_MIN_REPS_REQUIRED} || 5;

# properties vs normal-search is determined by the presence of 
# $dataD/genomes and $dataD/properties;

if ((! $rc) || (! -d $dataD) || (! $k))
{ 
    print STDERR $usage; exit ;
}
(((-s "$dataD/genomes") || (-s "$dataD/properties")) &&
 (! ((-s "$dataD/genomes") && (-s "$dataD/properties"))))
    || die "you need $dataD to contain exactly one of { genomes,properties }";

unlink("$dataD/final.kmers");
unlink("$dataD/kmer.table.mem_map");   ### force rebuild of memory map
# print STDERR "Building kmers of size $k\n";

my $distinct_signatures   = 0;     # number of distinct signature Kmers
my $distinct_functions    = {};    # functions having a signature
my $seqs_with_func        = {};    # hash giving number of seqs with a given function
my $seqs_with_a_signature = {};    # hash with keys of seqIDs - used to get number seqs with signatures

my $g_to_o = &build_otu_index($dataD); 

&build_function_index($dataD,$use_pub_seed,$rast_dirs);     # builds the function index 

my $sz = &build_reduced_kmers($dataD,
		     $use_pub_seed,
		     $rast_dirs,
		     $g_to_o,
		     \$distinct_signatures,
                     $distinct_functions,
                     $seqs_with_func,
                     $seqs_with_a_signature);

if ($ENV{KM_SKIP_WEIGHTED_SCORE_COMPUTATION})
{
    symlink("reduced.kmers", "$dataD/final.kmers") or die "Symlink reduced.kmers to $dataD/final.kmers failed: $!";
}
else
{
    &make_final($dataD,
		$seqs_with_a_signature,
		$distinct_signatures,
		$distinct_functions,
		$seqs_with_func);
}

open(SZ,">$dataD/size") || die "could not write size = $sz: $!";
print SZ "$sz\n";
close(SZ);

sub build_otu_index {
    my($dataD) = @_;
    my $g_to_o = {};
    my %counts;
    if (-s "$dataD/genomes") 
    {
	foreach $_ (`cat $dataD/genomes`)
	{
	    chomp;
	    my($name,$genome) = split(/\t/,$_);
	    if ($name =~ /^(\S+)\s(\S+)/)
	    {
		my $genus = $1;
		my $species = $2;
		if ($species !~ /^sp\.?$/i)
		{
		    if ($genus !~ /^(other|unclassified)/)
		    {
			$counts{"$genus $species"}++;
			$g_to_o->{$genome} = "$genus $species";
		    }
		}
	    }
	}
    }
    else   # if a property search
    {
	(-s "$dataD/properties") || die "you need a file $dataD/genomes or a file $dataD/properties";
	foreach $_ (`cat $dataD/properties`)
	{
	    chomp;
	    my($property,$genome) = split(/\t/,$_);
	    if (defined($genome) && defined($property))
	    {
		$counts{$property}++;
		$g_to_o->{$genome} = $property;
	    }
	}
    }
    open(WTS,">$dataD/otu.occurrences") || die "could not open $dataD/otu.occurrences: $!";
    my $nxt = 0;
    foreach my $gs_or_prop (sort keys(%counts))
    {
	print WTS join("\t",($nxt,$counts{$gs_or_prop},$gs_or_prop)),"\n";
	$nxt++;
    }
    close(WTS);

    system "cut -f1,3 $dataD/otu.occurrences > $dataD/otu.index";
    return $g_to_o;
}

sub build_function_index {
    my($dataD,$use_pub_seed,$rast_dirs) = @_;

    my %good_roles;
    my %good_funcs;
    if (open(my $fh, "<", "$dataD/subsystem.roles"))
    {
	while (<$fh>)
	{
	    chomp;
	    $good_roles{$_}++;
	}
	close($fh);
    }
    if (open(my $fh, "<", "$dataD/additional.funcs"))
    {
	while (<$fh>)
	{
	    chomp;
	    $good_funcs{$_}++;
	}
	close($fh);
    }

    open(FI,">$dataD/function.index") || die "could not open $dataD/function.index :$!";
    if (-s "$dataD/genomes")
    {
	my $cmd;
	if ($use_pub_seed)
	{
	    $cmd = "cut -f2 $dataD/genomes | svr_all_features peg | svr_function_of";
	}
	elsif ($rast_dirs)
	{
	    $cmd = "cat $rast_dirs/*/*functions";
	}
	else
	{
	    $cmd = "cut -f2 $dataD/genomes | genomes_to_fids | grep peg | fids_to_functions | cut -f2,3";
	}
	open(ASSIGNMENTS,"$cmd |")
	    || die "cannot access assignments: $!";

	my %funcs;
	while (defined($_ = <ASSIGNMENTS>))
	{
	    chomp;
	    my($fid, $func) = split(/\t/);
	    my $gid = SeedUtils::genome_of($fid);
	    if ($fid && $gid && $func)
	    {
		my $stripped = km_strip_func_comment($func);
		$funcs{$stripped}->{$gid}++;
	    }
	}
	close(ASSIGNMENTS);
	print STDERR Dumper(\%funcs);

	my $nxt = 0;
	foreach my $f (sort keys(%funcs))
	{
#	    if ($funcs{$f} > 5) #  && ((! &SeedUtils::hypo($f)) || ($f =~ /FIG/)))

	    my $gcount = keys %{$funcs{$f}};
	    my $ok;
	    if ($gcount >= $min_reps_required)
	    {
		$ok++;
		print STDERR "Adding function $f due to gcount>=$min_reps_required\n";
	    }
	    elsif ($f =~ /\#/)
	    {
		$ok++;
		print STDERR "Adding function $f due to containing #\n";
	    }
	    else
	    {
		my @r = SeedUtils::roles_of_function($f);
		if (grep { $_ }  @good_roles{@r})
		{
		    print STDERR "Adding function $f due to subsystem role membership\n";
		    $ok++;
		}

		if (!$ok)
		{
		    if ($good_funcs{$f})
		    {
			print STDERR "Adding function $f due to additional-funcs membership\n";
			$ok++;
		    }
		}
	    }
	    if ($ok)
	    {
		print FI "$nxt\t$f\n";
		$nxt++;
	    }
	}
    }
    else   # we are in a property Data directory
    {
	my $sets = {};
	(-s "$dataD/properties") || die "something is wrong";
	foreach $_ (`cat $dataD/properties`)
	{
	    if ($_ =~ /^(\S+)\t(\S+)$/)
	    {
		my($p,$g) = ($1,$2);
		$sets->{$p}->{$g} = 1;
	    }
	    else
	    {
		die "Bad entry in properties file: $_";
	    }
	}
	my $nxt = 0;
	foreach my $p (sort keys(%$sets))
	{
	    print FI join("\t",($nxt,$p)),"\n";
	    $nxt++;
	}
    }
    close(FI);
}

sub km_strip_func_comment
{
    my($func) = @_;
    my $stripped = $func;
    if (!$ENV{KM_BUILD_NO_STRIP})
    {
	$stripped = &SeedUtils::strip_func_comment($stripped);
    }
    return $stripped;
}

sub build_reduced_kmers {
    my($dataD,
       $use_pub_seed,
       $rast_dirs,
       $g_to_o,
       $distinct_signatures,
       $distinct_functions,
       $seqs_with_func,
       $seqs_with_a_signature) = @_;

    my %to_oI;
    foreach $_ (`cat $dataD/otu.occurrences`)
    {
	if ($_ =~ /^(\S+)\t(\S+)\t(\S.*\S)/)
	{
	    $to_oI{$3} = $1;
	}
    }

    my %to_fI;
    foreach $_ (`cat $dataD/function.index`)
    {
	if ($_ =~ /^(\S+)\t(\S.*\S)/)
	{
	    $to_fI{$2} = $1;
	}
    }

    my %id_to_fI;
    &load_id_to_fI($dataD,$use_pub_seed,$rast_dirs,\%to_fI,\%id_to_fI);
    my $properties = -s "$dataD/properties";
###
#   First, we build a file containing [kmer,fI,oI,off-set,sequence-ID]
#   sorted by kmer
###
    my $sort = "sort";
    if (-x "/scratch/olson/gnu-sort")
    {
	$sort = "/scratch/olson/gnu-sort --parallel 10";
    }
    open(RAW,"| $sort -S $sort_mem  > $dataD/sorted.kmers") || die "could not open $dataD/sorted.kmers: $!";
#    open(RAW,"| sort -T . -S 80G  > $dataD/sorted.kmers") || die "could not open $dataD/sorted.kmers: $!";
    my $seqID=0;
    my @genomes = map { chomp; $_ } ($properties ? `cut -f2 $dataD/properties` : `cut -f2 $dataD/genomes`);

    foreach my $g (sort @genomes)
    {
	my $cmd;
	if ($use_pub_seed)
	{
	    $cmd = "echo '$g' | svr_all_features peg | svr_translations_of";
	}
	elsif ($rast_dirs)
	{
	    $cmd = "flatten_fasta <  $rast_dirs/$g/Features/peg/fasta";
	}
	else
	{
	    $cmd = "echo '$g' | genomes_to_fids | grep peg | fids_to_protein_sequences -fasta 0 | cut -f2,3";
	}

	open(P, "$cmd |") or die "cannot open $cmd: $!";
	while (<P>)
	{
	    if ($_ =~ /^(fig\|\d+\.\d+\.peg\.\d+)\t(\S.*\S)$/)
	    {
		chomp;
		my $seq = $2;
		my $lseq = length($seq);
		$seq = uc($seq);
		my $id = $1;
		($id =~ /^fig\|(\d+\.\d+)/) || die "bad peg $_";
		my $oI = $to_oI{$g_to_o->{$g}};
		my $fI = $id_to_fI{$id};
		my $show;
		    # if ($id eq 'fig|83334.1.peg.2108')
		    # {
		    # 	print Dumper($id, $oI, $fI);
		    # 	$show = 1;
		    # }
		if (! defined($fI)) { $fI = $id_to_fI{"hypothetical protein"} }  # default is hypothetical
		if (defined($fI))
		{
		    $seqs_with_func->{$fI}++;                                                                # NFj
			    
		    for (my $i=0; ($i < (length($seq) - $k + 1)); $i++)
		    {
			my $kmer = substr($seq,$i,$k);

			# if ($show)
			# {
			#     print "GOT $kmer $fI $id\n";
			# }
			if ($kmer !~ /[^ACDEFGHIKLMNPQRSTVWY]/)
			{
			    print RAW join("\t",($kmer,
						 $fI,
						 $oI,
						 $lseq-$i,
						 $seqID)) . "\n";
			}
			else
			{
			    # print "bad kmer $kmer in $id\n";
			}
		    }
		    $seqID++;
		}
		else
		{
		    # warn "No fI for $id \n";
		}
	    }
	    else
	    {
		# print "Bad data $_\n";
	    }
	}
	close(P);
    }
    if (!close(RAW))
    {
	die "Sort failed: close returns $?";
    }

### Now, we reduce sets of adjacent lines with the same kmer to a single line (or none)
###
    open(RAW,"<$dataD/sorted.kmers") || die "could not open sorted.kmers: $!";
    
    open(REDUCED,">$dataD/reduced.kmers") || die "could not open reduced kmers: $!";
    open(REJECTED,">$dataD/rejected.kmers") || die "could not open rejected kmers: $!";
#     my $last = <RAW>;

#     while ($last && ($last =~ /^(\S+)/))
#     {
# 	my $curr = $1;
# 	my @set;
# 	while ($last)
# 	{
# 	    if ($last =~ /^(\S+)\t(\S*)\t(\S*)\t(\S*)\t(\S+)$/)
# 	    {
# 		if ($1 eq $curr)
# 		{
# 		    push(@set,[$2,$3,$4,$5]);
# 		    $last = <RAW>;
# 		}
# 		else
# 		{
# 		    last;
# 		}
# 	    }
# 	    else
# 	    {
# 		die "Invalid line $. in merged data\n$last\n";
# 	    }
# 	}
# 	&process_set($curr,\@set,\*REDUCED,$seqs_with_a_signature,$distinct_signatures,$distinct_functions,\*REJECTED);
#     }

    my @set;
    my $cur;
    my $size = 0;

    my $write_kmers;
    if ($ENV{KM_SKIP_WEIGHTED_SCORE_COMPUTATION})
    {
	$write_kmers = \&write_kmers_for_final_output;
    }
    else
    {
	$write_kmers = \&write_kmers_for_weighted_score_computation;
    }
       
    while (defined(my $row = <RAW>))
    {
	chomp $row;
	my($kmer, $fI, $oI, $off, $seqID) = split(/\t/, $row);

	if ($kmer ne $cur)
	{
	    &process_set($cur,\@set,\*REDUCED,$seqs_with_a_signature,$distinct_signatures,$distinct_functions,\*REJECTED, $write_kmers, \$size) if @set;
	    @set = ();
	    $cur = $kmer;
	}
	push(@set,[$fI, $oI, $off, $seqID]);
    }
    &process_set($cur,\@set,\*REDUCED,$seqs_with_a_signature,$distinct_signatures,$distinct_functions,\*REJECTED, $write_kmers, \$size) if @set;
    close(REDUCED);
    close(REJECTED);
    if (!$ENV{KM_KEEP_INTERMEDIATES})
    {
	unlink("$dataD/sorted.kmers");
    }
    return $size;
}
sub process_set {
    my($kmer,$set,$fh,$seqs_with_a_signature,$distinct_signatures,$distinct_functions,
       $rejected_fh, $write_kmers, $sizeP) = @_;

    my %funcs;
    my $tot = 0;
    foreach my $tuple (@$set)
    {
	my($fI,$oI,$off,$seqID) = @$tuple;
	$tot++;
	$funcs{$fI}++;
    }

    my $seqID;
    my @tmp = sort { $funcs{$b} <=> $funcs{$a} } keys(%funcs);

    # if ($kmer eq 'AANVIITD')
    # {
    # 	my $thresh = 0.8 * $tot;
    # 	print Dumper($kmer, $set, \@tmp, \%funcs, $funcs{$tmp[0]}, $tot, $thresh);
    # }
	    
    if (defined($tmp[0]) && ($funcs{$tmp[0]} >= (0.8 * $tot)))
    {
	my $best_fI = $tmp[0];
#	my $func_wt = $funcs{$tmp[0]};

#	print "$kmer: best=$best_fI with $tmp[1]\n";

	my %otus;
	my $otu = '';
	my $seqs_containing_func = 0;
	my $seqs_containing_sig  = @$set;
	foreach my $tuple (@$set)
	{
	    my($fI,$oI,$off,$seqID) = @$tuple;
	    $seqs_with_a_signature->{$seqID} = 1;
	    if ($fI == $best_fI)
	    {
		#
		# Is this a bug - shouldn't seqs_containing_func increment even if no otu?
		# (Fixed below)
		#
		$seqs_containing_func++;
		if (defined($oI))
		{
		    $otus{$oI}++;
		}
	    }
	}
	@tmp = sort { $otus{$b} <=> $otus{$a} } keys(%otus);
	if (defined($tmp[0]) && ($otus{$tmp[0]} >= (0.8 * $tot)))
	{
	    $otu = $tmp[0];
	}
	my @offsets = sort { $a <=> $b } map { $_->[2] } @$set;
	my $median_off = $offsets[int(scalar @offsets / 2)];

	$$sizeP++;
	&$write_kmers($fh,
		      $kmer,
		      $median_off,
		      $best_fI,
		      $otu,
		      $seqs_containing_sig,
		      $seqs_containing_func);
# 	print $fh join("\t",($kmer,
# 			     $median_off,
# 			     $best_fI,
# 			     $otu,
# 			     $seqs_containing_sig,
# 			     $seqs_containing_func)),"\n";
	$$distinct_signatures++;
	$distinct_functions->{$best_fI} = 1;
    }
    else
    {
	print $rejected_fh join("\t", $kmer, $tot, map { ($_, $funcs{$_}) } @tmp,), "\n";
    }
}

sub write_kmers_for_weighted_score_computation
{
    my($fh, $kmer,
       $median_off,
       $best_fI,
       $otu,
       $seqs_containing_sig,
       $seqs_containing_func) = @_;

    print $fh join("\t",($kmer,
			 $median_off,
			 $best_fI,
			 $otu,
			 $seqs_containing_sig,
			 $seqs_containing_func)),"\n";
}

sub write_kmers_for_final_output
{
    my($fh, $kmer,
       $median_off,
       $best_fI,
       $otu,
       $seqs_containing_sig,
       $seqs_containing_func) = @_;

    print $fh join("\t",($kmer,
			 $median_off,
			 $best_fI,
			 100,
			 $otu)), "\n";
}

#
# Create mapping from peg id to function index.
#
sub load_id_to_fI {
    my($dataD,$use_pub_seed,$rast_dirs,$to_fI,$id_to_fI) = @_;

    my $properties = -s "$dataD/properties";
    my %genomes = map { ($_ =~ /^(\S[^\t]*\S)\t(\S+)/) ? ($2 => $1) : () } 
                  ($properties ? `cat $dataD/properties` : `cat $dataD/genomes`);
    foreach my $genome (keys(%genomes))
    {


	my $proc_elt = sub {
	    my($x) = @_;
	    if ($x =~ /^(\S+)\t(\S.*\S)/)
	    {
		my $peg = $1;
		my $func_or_prop = $properties ? $genomes{$genome} : $2;
		$func_or_prop    = km_strip_func_comment($func_or_prop);
		if (defined($to_fI->{$func_or_prop}))
		{
		    $id_to_fI->{$peg} = $to_fI->{$func_or_prop};
		}
	    }
	};

	my $cmd;
	my @files;
	if ($use_pub_seed)
	{
	    $cmd = "echo '$genome' | svr_all_features peg | svr_function_of";
	}
	elsif ($rast_dirs)
	{
	    @files = <$rast_dirs/$genome/*functions>;
	}
	else
	{
	    $cmd = "echo '$genome' | genomes_to_fids | grep peg | fids_to_functions | cut -f2,3";
	}
	if ($cmd)
	{
	    open(P, "$cmd|") or die "Cannot open $cmd: $!";
	    while (<P>)
	    {
		$proc_elt->($_);
	    }
	    close(P);
	}
	elsif (@files)
	{
	    for my $f (@files)
	    {
		open(P, "<", $f) or die "Cannot open $f: $!";
		while (<P>)
		{
		    $proc_elt->($_);
		}
		close(P);
	    }
	}

    }
}

sub make_final {
    my($dataD,
       $seqs_with_a_signature,
       $distinct_signatures,
       $distinct_functions,
       $seqs_with_func) = @_;

    my $size = 0;
    my %to_oI;
    foreach $_ (`cat $dataD/otu.occurrences`)
    {
	if ($_ =~ /^(\d+)\t\d+\t(\S+)/)
	{
	    $to_oI{$2} = $1;
	}
    }

    my $num_seqs_with_a_signature = keys(%$seqs_with_a_signature);
    my $num_distinct_functions = keys(%$distinct_functions);
    open(IN,"<$dataD/reduced.kmers") || die "could not open final.kmers: $!";
    open(OUT,">$dataD/final.kmers") || die "could not open $dataD/final.signatures: $!";

    print STDERR "num_seqs_with_a_signature=$num_seqs_with_a_signature\n";
    print STDERR "num_distinct_functions=$num_distinct_functions\n";
    print STDERR "distinct_signatures=$distinct_signatures\n";

    while (defined($_ = <IN>))
    {
	chop;
	my($kmer,$off,$fI,$otu,$seqs_containing_sig,$seqs_containing_func) = split(/\t/,$_);
	my $fI_wt = &compute_weight_of_signature($num_seqs_with_a_signature, 	# NSF
						 $distinct_signatures,		# KS
						 $num_distinct_functions,	# KF
						 $seqs_containing_sig,		# NSi
						 $seqs_with_func->{$fI},	# NFj
						 $seqs_containing_func);	# NSiFj
	print OUT join("\t",($kmer,
			     $off,
			     $fI,
			     $fI_wt,
			     $otu)), "\n";
		       # , $seqs_containing_sig, $seqs_with_func->{$fI})),"\n";
	$size++;
    }
    close(IN);
    if (!$ENV{KM_KEEP_INTERMEDIATES})
    {
	unlink("$dataD/reduced.kmers");
    }
    close(OUT);
    return $size;
}

 
# The meaning of the variables is now as follows:

# $NSF   := number of training sequences with at least one signature and an assigned function.
# $KS    := number of signatures of function.
# $KF    := number of distinct functions with at least one signature.
# $NSi   := number of training sequences with at least one occurrence of signature Si.
# $NFj   := number of training sequences having function Fj.
# $NSiFj := number of sequences with at least one occurrence of Si and having function Fj

# Logical consistency requires that the following identities be satisfied:

#     $NSF = Sum(i=1..$KS, Sum(j=1..$KF, $NSiFj));
#     $NSi = Sum(j=1..$KF, $NSiFj);
#     $NFj = Sum(i=1..$KS, $NSiFj);

sub compute_weight_of_signature {
    my ($NSF, $KS, $KF, $NSi, $NFj, $NSiFj) = @_;

    return sprintf("%0.4f",(log( ($NSiFj + 1.0) / ($NSi - $NSiFj + 1.0) ) +
                            log( ($NSF - $NFj + $KS) / ($NFj   + $KS) )));
}

sub compute_weight_of_prior {
    my ($NSF, $KF, $NFj) = @_;

    return (log( ($NFj + 1.0) / ($NSF + $KF - $NFj - 1.0) ) );
}
