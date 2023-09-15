# -*- perl -*-
#########################################################################
# Copyright (c) 2003-2008 University of Chicago and Fellowship
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
#########################################################################

package FIGV;

use Carp qw(carp cluck confess);
use strict;
use FIG;
use FIG_Config;
use SFXlate;
use SproutFIG;
use Tracer;
use Data::Dumper;
use vars qw($AUTOLOAD);
use StandaloneSims;
use DB_File;
use File::Basename;
use FileHandle;
use FileLocking qw(lock_file unlock_file);
use SeedUtils ();

sub new {
    my ($class, $org_dir, $low_level, $fig) = @_;
    print STDERR "FIGV::new => (", join(qq(, ), ($class, $org_dir, $low_level)), ")\n"
	if $ENV{FIG_VERBOSE};
    
    if (!ref($fig)) {
	if ($low_level && ($low_level =~ /sprout/i)) {
	    warn "New FIGV using SPROUT low-level" if $ENV{FIG_DEBUG};
	    $fig = new SproutFIG($FIG_Config::sproutDB, $FIG_Config::sproutData);
	}
	elsif (! $org_dir) {
	    warn "New FIGV using FIG low-level" if $ENV{FIG_DEBUG};
	    $fig = new FIG;
	    return $fig;
	}
	else {
	    warn "New FIGV, defaulting to FIG low-level" if $ENV{FIG_DEBUG};
	    $fig = new FIG;
	}
    }
    
    if (defined($org_dir)) {
	warn "New FIGV will use organism-directory $org_dir" if $ENV{FIG_VERBOSE};
	if (!-d $org_dir) {
	    confess "ERROR: Organism directory $org_dir does not exist";
	}
    }
    else {
	warn "New FIGV using system organism directories in $FIG_Config::organsims" if $ENV{FIG_VERBOSE};
    }
    
    my $self         = {};
    $self->{_fig}    = $fig;
    $self->{_orgdir} = $org_dir;
    $self->{_peer_sims_cache} = {};
    
    #
    # Determine genome-id. Use the GENOME_ID file if it exists; otherwise
    # we assume the directory name of the org_dir is the genome id.
    #
    my $genome_id;
    if (open(my $fh, "<", "$org_dir/GENOME_ID")) {
	$genome_id = <$fh>;
	chomp $genome_id;
	if ($genome_id !~ /^\d+\.\d+/) {
	    warn "Invalid genome id '$genome_id' in $org_dir/GENOME_ID";
	    return undef;
	}
	close($fh);
    }
    else {
	if ($org_dir =~ /(\d+\.\d+)$/) {
	    $genome_id = $1;
	}
	else {
	    warn "No GENOME_ID file found in $org_dir and the directory name is not a valid genome id";
	}
    }
    $self->{_genome} = $genome_id;
    
    return bless $self, $class;
}

sub DESTROY {
    my($self) = @_;
    print STDERR "FIGV destroy $self $$\n";
    $self->SUPER::DESTROY if $self->can("SUPER::DESTROY");
}

sub genome_info {
    my ($self) = @_;
    
    my $info = $self->{_fig}->genome_info;
    my $id = $self->genome_id;
    push(@$info, [ $id, "Private: ".$self->genus_species($id), $self->genome_szdna($id), $self->genome_domain($id), $self->genome_pegs($id), $self->genome_rnas($id), $self->is_complete($id), $self->taxonomy_of($id) ]);
    
    return $info;
}

sub genome_id {
    return $_[0]->{_genome};
}

# return the path of the organism directory
sub organism_directory {
    my ($self,$genome) = @_;
    
    if (not $genome) {
	#...Grandfather in the old (and wrong!) behavior of returning the virtual orgdir
	return $self->{_orgdir};
    }
    elsif ($genome eq $self->{_genome}) {
	return $self->{_orgdir};
    }
    
    return $self->{_fig}->organism_directory($genome);
}

sub is_complete {
    my($self,$genome) = @_;
    my $fig = $self->{_fig};

    my $newG = $self->{_genome};
    if ($genome ne $newG) {
	return $fig->is_complete($genome);
    }
    return 1;
}

#
# Redirect any method invocations that we don't implement out to the
# underlying FIG object.
#
sub AUTOLOAD {
    my($self, @args) = @_;
    
    if (ref($self) ne "FIGV") {
	confess "BAD FIGV object passed to AUTOLOAD";
    }
    
    no strict 'refs';
    
    my $meth = $AUTOLOAD;
    $meth =~ s/.*:://;
    my $fmeth = "FIG::$meth";
    
    my $fig = $self->{_fig};
#    my $args = Dumper(\@args);
    if (wantarray) {
	my @res = $fig->$meth(@args);
#	my @res = &$fmeth($self, @args);
#	warn "FIGV invoke $meth($args) returns\n", Dumper(\@res);
	return @res;
    }
    else {
	my $res = $fig->$meth(@args);
#	my $res = &$fmeth($self, @args);
#	warn "FIGV invoke $meth($args) returns\n", Dumper($res);
	return $res;
    }
}

sub FIG {
    my($self) = @_;
    return $self;
}

#
# This should be hoisted to FIGM
#
sub sort_fids_by_taxonomy {
    my($self,@fids) = @_;
    
    return map     { $_->[0] }
           sort    { $a->[1] cmp $b->[1] }
           map     { [$_,$self->taxonomy_of($self->genome_of($_))] }
           grep { ! $self->is_deleted_fid($_) }
           @fids;
}

sub get_basic_statistics {
    my($self, $genome) = @_;
    
    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};
    
    if ($genome ne $newG) {
	return $fig->get_basic_statistics($genome);
    }
    
    #
    # Check cache.
    #
    my $cache = "$newGdir/cache.basic_statistics";
    my $fh = new FileHandle($cache);
    if ($fh)
    {
	my $stats = {};
	while (<$fh>)
	{
	    chomp;
	    my($k, $v) = split(/\t/);
	    $stats->{$k} = $v;
	}
	close($fh);
	return $stats;
    }
    
    my $subsystem_data = $self->get_genome_subsystem_data($genome);
    
    my %sscount = map { $_->[0] => 1 } @$subsystem_data;
    my $nss=scalar(keys(%sscount));
    
    my $statistics = {
	num_subsystems => $nss,
	num_contigs    => scalar($self->all_contigs($genome)),
	num_basepairs  => $self->genome_szdna($genome),
	genome_name    => $self->genus_species($genome),
	genome_domain  => $self->genome_domain($genome),
	genome_pegs    => $self->genome_pegs($genome),
	genome_rnas    => $self->genome_rnas($genome),
	genome_version => $self->genome_version($genome)
	};

    $fh = new FileHandle(">$cache");
    if ($fh) {
	while (my($k, $v) = each %$statistics) {
	    print $fh join("\t", $k, $v), "\n";
	}
	close($fh);
    }
    
    return $statistics;
}


sub get_peg_statistics {
    my ($self, $genome) = @_;
    
    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};
    
    if ($genome ne $newG) {
	return $fig->get_peg_statistics($genome);
    }
    
#... Check cache.
    my $cache = "$newGdir/cache.peg_statistics";
    my $fh = new FileHandle($cache);
    if ($fh) {
	my $stats = {};
	while (<$fh>) {
	    chomp;
	    my($k, $v) = split(/\t/);
	    $stats->{$k} = $v;
	}
	close($fh);
	return $stats;
    }
    
    my $subsystem_data = $self->get_genome_subsystem_data($genome);
    my $assignment_data = $self->get_genome_assignment_data($genome);
    
    my $hypo_sub = 0;
    my $hypo_nosub = 0;
    my $nothypo_sub = 0;
    my $nothypo_nosub = 0;
    my %in = map { $_->[2] => 1 } @$subsystem_data;
    my $in = keys(%in);
    
    my %sscount = map { $_->[0] => 1 } @$subsystem_data;
    
    foreach $_ (@$assignment_data) {
	my($peg,$func) = @$_;
	if (! $self->is_deleted_fid($peg)) {
	    my $is_hypo = &FIG::hypo($func);
	    
	    if    ($is_hypo && $in{$peg})           { $hypo_sub++ }
	    elsif ($is_hypo && ! $in{$peg})         { $hypo_nosub++ }
	    elsif ((! $is_hypo) && (! $in{$peg}))   { $nothypo_nosub++ }
	    elsif ((! $is_hypo) && $in{$peg})       { $nothypo_sub++ }
	}
    }
    my $tot = $hypo_sub + $nothypo_sub + $hypo_nosub + $nothypo_nosub;
    
    my ($fracHS, $fracNHS, $fracHNS, $fracNHNS);
    
    if ($tot == 0) {
	$fracHS = sprintf "%.2f", 0.0;
	$fracNHS = sprintf "%.2f", 0.0;
	$fracHNS = sprintf "%.2f", 0.0;
	$fracNHNS = sprintf "%.2f", 0.0;
    }
    else {
	$fracHS = sprintf "%.2f", $hypo_sub / $tot * 100;
	$fracNHS = sprintf "%.2f", $nothypo_sub / $tot * 100;
	$fracHNS = sprintf "%.2f", $hypo_nosub / $tot * 100;
	$fracNHNS = sprintf "%.2f", $nothypo_nosub / $tot * 100;
    }
    
    my $statistics = {
	hypothetical_in_subsystem => $hypo_sub,
	hypothetical_not_in_subsystem => $hypo_nosub,
	non_hypothetical_in_subsystem => $nothypo_sub,
	non_hypothetical_not_in_subsystem => $nothypo_nosub,
	hypothetical_in_subsystem_percent => $fracHS,
	hypothetical_not_in_subsystem_percent => $fracHNS,
	non_hypothetical_in_subsystem_percent => $fracNHS,
	non_hypothetical_not_in_subsystem_percent => $fracNHNS
	};
    
    $fh = new FileHandle(">$cache");
    if ($fh) {
	while (my($k, $v) = each %$statistics) {
	    print $fh join("\t", $k, $v), "\n";
	}
	close($fh);
    }
    
    return $statistics;
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# To retrieve a subsystem in FIGV, we create the subsystem as normal via $fig->get_subsystem,
# then insert the row for the virtual org dir we are processing.
#-----------------------------------------------------------------------
sub get_subsystem {
    my($self,$ssa) = @_;
    
    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};
    
    my $ss = $fig->get_subsystem($ssa);
    return undef unless $ss;
    
    $self->load_ss_data();
    
    my $bindings = $self->{_ss_bindings}->{$ssa};
    my $variant = $self->{_ss_variants}->{$ssa};
    
    unless ($bindings) {
	$variant = '*-1';
	$bindings = {};	
    }
    
    $ss->add_virtual_genome($self->genus_species(), $newG, $variant, $bindings);
    
    return $ss;
}

sub get_variant_and_bindings {
    my($self, $ssa) = @_;
    
    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};
    
    $self->load_ss_data();
    
    my $bindings = $self->{_ss_bindings}->{$ssa};
    my $variant = $self->{_ss_variants}->{$ssa};
    
    unless ($bindings) {
	$variant = '*-1';
	$bindings = {};	
    }
    
    return ($variant, $bindings);
}

sub active_subsystems {
    my($self, $genome, $all) = @_;
    
    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};
    
    if ($genome ne $newG) {
	return $fig->active_subsystems($genome, $all);
    }
    
    $self->load_ss_data();
    
    my $slist = {};
    
    if ($self->{_ss_variants}) {
	%{$slist} = %{$self->{_ss_variants}};
    }

    if (not $all) {
	for my $ss (keys %$slist) {
	    my $var = $slist->{$ss};
	    delete $slist->{$ss} if $var eq 0 or $var eq -1;
	}
    }
    return $slist;
}

sub peg_to_subsystems {
    my($self, $peg) = @_;
    
    if ($self->is_deleted_fid($peg)) { return () }
    
    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};
    
    if ($peg !~ /fig\|$newG\.peg/) {
	return $fig->peg_to_subsystems($peg);
    }
    
    $self->load_ss_data();
    
    my $variant = $self->{_ss_variants};
    my %sub = map { $_ => 1 }
              grep { $variant->{$_} !~ /^(-1)|0$/ }
              map { $_->[0] }
              $self->peg_to_roles_in_subsystems($peg);
    return sort keys(%sub);
}

sub peg_to_roles_in_subsystems {
    my($self,$peg) = @_;
    
    if ($self->is_deleted_fid($peg)) { return () }
    
    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};
    
    if ($peg !~ /fig\|$newG\.peg/) {
	return $fig->peg_to_roles_in_subsystems($peg);
    }
    
    $self->load_ss_data();
    
    my $ret  = $self->{_ss_peg_index}->{$peg};

    return $ret ? @$ret : ();
}

sub subsystems_for_peg {
    my($self,$peg) = @_;
    return $self->peg_to_roles_in_subsystems($peg);
}

sub subsystems_for_peg_complete {
    my ($self, $peg) = @_;
    
    my $newG = $self->{_genome};
    my $fig  = $self->{_fig};
    
    if ($peg !~ /fig\|$newG\.peg/) {
	return $fig->subsystems_for_peg_complete($peg);
    }
    
    $self->load_ss_data();
    
    my $ret = $self->{_ss_peg_index}->{$peg};
    if ($ret) {
	return map { [ $_->[0], $_->[1], $self->{_ss_variants}->{$_->[0] }, $peg ] } @$ret;
    } else {
	return ();
    }
}

sub genomes {
    my($self,$complete) = @_;
    my $fig = $self->{_fig};
    
    return ($self->{_genome},$fig->genomes($complete));
}

sub genus_species {
    my($self,$genome) = @_;
    
    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};
    
    if (($genome eq $newG) && open(GENOME,"<$newGdir/GENOME")) {
	my $x = <GENOME>;
	close(GENOME);
	chop $x;
	return $x;
    }
    else {
	return $fig->genus_species($genome);
    }
}

sub get_genome_assignment_data {
    my($self,$genome) = @_;
    
    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};
    
    if ($genome eq $newG) {
	my @fnfiles = <$newGdir/proposed*functions>;
	if (@fnfiles == 0) {
	    @fnfiles = <$newGdir/assigned*functions>;
	}
	
	if (@fnfiles == 0) {
	    return [];
	}
	my %assign;
	foreach $_ (`cat @fnfiles`) {
	    if ( $_ =~ /^(fig\|\d+\.\d+\.peg\.\d+)\t(\S.*\S)/) {
		my($fid,$func) = ($1,$2);
		if (! $self->is_deleted_fid($fid)) {
		    $assign{$fid} = $func;
		}
	    }
	}
	return [map { [$_,$assign{$_}] } sort { &FIG::by_fig_id($a,$b) } keys(%assign)];
    }
    else {
	return $fig->get_genome_assignment_data($genome);
    }
}

sub org_of {
    my($self,$peg) = @_;
    
    if ($peg =~ /^fig\|(\d+\.\d+)\.peg\.\d+/) {
	return $self->genus_species($1);
    }
    return "";
}

sub get_genome_subsystem_data {
    my($self,$genome) = @_;
    
    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};
    
    if ($genome eq $newG) {
	my %operational = map { (($_ =~ /^(\S.*\S)\t(\S+)/) && (($2 ne '-1') && ($2 ne '0'))) ? ($1 => 1) : () }
		          `cat $newGdir/Subsystems/subsystems`;
	           
	return [grep { ! $self->is_deleted_fid($_->[2]) }
	        map { (($_ =~ /^(\S[^\t]+\S)\t(\S[^\t]*\S)\t(\S+)/) && $operational{$1} ) ? [$1,$2,$3] : () }
		`cat $newGdir/Subsystems/bindings`];
    }
    else {
	return $fig->get_genome_subsystem_data($genome);
    }
}

sub get_genome_subsystem_count {
    my($self,$genome) = @_;
    
    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};
    
    if ($genome eq $newG) {
	my $count = 0;
	my ($entry, $vc);
	open(SUBSYSTEMS, "<$newGdir/Subsystems/subsystems");
	while (defined($entry = <SUBSYSTEMS>)) {
	    chomp $entry;
	    (undef, $vc) = split /\t/, $entry;
	    if ($vc != -1) { ++$count; }
	}
	close(SUBSYSTEMS);
	return $count;
    }
    else {
	return $fig->get_genome_subsystem_count($genome);
    }
}

sub orgname_of_orgid {
    my($self,$genome) = @_;
    
    return $self->genus_species($genome);
}

sub orgid_of_orgname {
    my($self,$genome_name) = @_;
    
    my @genomes = $self->genomes('complete');
    my $i;
    for ($i=0; ($i < @genomes) && ($genome_name ne $self->genus_species($genomes[$i])); $i++) {}
    return ($i == @genomes) ? undef : $genomes[$i];
}

sub genus_species_domain {
    my($self,$genome) = @_;
    
    return [$self->genus_species($genome),$self->genome_domain($genome)];
}

sub protein_subsystem_to_roles {
    my ($self,$peg,$subsystem) = @_;
    
    if ($self->is_deleted_fid($peg)) { return [] }
    
    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};
    
    if (&FIG::genome_of($peg) ne $newG) {
	return $fig->protein_subsystem_to_roles($peg,$subsystem);
    }
    else {
	my @roles = map { (($_ =~ /^([^\t]+)\t([^\t]+)\t(\S+)$/) && ($1 eq $subsystem) && ($3 eq $peg)) ?
			  $2 : () } `cat $newGdir/Subsystems/bindings`;
	my %roles = map { $_ => 1 } @roles;
	return [sort keys(%roles)];
    }
}

sub contig_lengths {
    my ($self, $genome) = @_;
    
    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};
    
    if ($genome ne $newG) {
	return $fig->contig_lengths($genome);
    }
    else {
	my $contig_lengths = $self->{_contig_len_index};
	
	$self->load_contigs_index();
	if (!tied(%$contig_lengths)) {
	    $self->load_contig_len_cache();
	    return $self->{_contig_len_cache};
	}
	return $contig_lengths;
    }
}

sub load_contig_len_cache {
    my($self) = @_;
    
    return if ref $self->{_contig_len_cache};
    
    my $newGdir = $self->{_orgdir};
    
    my $contig_lengths = {};
    if (open(CONTIGS,"<$newGdir/contigs")) {
	local $/ = "\n>";
	while (defined(my $x = <CONTIGS>)) {
	    chomp $x;
	    if ($x =~ />?(\S+)[^\n]*\n(.*)/s) {
		my $id = $1;
		my $seq = $2;
		$seq =~ s/\s//gs;
		$contig_lengths->{$id} = length($seq);
	    }
	}
	close(CONTIGS);
    }
    $self->{_contig_len_cache} = $contig_lengths;
}

sub contig_ln {
    my ($self, $genome, $contig) = @_;
    
    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};
    
    if ($genome ne $newG) {
	return $fig->contig_ln($genome, $contig);
    }
    else {
	my $contig_len = $self->{_contig_len_index};
	
	$self->load_contigs_index();
	if (tied(%$contig_len)) {
	    return $contig_len->{$contig};
	}
	
	$self->load_contig_len_cache();
	return $self->{_contig_len_cache}->{$contig};
    }
}

sub contigs_of {
    my ($self, $genome) = @_;
    
    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};
    
    if ($genome ne $newG) {
	return $fig->contigs_of($genome);
    }
    else {
	my @out;
	$self->load_contigs_index();
	
	my $contigs = $self->{_contigs_index};
	if (tied(%$contigs)) {
	    return keys %$contigs;
	}
	
	$self->load_contig_len_cache();
	
	return keys %{$self->{_contig_len_cache}};
    }
}

sub get_dna_seq {
    my ($self, $fid) = @_;

    my $genome    = $self->genome_of( $fid );
    my @locations = $self->feature_location( $fid );

    my $seq = $self->dna_seq($genome, @locations);

    return $seq;
}

=head3 dna_seq

usage: $seq = dna_seq($genome,@locations)

Returns the concatenated subsequences described by the list of locations.  Each location
must be of the form

    Contig_Beg_End

where Contig must be the ID of a contig for genome $genome.  If Beg > End the location
describes a stretch of the complementary strand.

=cut
#: Return Type $;
sub dna_seq {
    my($self,$genome,@locations) = @_;
    
    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};
    
    if ($genome ne $newG) {
	return $fig->dna_seq($genome, @locations);
    }
    
    my $contigs = $self->{_contigs_index};
    if (!tied %$contigs) {
	$self->load_contig_seq();
	$contigs = $self->{_contigs_seq_cache};
    }
    
    my(@pieces,$loc,$contig,$beg,$end,$ln,$rdbH);
    
    @locations = map { split(/,/,$_) } @locations;
    @pieces = ();
    foreach $loc (@locations) {
        if ($loc =~ /^(\S+)_(\d+)_(\d+)$/) {
            ($contig,$beg,$end) = ($1,$2,$3);
	    my $seq = $contigs->{$contig};
	    
            $ln = length($seq);
	    
            if (! $ln) {
                print STDERR "$genome/$contig: could not get length\n";
                return "";
            }
	    
            if (&FIG::between(1,$beg,$ln) && &FIG::between(1,$end,$ln)) {
                if ($beg < $end) {
                    push(@pieces, substr($seq, $beg - 1, ($end - $beg) + 1));
                }
                else {
                    push(@pieces, &FIG::reverse_comp(substr($seq, $end - 1, ($beg - $end) + 1)));
                }
            }
        }
    }
    return lc(join("",@pieces));
}

sub load_contig_seq {
    my($self) = @_;
    
    return if ref($self->{_contigs_seq_cache});
    
    my $newGdir = $self->{_orgdir};
    
    my $contigs = {};
    
    if (open(CONTIGS,"<$newGdir/contigs")) {
	local $/ = "\n>";
	while (defined(my $x = <CONTIGS>)) {
	    chomp $x;
	    if ($x =~ />?(\S+)[^\n]*\n(.*)/s) {
		my $id = $1;
		my $seq = $2;
		$seq =~ s/\s//gs;
		$contigs->{$id} = $seq;
	    }
	}
	close(CONTIGS);
    }
    $self->{_contigs_seq_cache} = $contigs;
}

sub genome_szdna {
    my ($self, $genome) = @_;
    
    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};
    
    if ($genome ne $newG) {
	return $fig->genome_szdna($genome);
    }
    else {
	my $contig_lens = $self->contig_lengths($genome);
	my $tot = 0;
	while ( my($contig,$len) = each %$contig_lens) {
	    $tot += $len;
	}
	return $tot;
    }
}

sub genome_version {
    my ($self, $genome) = @_;
    
    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};
    
    if ($genome ne $newG) {
	return $fig->genome_version($genome);
    }
    else {
	return "$genome.0";
    }
}

sub genome_pegs {
    my ($self, $genome) = @_;
    
    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};
    
    if ($genome ne $newG) {
	return $fig->genome_pegs($genome);
    }
    else {
	my @tmp = $self->all_features($genome,"peg");
	my $n = @tmp;
	return $n;
    }
}

sub genome_rnas {
    my ($self, $genome) = @_;
    
    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};
    
    if ($genome ne $newG) {
	return $fig->genome_rnas($genome);
    }
    else {
	my @tmp = $self->all_features($genome,"rna");
	my $n = @tmp;
	return $n;
    }
}

sub genome_domain {
    my ($self, $genome) = @_;
    
    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};
    
    if ($genome ne $newG) {
	return $fig->genome_domain($genome);
    }
    else {
	my $tax = $self->taxonomy_of($genome);
	return ($tax =~ /^([^ \t;]+)/) ? $1 : "unknown";
    }
}

sub genes_in_region {
    my($self,$genome,$contig,$beg,$end) = @_;
    
    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};
    
    if ($genome ne $newG) {
	return $fig->genes_in_region($genome,$contig,$beg,$end);
    }
    else {
	$self->load_feature_indexes();
	
	#... Use the recno index if exists.
	my $maxV = 0;
	my $minV = 1000000000;
	my $genes = [];
	
	my $recnos = $self->{_feat_recno};
	
	if (ref($recnos)) {
	    while (my($ftype, $list) = each (%$recnos)) {
		#... Look up start/end of this contig in the btree index.
		
		my $inf = $self->{_feat_btree}->{$ftype}->{$contig};
		my($istart, $iend);
		if ($inf) {
		    ($istart, $iend) = split(/$;/, $inf);
		}
		else {
		    $istart = 0;
		    $iend = $#$list;
		}
		
		for (my $idx = $istart; $idx <= $iend; $idx++) {
		    my($fid, $fcontig, $fbeg, $fend, $fstrand) = split(/$;/, $list->[$idx]);
		    if ($contig eq $fcontig and &overlaps($beg, $end, $fbeg, $fend)) {
			$minV = &FIG::min($minV,$fbeg,$fend);
			$maxV = &FIG::max($maxV,$fbeg,$fend);
			push(@$genes,$fid);
		    }
		}
	    }
	}
	else {
	    &load_tbl($self);
	    my $tblH = $self->{_tbl};
	    while ( my($fid,$tuple) = each %$tblH) {
		if (($tuple->[0]->[0] =~ /^(\S+)_(\d+)_\d+$/) && ($1 eq $contig)) {
		    my $beg1 = $2;
		    my $last = @{$tuple->[0]} - 1;
		    if (($tuple->[0]->[$last] =~ /^(\S+)_\d+_(\d+)$/) && ($1 eq $contig)) {
			my $end1 = $2;
			if (&overlaps($beg,$end,$beg1,$end1)) {
			    $minV = &FIG::min($minV,$beg1,$end1);
			    $maxV = &FIG::max($maxV,$beg1,$end1);
			push(@$genes,$fid);
			}
		    }
		}
	    }
	}
	return ($genes,$minV,$maxV);
    }
}

sub overlaps {
    my($b1,$e1,$b2,$e2) = @_;
    
    if ($b1 > $e1) { ($b1,$e1) = ($e1,$b1) }
    if ($b2 > $e2) { ($b2,$e2) = ($e2,$b2) }
    return &FIG::between($b1,$b2,$e1) || &FIG::between($b2,$b1,$e2);
}

sub all_contigs {
    my($self,$genome) = @_;
    return $self->contigs_of($genome);
}

sub all_features {
    my($self,$genome,$type) = @_;
    if (not defined($type)) { $type = q(); }
    
    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};
    
    if ($genome ne $newG) {
	warn "SEED genome $genome\n" if $ENV{FIG_VERBOSE};
	return $fig->all_features($genome,$type);
    }
    else {
	warn "Virtual genome $genome\n" if $ENV{FIG_VERBOSE};
	warn "Loading feature indices" if $ENV{FIG_VERBOSE};
	$self->load_feature_indexes();
	
	my %contigs;
	my $btrees = $self->{_feat_btree};
	
	if (ref($btrees)) {
	    warn "B-tree already loaded" if $ENV{FIG_VERBOSE};
	    
	    my $btree = $btrees->{$type};
	    return sort { &FIG::by_fig_id($a, $b) }
	    	grep { /^fig/ } keys %$btree;
	}
	else {
	    warn "Loading contig B-tree" if $ENV{FIG_VERBOSE};
	    
	    &load_tbl($self);
	    my $tblH = $self->{_tbl};
	    
	    return sort { 
		&FIG::by_fig_id($a,$b)
		} grep {
		    #...NOTE: Matches all feature types if $type is the null string
		    (not $type) || (($_ =~ /^fig\|\d+\.\d+\.([^\.]+)/) && ($1 =~ m/$type/))
		    } keys(%$tblH);
	}
    }
}

sub all_features_detailed_fast {
    my($self,$genome, $regmin, $regmax, $contig) = @_;
    
    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};
    
    if ($genome ne $newG) {
	return $fig->all_features_detailed_fast($genome, $regmin, $regmax, $contig);
    }
    else {
	$self->load_feature_indexes();
	my $feat_details = [];
	my $recnos = $self->{_feat_recno};
	
	if (ref($recnos)) {
	    while (my($ftype, $list) = each (%$recnos)) {
		#
		# Look up start/end of this contig in the btree index.
		#
		
		my $inf = $self->{_feat_btree}->{$ftype}->{$contig};
		my($istart, $iend);
		if ($inf) {
		    ($istart, $iend) = split(/$;/, $inf);
		}
		else {
		    $istart = 0;
		    $iend = $#$list;
		}
		
		for (my $idx = $istart; $idx <= $iend; $idx++) {
		    my($fid, $fcontig, $fbeg, $fend, $fstrand) = split(/$;/, $list->[$idx]);
		    
		    if (not defined($regmin) or not defined($regmax) or not defined($contig) or
			(($contig eq $fcontig) and
			 ($fbeg < $regmin and $regmin < $fend) or ($fbeg < $regmax and $regmax < $fend) or ($fbeg > $regmin and $fend < $regmax)))
		    {
			my($loc, $index, @aliases) = split(/$;/, $self->{_feat_btree}->{$ftype}->{$fid});

			my $function = $self->function_of($fid);			
			push(@$feat_details,[$fid, $loc, join(",", @aliases), $ftype, $fbeg, $fend, $function,'master','']);
		    }
		}
	    }
	}
	else {
	    &load_tbl($self);
	    my $tblH = $self->{_tbl};
	    while ( my($fid,$tuple) = each %$tblH) {
		if ($fid =~ /^fig\|\d+\.\d+\.(\S+)\.\d+/) {
		    my $type = $1;
		    my($ctg, $min, $max, $strand) = &SeedUtils::boundaries_of($tuple->[0]);

		    next if (defined($contig) and $contig ne $ctg);
			
		    if (not defined($regmin) or not defined($regmax) or
			($min < $regmin and $regmin < $max) or ($min < $regmax and $regmax < $max) or ($min > $regmin and $max < $regmax))
		    {
			my $function = $self->function_of($fid);
			push(@$feat_details,[$fid,join(",", @{$tuple->[0]}),join(",",@{$tuple->[1]}),$type,$min,$max,$function,'master','']);
		    }
		}
	    }
	}
	return $feat_details;
    }
}

sub compute_clusters {
    # Get the parameters.
    my ($self, $pegList, $subsystem, $distance) = @_;
    if (! defined $distance) {
        $distance = 5000;
    }

    my($peg,%by_contig);
    foreach $peg (@$pegList)
    {
        my $loc;
        if ($loc = $self->feature_location($peg))
        {
            my ($contig,$beg,$end) = &FIG::boundaries_of($loc);
            my $genome = &FIG::genome_of($peg);
            push(@{$by_contig{"$genome\t$contig"}},[($beg+$end)/2,$peg]);
        }
    }

    my @clusters = ();
    foreach my $tuple (keys(%by_contig))
    {
        my $x = $by_contig{$tuple};
        my @pegs = sort { $a->[0] <=> $b->[0] } @$x;
        while ($x = shift @pegs)
        {
            my $clust = [$x->[1]];
            while ((@pegs > 0) && (abs($pegs[0]->[0] - $x->[0]) <= $distance))
            {
                $x = shift @pegs;
                push(@$clust,$x->[1]);
            }

            if (@$clust > 1)
            {
                push(@clusters,$clust);
            }
        }
    }
    return sort { @$b <=> @$a }  @clusters;
}

sub boundaries_of {
    my($self,@args) = @_;

    my $fig     = $self->{_fig};
    return $fig->boundaries_of(join(",", @args));
}



sub feature_location {
    my($self,$fid) = @_;

    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};

    if (($fid =~ /^fig\|(\d+\.\d+)\.([^.]+)/) && ($1 eq $newG))
    {
	my $ftype = $2;
	
	$self->load_feature_indexes();

	my $btree = $self->{_feat_btree}->{$ftype};
	if ($btree)
	{
	    my($loc, $idx, @aliases) = split(/$;/, $btree->{$fid});
	    return wantarray ? split(/,/, $loc) : $loc;
	}
	else
	{
	    &load_tbl($self);
	    if (my $x = $self->{_tbl}->{$fid})
	    {
		if (wantarray)
		{
		    return @{$x->[0]};
		}
		else
		{
		    return join(",",@{$x->[0]});
		}
	    }
	    else
	    {
		return undef;
	    }
	}
    }
    else
    {
	return scalar $fig->feature_location($fid);
    }
}

sub feature_location_bulk {
    my($self,$fids) = @_;
    
    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};
    
    my @ids;
    my @out;
    for my $fid (@$fids) {
	if ($self->is_virtual_feature($fid)) {
	    push(@out, [$fid, scalar $self->feature_location($fid)]);
	}
	else {
	    push(@ids, $fid);
	}
    }
    push(@out, $fig->feature_location_bulk(\@ids));
    return @out;
}


sub add_feature {
    my( $self, $user, $genome, $type, $location, $aliases, $sequence, $called_by_file, $calling_method) = @_;
    #... optional $called_by_file and $calling_method arguments help support RAST.
    
    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};
    
    if ($genome && $genome ne $newG) {
	return $fig->add_feature($user, $genome, $type, $location, $aliases, $sequence);
    }
    
    # perform sanity checks
    unless ($user && $genome && $type && $location && $sequence) {
	print STDERR "SEED error: add_feature failed due to missing parameter\n";
	return undef;
    }
    
    if ( $genome !~ /^\d+\.\d+$/ ) {
        print STDERR "SEED error: add_feature failed due to bad genome id: $genome\n";
        return undef;
    }
    
    if ( $type !~ /^[0-9A-Za-z_]+$/ ) {
        print STDERR "SEED error: add_feature failed due to bad type: $type\n";
        return undef;
    }
    
    if ( length ( $location ) > 5000 ) {
        print STDERR "SEED error: add_feature failed because location is over 5000 char:\n";
        print STDERR "$location\n";
        return undef;
    }
    
    my @loc  = split( /,/, $location );
    my @loc2 = grep { defined($_->[0]) && $_->[1] && $_->[2] }    
               map  { [ $_ =~ m/^(.+)_(\d+)_(\d+)$/ ] }
               @loc;
    
    if ( (@loc2 == 0) || ( @loc != @loc2 ) ) {
        print STDERR "SEED error: add_feature failed because location is missing or malformed:\n";
        print STDERR "$location\n";
        return undef;
    }
    
    if ( my @bad_names = grep { length( $_->[0] ) > 96 } @loc2 ) {
        print STDERR "SEED error: add_feature failed because location contains a contig name of over 96 char:\n";
        print STDERR join( ", ", @bad_names ) . "\n";
        return undef;
    }
    
    # create type directory if it does not exist
    unless (-d "$newGdir/Features/$type") {
	$self->{_fig}->verify_dir("$newGdir/Features/$type");
	(-d "$newGdir/Features/$type")
	    || die qq(Feature directory \'$newGdir/Features/$type\' does not exist, and could not be created);
	
	open(TMP, qq(>$newGdir/Features/$type/tbl))
	    || die qq(Could not create empty \'$newGdir/Features/$type/tbl\');
	close(TMP);
	
	open(TMP, qq(>$newGdir/Features/$type/fasta))
	    || die qq(Could not create empty \'$newGdir/Features/$type/fasta\');
	close(TMP);
    }
    
    # create an id
    my $id = "fig|$genome.$type.";
    my $feature_dir = "$newGdir/Features/$type";
    my $file = "$feature_dir/tbl";
    if (-f $file) {
	unless (open(FILE, "<$file")) {
	    print STDERR "SEED error: could not open tbl file: $@\n";
	    return undef;
	}
	my $entry;
        my $max = 0;
	while (defined($entry = <FILE>)) {
	    chomp $entry;
	    if ($entry =~ /^fig\|$genome\.$type\.(\d+)/) {
		my $curr_id = $1;
		if ($curr_id > $max) {
		    $max = $curr_id;
		}
	    }
	    else {
		confess qq(Could not parse $type tbl entry: $entry);
	    }
	}
	close FILE;
	++$max;
	$id .= $max;
    } else {
	$id .= "1";
    }
    
    # append to tbl file
    unless (open(FILE, ">>$file")) {
	print STDERR "SEED error: could not open tbl file: $@\n";
	return undef;
    }
    
    lock_file(\*FILE) || confess "cannot lock tbl file";
    seek(FILE,0,2)      || confess "failed to seek to the end of the file";
    $aliases =~ s/,/\t/g;
    print FILE "$id\t$location\t$aliases\n";
    close FILE;	
    chmod(0777,$file);
    
    my $typeH = $self->{_features}->{$type};
    if ($typeH) {
	$typeH->{$id} = 1;
    }
    
    my $tbl = $self->{_tbl};
    if ($tbl) {
	$tbl->{$id} = [[@loc],[split(/\t/,$aliases)]];
    }
    
    # append to fasta file
    $sequence =~ s/\s//g;
    
    $file = "$feature_dir/fasta";
    unless (open(FILE, ">>$file")) {
	print STDERR "SEED error: could not open fasta file: $@\n";
	return undef;
    }
    lock_file(\*FILE) || confess "cannot lock fasta file";
    seek(FILE,0,2)      || confess "failed to seek to the end of the file";
    print FILE ">$id\n$sequence\n";
    close FILE;
    chmod(0777,$file);
    
    if ($called_by_file && $calling_method) {
	#... append to called_by file
	$file = "$newGdir/called_by";
	unless (open(FILE, ">>$file")) {
	    print STDERR "SEED error: could not open called_by file: $@\n";
	    return undef;
	}
	else {
	    lock_file(\*FILE) || confess "cannot lock called_by file";
	    seek(FILE,0,2)      || confess "failed to seek to the end of the file";
	    print FILE "$id\t$calling_method\n";
	    close FILE;
	    chmod(0777,$file);
	}
    }
    
    #
    # If a btree was created for this, clear the ref to it and delete the files since
    # they are now invalid.
    #
    
    my $tie = $self->{_feat_tie}->{$type};
    my $btree = $self->{_feat_btree}->{$type};
    
    if ($tie) {
	untie $tie;
	delete $self->{$_}->{$type} for qw(_feat_tie _feat_btree _feat_ltie _feat_recno _feat_fasta);

	unlink("$feature_dir/tbl.btree");
	unlink("$feature_dir/tbl.recno");
    }
    
    if (-f "$feature_dir/fasta.norm.phr") {
	unlink(<$feature_dir/fasta.norm.*>);
    }
    
    # declare success
    return $id;
}

sub assign_function {
    my ($self, $fid, $user, $function, $confidence, $functions_file, $reason) = @_;
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#... optional $functions_file argument used to support RAST
#    (proposed_functions, proposed_non_ff_functions, etc.)
#-----------------------------------------------------------------------
    if (!$self->is_virtual_feature($fid)) {
	return $self->{_fig}->assign_function($fid, $user, $function, $confidence);
    }
    
    my $is_master = ($user =~ s/^master://io) ? 1 : 0;
    my $realuser = $user;  #  For annotation
    $user = 'master';      #  Actual assignments are treated as master assignments
    
    $confidence = $confidence ? $confidence : "";
    
    if (not $functions_file) {
	$functions_file = join(q(/), ($self->{_orgdir}, ($is_master ? q(assigned_functions) : q(proposed_user_functions))));
    }
    
    if ((! $self->is_real_feature($fid)) || (! $user)) { return 0 }
    my $genome = $self->genome_of($fid);
    
    $function =~ s/\s+/ /sg;  # No multiple spaces
    $function =~ s/^\s+//;    # No space at begining
    $function =~ s/\s+$//;    # No space at end
    $function =~ s/ ; /; /g;  # No space before semicolon
    
    if ($function =~ /^(.*?)[\!](.*)/) {
	my $kvs;
        ($function,$kvs) = ($1,$2);
	#
	# We don't have support for these any more
	#
	warn "Ignoring key/value pairs in assignment to $fid\n";
        if (0 and $kvs) {
            $kvs =~ s/^\s+//;
            $kvs =~ s/\s+$//;
            foreach my $kv (split(/\s+[\!\#]\s+/,$kvs)) {
                if ($kv =~ /^([A-Za-z0-9._\-\+\%]+)\s+\^\s+(.*)$/) {
                    my ($k,$v) = ($1,$2);
                    if ($v !~ /\S/) {
                        &replace_peg_key_value($self,$fid,$k,"");
                    }
                    else {
                        &replace_peg_key_value($self,$fid,$k,$v);
                    }
                }
                elsif ($kv =~ /^([A-Za-z0-9._\-\+\%]+)$/) {
                    &replace_peg_key_value($self,$fid,$1,1);
                }
            }
        }
    }
    
    my $status = 1;
    if ( open( TMP, ">>$functions_file" ) ) {
        lock_file(\*TMP)   || confess "cannot lock file \'$functions_file\'";
        seek(TMP,0,2)      || confess "failed to seek to the end of file \'$functions_file\'";
        print TMP "$fid\t$function\t$confidence\n";
        close(TMP);
        chmod(0777,$functions_file);
    }
    else {
        print STDERR "FAILED ASSIGNMENT: fid=$fid\tfunc=$function\tconfidence=$confidence\tfile=$functions_file\n";
        $status = 0;
    }
    
    # mdj:  force reload of functions to pick up new assignment
    $self->load_functions(1);
    
    #  We are not getting annotations logged.  So, we will impose it here.
    my $annotation = "Set master function to\n$function\n";
    $annotation .= $reason ? $reason.qq(\n) : q();
    $self->add_annotation( $fid, $realuser, $annotation );
    
    #
    # Mark the genome directory as in need of having bindings recomputed.
    #
    if (open(SS, "<$self->{_orgdir}/Subsystems/subsystems")) {
	while (<SS>) {
	    chomp;
	    my($sname, $v) = split(/\t/);
	    open(SS_FILE, ">$self->{_orgdir}/Subsystems/${sname}_bindings_need_recomputation");
	    close(SS_FILE);
	}
	close(SS);
    }
    
    return $status;
}

sub need_bindings_recomputed
{
    my($self, $ssa) = @_;
    return -f "$self->{_orgdir}/Subsystems/${ssa}_bindings_need_recomputation";
}

sub clear_need_bindings_recomputed
{
    my($self, $ssa) = @_;
    unlink "$self->{_orgdir}/Subsystems/${ssa}_bindings_need_recomputation";
}

sub add_annotation {
    my($self,$feature_id,$user,$annotation, $time_made) = @_;
    my($genome);

    if (!$self->is_virtual_feature($feature_id))
    {
	return $self->{_fig}->add_annotation($feature_id,$user,$annotation, $time_made);
    }

    $time_made = time unless $time_made =~ /^\d+$/;

    if ($self->is_deleted_fid($feature_id)) { return 0 }

#   print STDERR "add: fid=$feature_id user=$user annotation=$annotation\n";
    if ($genome = $self->genome_of($feature_id))
    {
	my $file = "$self->{_orgdir}/annotations";
        my $ma   = ($annotation =~ /^Set master function to/);

        if (open(TMP,">>$file"))
        {
            lock_file(\*TMP) || confess "cannot lock annotations";
            seek(TMP,0,2)      || confess "failed to seek to the end of the file";

            my $dataLine = "$feature_id\n$time_made\n$user\n$annotation" . ((substr($annotation,-1) eq "\n") ? "" : "\n");
            print TMP $dataLine . "//\n";
            close(TMP);
            chmod 0777, $file;

	    #
	    # Update local cache.
	    #
	    my $ann = $self->{_ann};
	    push(@{$ann->{$feature_id}}, [$feature_id, $time_made, $user, $annotation . "\n"]);
        }
    }
    return 0;
}

sub function_of {
    my($self,$fid) = @_;

    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};

    if (($fid =~ /^fig\|(\d+\.\d+)/) && ($1 eq $newG))
    {
	&load_functions($self);
	
	my $fn = $self->{_functions}->{$fid};
	if (wantarray)
	{
	    return ['master', $fn];
	}
	else
	{
	    return $fn;
	}
    }
    else
    {
	return $fig->function_of($fid);
    }
}

sub function_of_bulk {
    my($self,$fid_list) = @_;

    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};

    &load_functions($self);

    my $out = {};
    for my $fid (@$fid_list)
    {
	my $fn;
	if (($fid =~ /^fig\|(\d+\.\d+)/) && ($1 eq $newG))
	{
	    $fn = $self->{_functions}->{$fid};
	}
	else
	{
	    $fn = $fig->function_of($fid);
	}
	$out->{$fid} = $fn if defined($fn);
    }
    return $out;
}


=pod

find_features_by_annotation

Takes a reference to a hash of functions to find and an optional case boolean, and returns a hash with keys are the function and values are a reference to an array of the IDs that have that function.

If the case boolean is set the search is case insensitive. Otherwise case sensitive.

=cut

sub find_features_by_annotation {
	my($self,$anno_hash, $case)=@_;
	$self->load_functions;

	if ($case) {map {$anno_hash->{uc($_)}=1} keys %$anno_hash}
	
	my $res={};
	foreach my $id (keys %{$self->{_functions}})
	{
		my $fn = $self->{_functions}->{$id};
		$case ? $fn = uc($fn) : 1;
		if ($anno_hash->{$fn}) {push @{$res->{$fn}}, $id}
	}
	
	return $res;
}


=pod 

search_features_by_annotation

Takes a string to find and an optional case boolean, and returns a hash with keys are the function and values are a reference to an array of the IDs that have that function.

If the case boolean is set the search is case insensitive. Otherwise case sensitive.

Note that this was originally based on the find above, but this uses a regexp for matching. Will likely be a lot slower which is why I only search for a single term. There may be an SQL way of doing this, if so let Rob know how to do it and I'll replace this method.


=cut

sub search_features_by_annotation {
	my($self,$term, $case)=@_;
	$self->load_functions;

	# to make case insensitive convert everything to uppercase
	# alternative is to use two regexps, one for case insens and one for not case insens
	# but the bad thing about that approach is that if you have a case insensitive search
	# you do two regexps for each failed match
	
	$case ? $term = uc($term) : 1;

	my $res={};
	foreach my $id (keys %{$self->{_functions}})
	{
		# we set two variables, one that has case changed for case insensitive searches
		my $fn = my $fnc = $self->{_functions}->{$id};
		$case ? $fn = uc($fn) : 1;
		if ($fn =~ m/$term/) {push @{$res->{$fnc}}, $id}
	}
	
	return $res;
}


sub feature_aliases {
    my($self,$fid) = @_;

    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};

    my @aliases;

    if (($fid =~ /^fig\|(\d+\.\d+)/) && ($1 eq $newG))
    {
	&load_tbl($self);
	if (my $x = $self->{_tbl}->{$fid})
	{
	    @aliases = @{$x->[1]};
	}
	else
	{
	    @aliases = ();
	}
    }
    else
    {
	@aliases = $fig->feature_aliases($fid);
    }
    return wantarray() ? @aliases : join(",",@aliases);
}

sub get_corresponding_ids {
    my($self, $id, $with_type_info) = @_;
    
    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};

    if (($id =~ /^fig\|(\d+\.\d+)/) && ($1 eq $newG)) {
	my @aliases = $self->feature_aliases($id);
	my @corresponding_ids = ();
	foreach my $alias (@aliases) {
	    if ($alias =~ /^gi\|/) {
		if ($with_type_info) {
		    push(@corresponding_ids, [$alias, 'NCBI']);
		} else {
		    push(@corresponding_ids, $alias);
		}
		last;
	    }
	}
	return @corresponding_ids;
    } else {
	return $fig->get_corresponding_ids($id, $with_type_info);
    }

}

sub feature_annotations {
    my($self,$fid,$rawtime) = @_;

    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};

    my @annotations;
    if (($fid =~ /^fig\|(\d+\.\d+)/) && ($1 eq $newG))
    {
	&load_ann($self);
	if (my $x = $self->{_ann}->{$fid})
	{
	    @annotations = @{$x};
	}
	else
	{
	    @annotations = ();
	}

	if ($rawtime)
	{
	    return @annotations;
	}
	else
	{
	    return map { my $r = [@$_]; $r->[1] = localtime($r->[1]); $r } @annotations;
	}
    }
    else
    {
	return $fig->feature_annotations($fid);
    }
}

sub get_translation {
    my($self,$peg) = @_;

    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};

    if (($peg =~ /^fig\|(\d+\.\d+)/) && ($1 eq $newG))
    {
	&load_feature_indexes($self);

	my $out = $self->{_feat}->{peg}->{$peg};
	if (defined($out))
	{
	    return $out;
	}

	#
	# If we have a blast-formatted fasta, use fastacmd to
	# do the lookup, and cache the output for later use.
	#

	if ($self->{_feat_fasta}->{peg})
	{
	    my $id = "gnl|$peg";
	    my $cmd = "$FIG_Config::ext_bin/fastacmd -d $self->{_feat_fasta}->{peg} -s '$id'";
	    open(P, "$cmd|") or die "get_translation: cmd failed with $!: $cmd";
	    $_ = <P>;
	    my $out;
	    while (<P>)
	    {
		s/\s+//g;
		$out .= $_;
	    }
	    close(P);
	    $self->{_feat}->{$peg} = $out;
	    return $out;
	}
	else
	{
	    return $self->{_feat}->{$peg};
	}
    }
    else
    {
	return $fig->get_translation($peg);
    }
}

sub translation_length
{
    my($self, $peg) = @_;

    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};

    if (($peg =~ /^fig\|(\d+\.\d+)/) && ($1 eq $newG))
    {
	my $t = $self->get_translation($peg);
	return length($t);
    }
    else
    {
	return $fig->translation_length($peg);
    }
}

sub translatable
{
    my($self, $peg) = @_;

    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};

    if (($peg =~ /^fig\|(\d+\.\d+)/) && ($1 eq $newG))
    {
	return $self->translation_length($peg) > 0 ? 1 : 0;
    }
    else
    {
	return $fig->translatable($peg);
    }
}

sub pick_gene_boundaries {
    return &FIG::pick_gene_boundaries(@_);
}

sub call_start {
    return &FIG::call_start(@_);
}

sub is_real_feature
{
    my($self, $fid) = @_;

    if ($self->is_virtual_feature($fid) && ($fid =~ /^fig\|\d+\.\d+\.([^\.]+)/))
    {
	my $type = $1;
	my $typeH;
	$typeH = $self->{_features}->{$type};
	if (! $typeH)
	{
	    $typeH = &load_feature_hash($self,$type);
	    $self->{_features}->{$type} = $typeH;
	}
	return $typeH->{$fid} ? 1 : 0;
    }
    else
    {
	return $self->{_fig}->is_real_feature($fid);
    }

}

sub is_deleted_fid
{
    my($self, $fid) = @_;
    my $newG    = $self->{_genome};

    if (($fid =~ /^fig\|(\d+\.\d+)/) && ($1 eq $newG))
    {
	my $delH;
	$delH = $self->{_deleted_features};
	if (! $delH)
	{
	    $delH = &FIGV::load_deleted_features_hash($self);
	}
	return $delH->{$fid} ? 1 : 0;
    }
    else
    {
	return $self->{_fig}->is_deleted_fid($fid);
    }
}

sub load_deleted_features_hash {
    my($self) = @_;

    my $newGdir = $self->{_orgdir};
    my $deletedH = {};

    my @del_files = (<$newGdir/deleted.fids>, <$newGdir/Features/*/deleted.features>);

    for my $del (@del_files)
    {
	if (open(DEL, "<", $del))
	{
	    my $eol = $/;
	    $/ = "\n";
	    while (<DEL>)
	    {
		if ($_ =~ /^(\S+)/)
		{
		    $deletedH->{$1} = 1;
		}
	    }
	    close(DEL);
	    $/ = $eol;
	}
    }
    $self->{_deleted_features} = $deletedH;

    return $deletedH;
}

sub delete_feature {
    my($self,$user,$fid) = @_;
    my $fig     = $self->{_fig};
    
    my $genome = &FIG::genome_of($fid);
    my $newG    = $self->{_genome};
    if ($genome eq $newG) {
	my $newGdir = $self->{_orgdir};
	if (open(DEL,">>$newGdir/deleted.fids")) {
	    print DEL "$fid\t$user\n";
	    close(DEL);
	    &load_deleted_features_hash($self);
	}
	else {
	    carp "could not open $newGdir/deleted.fids: failed to delete $fid (user=$user)\n";
	}
    }
    else {
	return $fig->delete_feature($user,$fid);
    }
}

sub pegs_of {
    my($self, $genome) = @_;
    
    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};
    
    if ($genome ne $newG) {
	return $fig->pegs_of($genome);
    }

    $self->load_feature_indexes();

    my $t = $self->{_feat_btree}->{peg};
    
    if ($t) {
	return grep {/^fig\|/} keys %$t;
    }
    else {
	warn "Reverting to load_tbl for $newG\n";
	$self->load_tbl();
	
	return grep { /\.peg\./ } keys %{$self->{_tbl}};
    }
}


sub rnas_of {
    my($self, $genome) = @_;
    
    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};
    
    if ($genome ne $newG) {
	return $fig->rna_of($genome);
    }
    
    $self->load_feature_indexes();
    
    my $t = $self->{_feat_btree}->{rna};
    if ($t) {
	return grep {/^fig\|/} keys %$t;
    }
    else {
	warn "Reverting to load_tbl for $newG\n";
	$self->load_tbl();
	return grep { /\.rna\./ } keys %{$self->{_tbl}};
    }
}

sub is_virtual_feature {
    my($self, $fid) = @_;
    
    my $newG    = $self->{_genome};
    
    if (($fid =~ /^fig\|(\d+\.\d+)\.([^\.]+)/) && ($1 eq $newG)) {
	return 1;
    }
    else {
	return 0;
    }
}

sub load_feature_hash {
    my($self,$type) = @_;

    my $newGdir = $self->{_orgdir};
    my $typeH = {};
    if (open(FIDS,"<$newGdir/Features/$type/tbl"))
    {
	# cluck "load1 $newGdir\n";
	while (my $l = <FIDS>)
	{
	    if ($l =~ /^(\S+)/)
	    {
		my $fid = $1;
		if (! $self->is_deleted_fid($fid))
		{
		    $typeH->{$fid} = 1;
		}
	    }
	}
	close(FIDS);
    }
    return $typeH;
}
    

####################################################
#
# Following are some MG-RAST specific features. FIGV seems a good a place as any to include them.#
#
#

=head3 taxa_to_seed_ids

Given a prefix of a taxonomy string, return the list of metagenome fids that
mapped to SEED organisms in that taxonomy.

=cut

sub taxa_to_seed_ids
{
    my($self, $tax) = @_;

    my $btree_file = "$self->{_orgdir}/taxa_summary_by_blast.btree";
    my %btree;

    my $btree_tie = tie %btree, 'DB_File', $btree_file, O_RDONLY, 0666, $DB_BTREE;

    if (!$btree_tie)
    {
	warn "Cannot tie $btree_file: $!\n";
	return undef;
    }

    my $key = $tax;
    my ($res, $value);
    my @vals;
    my %seen;

    for ($res = $btree_tie->seq($key, $value, R_CURSOR);
	 $res == 0 and $key =~ /^$tax/;
	 $res = $btree_tie->seq($key, $value, R_NEXT))
    {
	print "$key\n";
	push(@vals, map { [$key, $_] } grep { not $seen{$_}++ } split(/$;/, $value));
    }
    return @vals;
}

sub load_feature_indexes
{
    my($self) = @_;

    if ($self->{_feat}) { return };

    my $newGdir = $self->{_orgdir};

    for my $fdir (<$newGdir/Features/*>)
    {
	my $ftype = basename($fdir);

	#
	# If we have a tbl.btree, tie that for our use.
	#
	
	my $tbl_idx = {};
	my $tie = tie %$tbl_idx, 'DB_File', "$fdir/tbl.btree", O_RDONLY, 0666, $DB_BTREE;
	if ($tie)
	{
	    $self->{_feat_tie}->{$ftype} = $tie;
	    $self->{_feat_btree}->{$ftype} = $tbl_idx;
	}

	my $tbl_list = [];
	my $ltie = tie @$tbl_list, 'DB_File', "$fdir/tbl.recno", O_RDONLY, 0666, $DB_RECNO;
	if ($tie)
	{
	    $self->{_feat_ltie}->{$ftype} = $ltie;
	    $self->{_feat_recno}->{$ftype} = $tbl_list;
	}

	#
	# If we have fasta.norm.phr, set _pseq_fasta to the fasta file to use with fastacmd.
	#
	
	my $pseq     = {};

	if (-f "$fdir/fasta.norm.phr")
	{
	    $self->{_feat_fasta}->{$ftype} = "$fdir/fasta.norm";

	}
	else
	{
	    #
	    # Otherwise, we need to load the data.
	    #

	    if (open(FASTA,"<$newGdir/Features/peg/fasta"))
	    {
		local $/ = "\n>";
		my $x;
		while (defined($x = <FASTA>))
		{
		    chomp $x;
		    if ($x =~ />?(\S+)[^\n]*\n(.*)/s)
		    {
			my $peg = $1;
			my $seq = $2;
			$seq =~ s/\s//gs;
			if (! $self->is_deleted_fid($peg))
			{
			    $pseq->{$peg} = $seq;
			}
		    }
		}
		close(FASTA);
	    }
	}
	$self->{_feat}->{$ftype} = $pseq;
    }
}

sub load_pseq {
    my($self) = @_;

    if ($self->{_pseq}) { return };

    my $newGdir = $self->{_orgdir};
    my $fdir = "$newGdir/Features/peg";

    #
    # If we have a tbl.btree, tie that for our use.
    #

    my $tbl_idx = {};
    my $tie = tie %$tbl_idx, 'DB_File', "$fdir/tbl.btree", O_RDONLY, 0666, $DB_BTREE;
    if ($tie)
    {
	$self->{_tbl_tie} = $tie;
	$self->{_tbl_btree} = $tbl_idx;
    }

    #
    # If we have fasta.norm.phr, set _pseq_fasta to the fasta file to use with fastacmd.
    #

    my $pseq     = {};

    if (-f "$fdir/fasta.norm.phr")
    {
	$self->{_pseq_fasta} = "$fdir/fasta.norm";
    }
    else
    {
	#
	# Otherwise, we need to load the data.
	#

	if (open(FASTA,"<$newGdir/Features/peg/fasta"))
	{
	    local $/ = "\n>";
	    my $x;
	    while (defined($x = <FASTA>))
	    {
		chomp $x;
		if ($x =~ />?(\S+)[^\n]*\n(.*)/s)
		{
		    my $peg = $1;
		    my $seq = $2;
		    $seq =~ s/\s//gs;
		    if (! $self->is_deleted_fid($peg))
		    {
			$pseq->{$peg} = $seq;
		    }
		}
	    }
	    close(FASTA);
	}
    }
    $self->{_pseq} = $pseq;
}

sub load_ss_data
{
    my($self, $force) = @_;
    
    return if defined($self->{_ss_bindings}) && (! $force);
    
    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};
    my $SS_dir  = "$newGdir/Subsystems";
    
    my $peg_index = {};
    my $bindings  = {};
    my $variant   = {};
    
    if (!-d $SS_dir) {
	print STDERR qq(WARNING: No directory \'$SS_dir\'\n) if $ENV{FIG_VERBOSE};
    }
    else {
	if (!-s "$SS_dir/bindings") {
	    print STDERR qq(WARNING: File \'$SS_dir/bindings\' does not exist or has zero size\n) if $ENV{FIG_VERBOSE};
	}
	else {
	    open(SS_FILE, "<$SS_dir/bindings")
		|| confess "Cannot read-open file=\'$SS_dir/bindings\': $!";
	    
	    while (<SS_FILE>) {
		chomp;
		my($sname, $role, $peg) = split(/\t/);
		if (! $self->is_deleted_fid($peg)) {
		    push(@{$bindings->{$sname}->{$role}}, $peg);
		    push(@{$peg_index->{$peg}}, [$sname, $role]);
		}
	    }
	    close(SS_FILE);
	}
	
	if (!-s "$SS_dir/subsystems") {
	    print STDERR qq(WARNING: File \'$SS_dir/subsystems\' does not exist or has zero size\n) if $ENV{FIG_VERBOSE};
	}
	else {
	    open(SS_FILE, "<$SS_dir/subsystems")
		|| confess "Cannot read-open file=\'$SS_dir/subsystems\': $!";
	    
	    while (<SS_FILE>) {
		chomp;
		my($sname, $v) = split(/\t/);
		$variant->{$sname} = $v;
	    }
	    close(SS_FILE);
	}
    }
    
    $self->{_ss_bindings}     = $bindings;
    $self->{_ss_variants}     = $variant;
    $self->{_ss_peg_index}    = $peg_index;
}

sub load_tbl {
    my($self) = @_;

    if ($self->{_tbl}) { return };

    my $newGdir = $self->{_orgdir};
    my $tbl     = {};

    foreach my $x (`cat $newGdir/Features/*/tbl`)
    {
	chomp $x;
	if ($x =~ /^(\S+)\t(\S+)(\t(\S.*\S))?/)
	{
	    my $fid = $1;
	    my $loc = [split(/,/,$2)];
	    my $aliases = $4 ? [split(/\t/,$4)] : [];
	    if (! $self->is_deleted_fid($fid))
	    {
		$tbl->{$fid} = [$loc,$aliases];
	    }
	}
	else {
	    warn "Bad feature line in $newGdir:$x:\n";
	}
    }
    print STDERR ("Loaded ", (scalar keys %$tbl), " features from $newGdir\n") if $ENV{FIG_VERBOSE};
    print STDERR Dumper($tbl) if ($ENV{FIG_VERBOSE} && ($ENV{FIG_VERBOSE} > 1));
    
    $self->{_tbl} = $tbl;
}

use Carp qw(cluck);

sub load_functions {
    my($self, $force) = @_;

    my @fns;

    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};

    my @check_fns;
    #
    # check for functions newer than our timestamp
    #
    if (-f "$newGdir/assigned_functions" && -f "$newGdir/RAST")
    {
	@fns = ("$newGdir/assigned_functions");
	@check_fns = @fns;
    }
    else
    {
	# order of "cat" is important - proposed_user_functions must be last
	# @fns = <$newGdir/*functions>;
	@fns = map { "$newGdir/$_" } qw(assigned_functions
					proposed_functions
					proposed_non_ff_functions
					proposed_user_functions);
	@check_fns = map { "$newGdir/$_" } qw(proposed_user_functions);
    }
    
    if (!$force && $self->{_functions})
    {
	my $cached_mod = $self->{_function_date};
	my $need_redo = 0;
	for my $fn (@check_fns)
	{
	    my $fnd = (stat($fn))[9];
	    if (defined($fnd) && $fnd > $cached_mod)
	    {
		$need_redo = 1;
		print STDERR "need_redo: cached=$cached_mod fnd=$fnd fn=$fn\n";
		last;
	    }
	}
	return unless $need_redo;
    }

    my $functions = {};
    my $roles     = {};


    #
    # if a file RAST exists, this is a RAST-processed genome taht has been
    # installed in a seed. Use assigned_functions, not the concatenation.
    #
    my $latest_mod;
    for my $file (@fns)
    {
	if (open(my $fh, "<", $file))
	{
	    # print STDERR "read $file\n";
	    my $m = (stat($fh))[9];
	    $latest_mod = $m if ($m > $latest_mod);
	    while (my $x = <$fh>)
	    {
		chomp $x;
		# print STDERR "x=$x\n";
		if (($x =~ /^(fig\|(\d+\.\d+)\.\S+)\t(\S[^\t]*\S)/) && ($2 eq $newG))
		{
		    my $peg = $1;
		    my $f = $3;

		    # print STDERR "'$peg' '$f'\n";
		    if (! $self->is_deleted_fid($peg))
		    {
			# user has overridden a function, so delete old roles
			if (exists($functions->{$peg})) {
			    my @roles = &FIG::roles_of_function($functions->{$peg});
			    foreach $_ (@roles)
			    {
				delete $roles->{$_};
			    }
			}
			
			# add new roles
			my @roles = &FIG::roles_of_function($f);
			foreach $_ (@roles)
			{
			    push(@{$roles->{$_}},$peg);
			}
			
			# set function
			$functions->{$peg} = $f;
		    }
		    else
		    {
			# print STDERR "Deleted $peg\n";
		    }
		}
		else
		{
		    # print STDERR "nomatch '$x' '$newG'\n";
		}
	    }
	    close($fh);
	}
	else
	{
	    # warn "Cannot open $file: $!\n";
	}
    }
    $self->{_functions} = $functions;
    $self->{_function_date} = $latest_mod;
    $self->{_roles} = $roles;
#   print STDERR "cache: set latest=$latest_mod\n";
}

sub seqs_with_role {
    my($self,$role,$who,$genome) = @_;

    my $newG    = $self->{_genome};
    my $fig     = $self->{_fig};
    if ($genome && $newG && ($genome eq $newG)) {
	&load_functions($self);
    	my $pegsL = $self->{_roles}->{$role};	
	return ($pegsL ? @{$pegsL} : ());
    } else {
	my @result = $fig->seqs_with_role($role,$who,$genome);
	if ($newG) {
	    &load_functions($self);
	    my $pegsL = $self->{_roles}->{$role};
	    push(@result, ($pegsL ? @{$pegsL} : ()));
	}
	return @result;
    }
}
 
sub bbhs
{
    my($self,$peg,$cutoff) = @_;
    if ($self->is_deleted_fid($peg)) { return () }
    
    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};

    $self->load_bbh_index();

    $cutoff = 1.0e-10 unless defined($cutoff);

    my $tie = $self->{_bbh_tie};

    my @bbhs = $tie->get_dup($peg);

    # @bbhs now a list of comma-separated tuples (id2, score, bitscore)
#    print Dumper(\@bbhs);

    @bbhs = grep { $_->[1] < $cutoff } map { [split(/,/)] } @bbhs;

    if (not ($peg =~ /^fig\|(\d+\.\d+)/ && ($1 eq $newG)))
    {
	#
	# If this isn't one of our pegs, we retrieve the bbhs from the underlying
	# SEED and merge in the bbhs that we have here.
	#

	my @fbbhs = $fig->bbhs($peg, $cutoff);
#	print Dumper(\@fbbhs);
	push(@bbhs, @fbbhs);
    }
    return sort { $a->[1] <=> $b->[1] } @bbhs;
}

sub sims
{
    my($self,$pegarg,$max,$maxP,$select, $max_expand, $filters) = @_;

    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};
    $max     = $max ? $max : 10000;
    $maxP    = $maxP ? $maxP : 1.0e-5;
    $select     = $select ? $select : "all";

    my $prefix = "_sims";
    my $flip_prefix = "_flips";

    if ($select eq 'raw')
    {
	$prefix = "_raw_sims";
	$flip_prefix = "_raw_flips";
    }

    #
    # Partition pegs into one set that is part of this organism
    # and another set that is not.
    #
    my @pegs = ref($pegarg) ? @$pegarg : ($pegarg);
    my(@part, @not_part, %part);
    my %sims;

    for my $peg  (@pegs)
    {
	if ($peg =~ /^fig\|(\d+\.\d+)/ and $1 eq $newG)
	{
	    push(@part, $peg);
	    $part{$peg}++;
	}
	else
	{
	    push(@not_part, $peg);
	    $sims{$peg} = [];
	}
    }

    #
    # Retrieve a batch of the sims from the normal SEED sims for the not-part pegs, and partition
    # into %sims.
    #

    if (@not_part)
    {
	for my $sim ($fig->sims(\@not_part, $max, $maxP, $select, $max_expand, $filters))
	{
	    push(@{$sims{$sim->id1}}, $sim);
	}
    }

    my @out;
    my $start_len;

    for my $peg (@pegs)
    {
	$start_len = @out;

	if (not $part{$peg})
	{
	    #
	    # For the ids that are not part of this organism, pull the flips that we
	    # computed - these are the sims between this org and the one in the requested
	    # peg.
	    #
	    my @flips = $self->retrieve_sims($peg, 	$flip_prefix, $maxP, $select);

	    @flips = sort { $b->bsc <=> $a->bsc } @flips;

	    # my @old = $fig->sims($peg,$max,$maxP,$select);

	    #
	    # Merge these sims together with any that we retrieved
	    # earlier - recall sims are returned in bsc order.
	    #

	    my @old = sort { $b->bsc <=> $a->bsc } @{$sims{$peg}};

	    # my @merged = ();
	    my $i1 = 0;
	    my $i2 = 0;
	    while (($i1 < @flips) || ($i2 < @old))
	    {
		if (($i1 == @flips) || (($i2 < @old) && ($flips[$i1]->[11] < $old[$i2]->[11])))
		{
		    # push(@merged,$old[$i2]);
		    push(@out,$old[$i2]);
		    $i2++;
		}
		else
		{
		    # push(@merged,$flips[$i1]);
		    push(@out,$flips[$i1]);
		    $i1++;
		}
	    }
	}
	else
	{
	    #
	    # For the ids that are in this organism, we just pull the sims directly.
	    #
	    my @sims = $self->retrieve_sims($peg, $prefix, $maxP, $select);
	    push(@out, @sims);
	}

	if (@out - $start_len > $max)
	{
	    $#out = $start_len + $max - 1;
	}
    }

    return @out;
}

sub sims_for_figm
{
    my($self, $peer_genomes, $pegarg,$max,$maxP,$select, $max_expand, $filters) = @_;

    my $fig     = $self->{_fig};
    my $newGdir = $self->{_orgdir};
    my $peer_cache = $self->{_peer_sims_cache};

    my @out = $self->sims($pegarg,$max,$maxP,$select, $max_expand, $filters);

    for my $peer_genome (@$peer_genomes)
    {
	# print STDERR "Get pegs for peer $peer_genome\n";
	my $ssims = $peer_cache->{$peer_genome};
	if (!defined($ssims))
	{
	    $ssims = StandaloneSimsSet->new("$newGdir/sims",
					    fwd_raw => $peer_genome,
					    fwd_exp => $peer_genome);
	    $peer_cache->{$peer_genome} = $ssims;
	}
	for my $peg (ref($pegarg) ? @$pegarg : ($pegarg))
	{
	    my @psims = $ssims->sims_fwd($peg, $max, $maxP, $select, $max_expand, $filters);

	    # print STDERR "Sims for $peg: ", join("\n\t", @psims), "\n";
	    unshift(@out, @psims);
	}
    }

    return @out;
}

sub retrieve_sims
{
    my($self, $peg, $prefix, $maxP, $select) = @_;

    if ($self->is_deleted_fid($peg)) { return () }

    $self->load_sims_index();

    my $info = $self->{"${prefix}_index"}->{$peg};

    if ($info !~ /^(\d+),(\d+)$/)
    {
	return;
    }
    my($seek, $len) = ($1, $2);

    my $sims_txt = &FIG::read_block($self->{"${prefix}_fh"}, $seek, $len);

    my @sims = map { bless $_, 'Sim' } grep { $_->[10] <= $maxP } map { [split(/\t/)] } @$sims_txt;

    if ($select =~ /^figx?$/)
    {
	@sims = grep { $_->[1] =~ /^fig/ } @sims;
    }

    return @sims;
}

sub sims_old {
    my($self,$peg,$max,$maxP,$select) = @_;

    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};
    $max     = $max ? $max : 10000;
    $maxP    = $maxP ? $maxP : 1.0e-5;
    $select     = $select ? $select : "all";

    if (($peg =~ /^fig\|(\d+\.\d+)/) && ($1 eq $newG) && (! $self->is_deleted_fid($peg)))
    {
	&load_sims($self);
	my @sims = ();
	my $raw_sims = $self->{_sims}->{$peg};
	if ($raw_sims)
	{
	    foreach my $sim ( grep { $_->[10] <= $maxP } @$raw_sims )
	    {
		my $id2 = $sim->id2;
		my $id1 = $sim->id1;
		my @relevant = ();

		my @maps_to = $fig->mapped_prot_ids( $id2 );
		my $ref_len = $maps_to[0]->[1];

		@maps_to = grep { $_->[0] !~ /^xxx\d+/ } @maps_to;

		if ( $select =~ /^figx?$/ )          # Only fig
		{
		    @relevant = grep { $_->[0] =~ /^fig/ } @maps_to;
		}
		else                                 # All
		{
		    @relevant = @maps_to;
		}

		my $seen = {};
		foreach my $x ( @relevant )
		{
		    my ( $x_id, $x_ln ) = @$x;

		    next if $seen->{$x_id};
		    $seen->{$x_id} = 1;

		    my $delta2  = $ref_len - $x_ln;   # Coordinate shift
		    my $sim1    = [ @$sim ];                  # Make a copy
		    $sim1->[1]  = $x_id;
		    $sim1->[8] -= $delta2;
		    $sim1->[9] -= $delta2;
		    bless( $sim1, "Sim" );
		    push( @sims, $sim1 );
		}
	    }
	}

	if (@sims > $max) { $#sims = $max-1; }
	return @sims;
    }
    else
    {
	return $fig->sims($peg,$max,$maxP,$select);
    }
}


sub coupled_to
{
    my($self,$peg) = @_;

    if ($self->is_deleted_fid($peg)) { return () }

    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};

    if ($peg =~ /^fig\|$newG\.peg/)
    {

	$self->load_coupling_index();

	my $tie = $self->{_pch_tie};

	if (not defined($tie))
	{
	    return;
	}

	my @dat = $tie->get_dup($peg);

	return map { [split(/$;/, $_)] } @dat;
    }
    else
    {
	return $fig->coupled_to($peg);
    }

}

sub coupled_to_batch
{
    my($self,@pegs) = @_;

    return () unless @pegs;

    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};

    # divide pegs into figv and fig pegs
    my $mypegs = [];
    my $otherpegs = [];
    foreach my $peg (@pegs) {
	if ($peg =~ /^fig\|$newG\.peg/) {
	    push(@$mypegs, $peg);
	} else {
	    push(@$otherpegs, $peg);
	}
    }

    my $ret = [];
    if (scalar(@$mypegs)) {
	$self->load_coupling_index();
	
	my $tie = $self->{_pch_tie};
	
	if (defined($tie))
	{
	    foreach my $peg (@$mypegs) {
		my @dat = $tie->get_dup($peg);
		push(@$ret, map { [$peg, split(/$;/, $_)] } @dat);
	    }
	}
    }
    if (scalar(@$otherpegs)) {
	push(@$ret, $fig->coupled_to_batch(@$otherpegs));
    }

    return @$ret;
}

sub coupling_evidence
{
    my($self,$peg1, $peg2) = @_;

    if ($self->is_deleted_fid($peg1) || $self->is_deleted_fid($peg2)) { return () }

    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};

    if ($peg1 =~ /^fig\|$newG\.peg/)
    {
	$self->load_coupling_index();

	my $tie = $self->{_pch_ev_tie};

	if (not defined($tie))
	{
	    return;
	}

	my @dat = $tie->get_dup("$peg1$;$peg2");

	my @a;
	return map { @a = split(/$;/, $_); [@a[0,1,4]] } @dat;
    }
    else
    {
	return $fig->coupling_evidence($peg1, $peg2);
    }

}

sub coupling_and_evidence
{
    my($self,$peg1) = @_;

    if ($self->is_deleted_fid($peg1)) { return () }
    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};

    if ($peg1 =~ /^fig\|$newG\.peg/)
    {
	$self->load_coupling_index();

	my $tie = $self->{_pch_tie};
	my $evtie = $self->{_pch_ev_tie};

	if (not(defined($tie) and defined($evtie)))
	{
	    return;
	}

	my @out;
	my @coupled_to = $tie->get_dup($peg1);
	for my $ent (@coupled_to)
	{
	    my ($peg2, $score) = split(/$;/, $ent);

	    my @ev = $evtie->get_dup("$peg1$;$peg2");

	    my $l = [];
	    for my $event (@ev)
	    {
		my($peg3, $peg4, $iden3, $iden4, $rep) = split(/$;/, $event);
		push(@$l, [$peg3, $peg4]);
	    }
	    push(@out, [$score, $peg2, $l]);
	}
	return @out;
    }
    else
    {
	return $fig->coupling_and_evidence($peg1);
    }

}

sub in_pch_pin_with
{
    my($self, $peg1, $diverse) = @_;

    my @all = $self->in_pch_pin_with_and_evidence($peg1);

    if ($diverse)
    {
	return map { $_->[0] } grep { $_->[1] == 1 } @all;
    }
    else
    {
	return map { $_->[0] } @all;
    }
}

sub in_pch_pin_with_and_evidence
{
    my($self,$peg1) = @_;

    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};

    if ($peg1 !~ /^fig\|$newG\.peg/)
    {
	return $fig->in_pch_pin_with_and_evidence($peg1);
    }
    if ($self->is_deleted_fid($peg1)) { return () }

    $self->load_coupling_index();

    my $evtie = $self->{_pch_ev_tie};

    my($key, $value, $st);

    my %max;

    $key = "$peg1$;";
    my $qkey = quotemeta($key);
    for ($st = $evtie->seq($key, $value, R_CURSOR); $st == 0 and $key =~ /^$qkey/; $st = $evtie->seq($key, $value, R_NEXT))
    {
#	print "key=$key value=$value\n";

	my($peg3, $peg4, $iden3, $iden4, $rep) = split(/$;/, $value);

	if (exists($max{$peg3}))
	{
	    $max{$peg3} = $rep if $rep > $max{$peg3};
	}
	else
	{
	    $max{$peg3} = $rep;
	}
    }

    return map { [$_, $max{$_}] } keys %max;
}


sub load_coupling_index
{
    my($self) = @_;

    return if defined($self->{_pch_index});

    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};

    my $pch_btree = "$newGdir/pchs.btree";
    my $ev_btree = "$newGdir/pchs.evidence.btree";

    my $pch_index = {};
    my $ev_index = {};

    my $tied = tie %$pch_index, 'DB_File', $pch_btree, O_RDONLY, 0666, $DB_BTREE;
    my $evtied = tie %$ev_index, 'DB_File', $ev_btree, O_RDONLY, 0666, $DB_BTREE;

    #
    # Set these even if failed so we don't keep trying to open and failing.
    #
    $self->{_pch_index} = $pch_index;
    $self->{_pch_tie} = $tied;
    $self->{_pch_ev_index} = $ev_index;
    $self->{_pch_ev_tie} = $evtied;

    if (not $tied)
    {
	warn "Cannot tie pch index $pch_btree: $!\n";
    }
}

sub load_bbh_index
{
    my($self) = @_;

    return if defined($self->{_bbh_index});

    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};

    my $bbh_file = "$newGdir/bbhs";
    my $bbh_index_file = "$bbh_file.index";
    my $bbh_index = {};

    my $tied = tie %$bbh_index, 'DB_File', $bbh_index_file, O_RDONLY, 0666, $DB_BTREE;

    #
    # Set these even if failed so we don't keep trying to open and failing.
    #
    $self->{_bbh_index} = $bbh_index;
    $self->{_bbh_tie} = $tied;

    if (not $tied)
    {
	warn "Cannot tie bbh index $bbh_index_file: $!\n";
    }
}

sub load_contigs_index
{
    my($self) = @_;

    return if defined($self->{_contigs_index});

    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};

    my $contig_index_file = "$newGdir/contigs.btree";
    my $contig_index = {};
    my $contig_len_index_file = "$newGdir/contig_len.btree";
    my $contig_len_index = {};

    my $tied = tie %$contig_index, 'DB_File', $contig_index_file, O_RDONLY, 0666, $DB_BTREE;
    if (not $tied)
    {
	# warn "Cannot tie contig index $contig_index_file: $!\n";
    }

    my $ltied = tie %$contig_len_index, 'DB_File', $contig_len_index_file, O_RDONLY, 0666, $DB_BTREE;
    if (not $ltied)
    {
	# warn "Cannot tie contig length index $contig_len_index_file: $!\n";
    }

    #
    # Set these even if failed so we don't keep trying to open and failing.
    #
    $self->{_contigs_index} = $contig_index;
    $self->{_contigs_tie} = $tied;
    $self->{_contig_len_index} = $contig_len_index;
    $self->{_contig_len_tie} = $ltied;

}

sub load_sims_index
{
    my($self) = @_;

    return if defined($self->{_sims_index});

    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};

    my $sims_file = "$newGdir/expanded_similarities";
    my $sims_index_file = "$sims_file.index";
    my $sims_index = {};

    my $flips_file = "$newGdir/expanded_similarities.flips";
    my $flips_index_file = "$flips_file.index";
    my $flips_index = {};

    my $raw_sims_file = "$newGdir/similarities";
    my $raw_sims_index_file = "$raw_sims_file.index";
    my $raw_sims_index = {};

    my $raw_flips_file = "$newGdir/similarities.flips";
    my $raw_flips_index_file = "$raw_flips_file.index";
    my $raw_flips_index = {};

    $self->open_sims($sims_index, $sims_file, $sims_index_file, "_sims");
    $self->open_sims($flips_index, $flips_file, $flips_index_file, "_flips");
    $self->open_sims($raw_sims_index, $raw_sims_file, $raw_sims_index_file, "_raw_sims");
    $self->open_sims($raw_flips_index, $raw_flips_file, $raw_flips_index_file, "_raw_flips");
}

#
# Open a sims file, tie it to a hash, and store the info into the $self obj.
#
sub open_sims
{
    my($self, $hash, $sims_file, $index_file, $prefix) = @_;

    my $tied = tie %$hash, 'DB_File', $index_file, O_RDONLY, 0666, $DB_BTREE;

    #
    # Set these even if failed so we don't keep trying to open and failing.
    #
    $self->{"${prefix}_index"} = $hash;
    $self->{"${prefix}_tie"} = $tied;

    if (not $tied)
    {
	warn "Cannot tie sims index $index_file: $!\n";
    }

    #
    # open the sims file as well.
    #

    $self->{"${prefix}_fh"} = new FileHandle("<$sims_file");

    if (!$self->{"${prefix}_fh"})
    {
	warn "Cannot open sims file $sims_file: $!\n";
    }
}

sub load_sims {
    my($self) = @_;

    if ($self->{_sims}) { return };

    my $newGdir = $self->{_orgdir};

    my $sims     = {};
    foreach my $x (`cat $newGdir/similarities`)
    {
	chomp $x;
	if ($x =~ /^(\S+)/)
	{
	    my $fid = $1;
	    if (! $self->is_deleted_fid($fid))
	    {
		push(@{$sims->{$1}},bless([split(/\t/,$x)],'Sim'));
	    }
	}
    }
    $self->{_sims} = $sims;
}

sub get_attributes {
    my($self, $id, $attr) = @_;

    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};

    if ((($id =~ /^fig\|(\d+\.\d+)\.peg\.\d+$/) and ($1 eq $newG)) or
	(($id =~ /^(\d+\.\d+)/  and $1 eq $newG)))
    {
	&load_attr($self);

	my $t = $self->{_attr_id_tie};
	if (!$t)
	{
	    #
	    # Tied index not present, bail back to simple attributes.
	    #

	    my $l = $self->{_attr_id}->{$id};
	    return $l ? @$l : ();
	}

	my @v = $t->get_dup($id);

	my @out;
	for my $v (@v)
	{
	    my @a = split(/$;/, $v);

	    if (!defined($attr) or $attr eq $a[1])
	    {
		push(@out, [@a]);
	    }
	}
	return @out;
    }
    else
    {
	my @f = $fig->get_attributes($id, $attr);
	if (!defined($id) and defined(my $t = $self->{_attr_key_tie}))
	{
	    #
	    # lookup locally for $attr matches and merge with other output.
	    #

	    my @mine = $t->get_dup($attr);
	    @mine = map { [split(/$;/, $_)] } @mine;

	    return (@mine, @f);
	}
	else
	{
	    return @f;
	}
    }
}

sub load_attr {
    my($self) = @_;

    if ($self->{_attr_id}) { return };

    my $newGdir = $self->{_orgdir};

    my $id = {};
    my $key = {};

    $self->{_attr_id_tie} = tie %$id, 'DB_File', "$newGdir/attr_id.btree", O_RDONLY, 0, $DB_BTREE;
    $self->{_attr_key_tie} = tie %$key, 'DB_File', "$newGdir/attr_key.btree", O_RDONLY, 0, $DB_BTREE;

    $self->{_attr_id} = $id;
    $self->{_attr_key} = $key;

    #
    # If the tie didn't work for ids, at least load up the evidence codes.
    #

    if (!$self->{_attr_id_tie})
    {
	foreach my $x (`cat $newGdir/evidence.codes`)
	{
	    if ($x =~ /^(\S+)\t(\S+)/)
	    {
		push(@{$id->{$1}},[$1,"evidence_code",$2,""]);
	    }
	}
    }

}

sub load_ann {
    my($self) = @_;

    if ($self->{_ann}) { return };

    my $newGdir = $self->{_orgdir};
    my $ann     = {};
    if (open(ANN,"<$newGdir/annotations"))
    {
	$/ = "\n//\n";
	while (defined(my $x = <ANN>))
	{
	    chomp $x;
	    if ($x =~ /^(\S+)\n([^\n]+)\n([^\n]+)\n(.*)/s)
	    {
		push(@{$ann->{$1}},[$1,$2,$3,"$4\n"]);
	    }
	}
	$/ = "\n";
	close(ANN);
    }
    $self->{_ann} = $ann;
}

sub taxonomy_of {
    my($self,$genome) = @_;
    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};

    if ($newG ne $genome)
    {
	return $fig->taxonomy_of($genome);
    }

    my $tax;
    if (open(TAX,"<$newGdir/TAXONOMY") && ($tax = <TAX>))
    {
	chop $tax;
	return $tax;
    }
    else
    {
	return "unknown";
    }
}

sub build_tree_of_complete {
    my($self,$min_for_label) = @_;
    return $self->build_tree_of_all($min_for_label, "complete");
}

sub build_tree_of_all {
    my($self, $min_for_label, $complete)=@_;
    my(@last,@tax,$i,$prefix,$lev,$genome,$tax);

    $min_for_label = $min_for_label ? $min_for_label : 10;
    open(TMP,">/tmp/tree$$") || die "could not open /tmp/tree$$";
    print TMP "1. root\n";

    @last = ();


    foreach $genome (grep { ! $self->is_environmental($_) } $self->sort_genomes_by_taxonomy($self->genomes($complete)))
    {
        $tax = $self->taxonomy_of($genome);
        @tax = split(/\s*;\s*/,$tax);
        push(@tax,$genome);
        for ($i=0; ((@last > $i) && (@tax > $i) && ($last[$i] eq $tax[$i])); $i++) {}
        while ($i < @tax)
        {
            $lev = $i+2;
            $prefix = " " x (4 * ($lev-1));
            print TMP "$prefix$lev\. $tax[$i]\n";
            $i++;
        }
        @last = @tax;
    }
    close(TMP);
    my $tree = &tree_utilities::build_tree_from_outline("/tmp/tree$$");
    $tree->[0] = 'All';
    &FIG::limit_labels($tree,$min_for_label);
    unlink("/tmp/tree$$");
    return ($tree,&tree_utilities::tips_of_tree($tree));
}

sub sort_genomes_by_taxonomy {
    my($self,@genomes) = @_;

    return map     { $_->[0] }
           sort    { $a->[1] cmp $b->[1] }
           map     { [$_,$self->taxonomy_of($_)] }
           @genomes;
}

sub taxonomic_groups_of_complete {
    my($self,$min_for_labels) = @_;

    my($tree,undef) = $self->build_tree_of_complete($min_for_labels);
    return &FIG::taxonomic_groups($tree);
}

=head2 Search Database

Searches the database for objects that match the query string in some way.

Returns a list of results if the query is ambiguous or an unique identifier
otherwise.

=cut

sub search_database {
    # get parameters
    my ($self, $query, $options) = @_;

    my $fig     = $self->{_fig};
    my $newG    = $self->{_genome};
    my $newGdir = $self->{_orgdir};

    #
    # Check for local ids and genome.
    #


    if ($query eq $newG or lc($query) eq lc($self->genus_species($newG)))
    {
	return { type => 'organism', result => $newG };
    }

    if ($query =~ /^fig\|(\d+\.\d+)\.peg/ and $1 eq $newG)
    {
	if ($self->deleted_fid($query)) { return {} }
	return { type => 'feature', result => $query };
    }

    #
    # Match on functions
    #

    &load_functions($self);

    my @protlist;
    my $fn = $self->{_functions};
    for my $id (keys %$fn)
    {
	if ($fn->{$id} =~ /$query/i)
	{
	    push @protlist, [$id, $fn->{$id}, $newG];
	}
    }

    #
    # Pull FIG lookups.
    #

    my $res = $fig->search_database($query, $options);

    #
    # And integrate our protein list
    #

    my $res_prot;
    if (ref($res) eq 'ARRAY')
    {
	for my $ent (@$res)
	{
	    if ($ent->{type} eq 'proteins')
	    {
		$res_prot = $ent->{result};
	    }
	}
	if (!$res_prot)
	{
	    $res_prot = {type => 'proteins', result => \@protlist}
	}
	else
	{
	    push(@$res_prot, @protlist);
	}
    }

    return $res;

}


=head3 scenario_directory

    FIG->scenario_directory($organism);

Returns the scenario directory of an organism.  If the organism is 'All', returns
the directory containing all possible paths through scenarios.

=over 4

=item $organism

The seed-taxonomy id of the organism, e.g. 83333.1, or 'All'.

=back

=cut

sub scenario_directory {
  my ($self, $organism) = @_;

  my $directory;

  if ($organism eq $self->{_genome})
  {
      $directory = $self->{_orgdir} . "/Scenarios";
  }
  elsif(!defined $organism && defined $self->{_genome})
  {
      $directory = $self->{_orgdir} . "/Scenarios";
  }
  else {
      $directory = $self->{_fig}->scenario_directory($organism);
  }

  return $directory;
}

sub update_bindings_and_variant {
    my ($self, $ssname, $bindings, $variant) = @_;

    # get the genome directory
    my $newGdir = $self->{_orgdir};

    # read existing bindings for all other subsystems
    my @old_bindings;
    open(S, "<$newGdir/Subsystems/bindings") or die "Cannot open $newGdir/Subsystems/bindings: $!";

    while (<S>) {
	chomp;
	my($sname, $role, $peg) = split(/\t/);
	if ($sname ne $ssname) {
	    push @old_bindings, "$sname\t$role\t$peg\n";
	}
    }
    close(S);

    # mdj: changed this to always rewrite bindings, because the logic for determining whether
    # mdj: there are changes was faulty and probably more time consuming than rewriting the file
    # mdj: seems like there should be a lock condition on reading/writing files here
    open(SS, ">$newGdir/Subsystems/bindings~") or die "Cannot open $newGdir/Subsystems/bindings~: $!";
    map { print SS $_ } @old_bindings;
    foreach my $role (keys(%{$bindings})) {
	foreach my $peg (@{$bindings->{$role}}) {
	    print SS "$ssname\t$role\t$peg\n";
	}
    }
    close(SS);
    rename("$newGdir/Subsystems/bindings~", "$newGdir/Subsystems/bindings");

    # update variant
    open(S, "<$newGdir/Subsystems/subsystems") or die "Cannot open $newGdir/Subsystems/subsystems: $!";
    open(SS, ">$newGdir/Subsystems/subsystems~") or die "Cannot open $newGdir/Subsystems/subsystems~: $!";
    while (<S>) {
	chomp;
	my($sname, $v) = split(/\t/);
	unless ($sname eq $ssname) {
	    print SS "$sname\t$v\n";
	}
    }
    print SS "$ssname\t$variant\n";
    close(SS);
    close(S);
    rename("$newGdir/Subsystems/subsystems~", "$newGdir/Subsystems/subsystems");
    $self->load_ss_data(1); # force a reload
}

sub find_role_in_org {
    my ($self, $role, $org, $user, $cutoff) = @_;

    my($id2,$psc,$col_hdrs,$tab,$peg,$curr_func,$id2_func,$seen);

    if (!$org)
    {
        return undef;
    }

    #
    # Create a list of candidates.
    #
    # These are the list of sequences that contain the given role,
    # sorted by the crude_estimate_of_distance from the given peg.
    #

    my @cand = map { $_->[0] }
               sort { $a->[1] <=> $b->[1] }
               map {
                      $peg = $_;
                      [$peg,$self->crude_estimate_of_distance($org,&FIG::genome_of($peg))]
                   }
               $self->seqs_with_role($role,$user);

    my $hits = {};
    $seen = {};

    #
    # Pick the top 10 hits if there are more than 10.
    #
    my $how_many0 = ((@cand > 10) ? 10 : scalar @cand) - 1;

    $self->try_to_locate($org,$hits,[@cand[0..$how_many0]],$seen, $cutoff);

    if (keys(%$hits) == 0)
    {
        splice(@cand,0,$how_many0+1);
        &try_to_locate($self,$org,$hits,\@cand,$seen, $cutoff);
    }

    #
    # At this point %$hits contains the pegs in our organism that
    # may have the given role. The key is the peg, the value
    # is a pair [score, similar-peg]
    #
    #
    # We reformat this into a list of entries
    # [ $psc, $peg-in-this-org, $length, $current-fn, $matched-protein, $matched-len, $matched-fun]
    #


    $col_hdrs = ["P-Sc","PEG","Ln1","Current Function", "Protein Hit","Ln2","Function"];

    my @ret;

    foreach $peg ( sort {$hits->{$a}->[0] <=> $hits->{$b}->[0]} keys(%$hits))
    {
        ($psc,$id2) = @{$hits->{$peg}};
        $curr_func = $self->function_of($peg,$user);
        $id2_func  = $self->function_of($id2,$user);

        push(@ret, [$psc, $peg, $self->translation_length($peg),
                    $curr_func, $id2, $self->translation_length($id2),$id2_func]);
    }
    return @ret;
}

sub try_to_locate {
    my($self,$genome,$hits,$to_try,$seen,$cutoff) = @_;
    my($prot,$id2,$psc,$id2a,$x,$sim);
    if (! $cutoff) { $cutoff = 1.0e-5 }

    foreach $prot (@$to_try)
    {
        if (! $seen->{$prot})
        {
            if (($prot =~ /^fig\|(\d+\.\d+)/) && ($1 eq $genome))
            {
                $hits->{$prot} = [0,$prot];
            }
            else
            {
                foreach $sim ($self->sims($prot,1000,$cutoff,"fig"))
                {
                    $id2 = $sim->id2;
		    next if $self->is_deleted_fid($id2);
                    $psc = $sim->psc;
                    if (($id2 =~ /^fig\|(\d+\.\d+)/) && ($1 eq $genome))
                    {
                        $x = $hits->{$id2};
                        if ((! $x) || ($x->[0] > $psc))
                        {
                            $hits->{$id2} = [$psc,$prot];
                        }
                    }
                    elsif (&neg_log($psc) > (2 * &neg_log($cutoff)))
                    {
                        $seen->{$id2} = 1;
                    }
                }
            }
        }
    }
}

sub neg_log {
    my($x) = @_;

    if ($x == 0)
    {
        return 200;
    }
    else
    {
        return -log($x) / log(10);
    }
}

sub crude_estimate_of_distance {
    my($self,$genome1,$genome2) = @_;
    my($i,$v,$d,$dist);

    if ($genome1 > $genome2) { ($genome1,$genome2) = ($genome2,$genome1) }
    $dist = $self->cached('_dist');
    if (! $dist->{"$genome1,$genome2"})
    {
        my @tax1 = split(/\s*;\s*/,$self->taxonomy_of($genome1));
        my @tax2 = split(/\s*;\s*/,$self->taxonomy_of($genome2));

        $d = 1;
        for ($i=0, $v=0.5; ($i < @tax1) && ($i < @tax2) && ($tax1[$i] eq $tax2[$i]); $i++, $v = $v/2)
        {
            $d -= $v;
        }
        $dist->{"$genome1,$genome2"} = $d;
    }
    return $dist->{"$genome1,$genome2"};
}

sub check_db_peg_to_fams {
    my ($self, $peg_to_fams_hash) = @_;

    # get the genome directory
    my $newGdir = $self->{_orgdir};
    
    if (open(FH, "$newGdir/found")) {
	while (<FH>) {
	    chomp;
	    my ($id, $figfam, $function) = split /\t/;
	    $peg_to_fams_hash->{$id} = $figfam;
	}
	close FH;
    }
    
    return $peg_to_fams_hash;
}

sub check_db_genome_to_fams {
    my ($self, $genome_to_fams_hash) = @_;
    
    # get the genome directory
    my $newGdir = $self->{_orgdir};
    
    my @all_fams;
    if (open(FH, "$newGdir/found")) {
	while (<FH>) {
	    chomp;
	    my ($id, $figfam, $function) = split /\t/;
	    push(@all_fams, $figfam);
	}
	close FH;
	$genome_to_fams_hash->{$self->{_genome}} = join("\t", @all_fams);
    }
    
    return $genome_to_fams_hash;
}

1;
