# -*- perl -*-
#########################################################################
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
#########################################################################

package FIGM;

use FIGV;
use Carp qw(confess carp cluck);
use strict;
use FIG;
use FIG_Config;
use SFXlate;
use SproutFIG;
use Tracer;
use Data::Dumper;
use vars qw($AUTOLOAD);
use DB_File;
use FileHandle;

our $persistent_figv_cache = {};

#
# Create a new FIGM.
# Since creating a FIGV is only a data structure manipulation
# we go ahead and create one for each orgdir listed. We need
# to poke at the orgdir to find the genome id that it represents
# anyway.
#
sub new {
    my($class, $fig, @org_dirs) = @_;
    # cluck Dumper(FIGM => \@org_dirs);
    if (!ref($fig))
    {
	$fig = new FIG;
    }

    my $self         = {};
    $self->{_fig}    = $fig;
    $self->{_org_dirs} =  [@org_dirs];

    $self->{_figv_cache} = {};
    $self->{_peer_org_dir} = {};

    bless $self, $class;

    for my $dir (@org_dirs)
    {
	my $figv;
	if ($FIG_Config::use_figv_singleton)
	{
	    if (my $fv = $persistent_figv_cache->{$dir})
	    {
		$figv = $fv;
	    }
	}

	if (!$figv)
	{
	    $figv = eval { new FIGV($dir, undef, $fig); };
	}
	
	if ($@)
	{
	    warn "Error creating FIGV for $dir: $@";
	}
	elsif ($figv)
	{
	    if ($FIG_Config::use_figv_singleton)
	    {
		$persistent_figv_cache->{$dir} = $figv;
	    }
	    $self->{_figv_cache}->{$figv->genome_id()} = $figv;
	    $self->{_peer_org_dir}->{$figv->genome_id()} = $dir;
	}
    }

    return $self;
}


sub DESTROY {
    my($self) = @_;
    print STDERR "FIGM destroy $self $$\n";
    $self->SUPER::DESTROY if $self->can("SUPER::DESTROY");
}

sub is_complete
{
    return 1;
}

#
# Redirect any method invocations that we don't implement out to the
# underlying FIG object.
#
sub AUTOLOAD
{
    my($self, @args) = @_;

    if (ref($self) ne "FIGM") {
	confess "BAD FIGM object passed to AUTOLOAD";
    }

    no strict 'refs';

    my $meth = $AUTOLOAD;
    $meth =~ s/.*:://;
    my $fmeth = "FIG::$meth";

    my $fig = $self->{_fig};
#    my $args = Dumper(\@args);
    if (wantarray)
    {
	my @res = $fig->$meth(@args);
#	warn "FIGV invoke $meth($args) returns\n", Dumper(\@res);
	return @res;
    }
    else
    {
	my $res = $fig->$meth(@args);
#	warn "FIGV invoke $meth($args) returns\n", Dumper($res);
	return $res;
    }
}

sub FIG
{
    my($self) = @_;
    return $self;
}

sub find_figv
{
    my($self, $genome) = @_;

    my $figv = $self->{_figv_cache}->{$genome};
    if (ref($figv))
    {
	return $figv;
    }
    else
    {
	return $self->{_fig};
    }
}

sub find_figv_for_fid
{
    my($self, $fid) = @_;
    if ($fid =~ /^fig\|(\d+.\d+)\./)
    {
	return $self->find_figv($1);
    }
    else
    {
	return $self->{_fig};
    }
}

sub sort_fids_by_taxonomy
{
    my($self,@fids) = @_;

    return map     { $_->[0] }
           sort    { $a->[1] cmp $b->[1] }
           map     { [$_,$self->taxonomy_of($self->genome_of($_))] }
           @fids;
}

sub genomes
{
    my($self, $complete) = @_;

    my $fig = $self->{_fig};
    my @base = $fig->genomes($complete);

    return @base, keys %{$self->{_figv_cache}};
}

sub genome_list
{
    my($self) = @_;
    
    my $genome_list = [];

    foreach my $id (keys %{$self->{_figv_cache}}) {
	push(@$genome_list, [ $id, $self->genus_species($id), $self->genome_domain($id) ]);
    }

    push(@$genome_list, @{$self->{_fig}->genome_list});

    return $genome_list;
}

sub get_basic_statistics
{
    my($self, $genome) = @_;

    my $figv = $self->find_figv($genome);

    return $figv->get_basic_statistics($genome);
}


sub get_peg_statistics {
    my ($self, $genome) = @_;

    my $figv = $self->find_figv($genome);
    return $figv->get_peg_statistics($genome);
}

sub genome_id_to_genome_object
{
    my ($self, $genome) = @_;

    my $figv = $self->find_figv($genome);
    return &FIG::genome_id_to_genome_object($figv, $genome);
}

#
# To retrieve a subsystem in FIGV, we create the subsystem as normal via $fig->get_subsystem,
# then insert the row for the virtual org dir we are processing.
#
# The FIGM solution needs work.
#

sub get_subsystem
{
    my($self,$ssa) = @_;

    my $fig     = $self->{_fig};

    # get the subsystem data from the seed
    my $ss = $fig->get_subsystem($ssa);

    # get the roles
    my @roles = $ss->get_roles();
    my @non_aux_roles = grep { ! $ss->is_aux_role($_) } @roles;

    my $vcodes;
    my $bindings;
    my $variant;
    
    # go through all figvs in this figm
    foreach my $org_id (keys(%{$self->{_figv_cache}}))
    {
	my $figv = $self->{_figv_cache}->{$org_id};
	$bindings = {};

	if ($figv->need_bindings_recomputed($ssa))
	{
	    if (!$vcodes)
	    {
		$vcodes = &collect_vcodes($ss, \@non_aux_roles);
	    }
	    # recalculate the bindings for this subsystem
	    foreach my $role (@roles) {
		my @pegs = $figv->seqs_with_role($role, undef, $org_id);
		$bindings->{$role} = \@pegs;
	    }
	    
	    # calculate the variant
	    my $variant = '-1';
	    my @roles_in_this_genome = ();
	    foreach my $role (@non_aux_roles) {
		if (scalar(@{$bindings->{$role}})) {
		    push(@roles_in_this_genome,$role);
		}
	    }
	    
	    #
	    # Use the key for this genome to find the best variant code in the
	    # SEED subsystem from %vcode.
	    #

	    my $key = join("\t",sort @roles_in_this_genome);
	    my $n;
	    my $bestN = 0;
	    my $bestK = undef;
	    my $matches = $vcodes->{$key};
	    
	    unless ($matches) {
		foreach my $key1 (keys(%$vcodes)) {
		    if (&not_minus_1($vcodes->{$key1}) && (length($key) > length($key1)) && ($n = &contains($key,$key1)) && ($n > $bestN)) {
			$bestN = $n;
			$bestK = $key1;
		    }
		}
		if ($bestK) {
		    $matches = $vcodes->{$bestK};
		}
	    }
	    
	    if ($matches) {
		my @vcs  = sort { ($vcodes->{$b} <=> $vcodes->{$a}) or ($b cmp $a) } keys(%$matches);
		$variant = $vcs[0];
	    }
	    
	    $figv->update_bindings_and_variant($ssa, $bindings, $variant);
	    $figv->clear_need_bindings_recomputed($ssa);
	}
	else
	{
	    ($variant, $bindings) = $figv->get_variant_and_bindings($ssa);
	}
	$ss->add_virtual_genome($figv->genus_species($org_id), $org_id, $variant, $bindings);
    }

    return $ss;
}

    #
    # determine the variants implemented in the subsystem
    #
    # collect into %vcodes.
    #   key is a tab-joined list of roles implemented in this genome
    #   value is a hash of variant codes for the key.
    #
sub collect_vcodes
{
    my($ss, $non_aux_roles) = @_;
    
    my $vcodes = {};
    foreach my $genome ($ss->get_genomes) {
	my @roles_in_genome = ();
	my $vcode = $ss->get_variant_code($ss->get_genome_index($genome));
	next if (($vcode eq '0') || ($vcode =~ /\*/));
	
	foreach my $role (@$non_aux_roles) {
	    my @pegs = $ss->get_pegs_from_cell($genome,$role);
	    if (@pegs > 0) {
		push(@roles_in_genome,$role);
	    }
	}
	
	my $key = join("\t",sort @roles_in_genome);
	$vcodes->{$key}->{$vcode}++;
    }
    return $vcodes;
}

# returns undef if $k2 is not a subset of $k1. If it is, it returns the size of $k2
sub contains {
    my($k1,$k2) = @_;

    my %s1 = map { $_ => 1 } split(/\t/,$k1);
    my @s2 = split(/\t/,$k2);
    my $i;
    for ($i=0; ($i < @s2) && $s1{$s2[$i]}; $i++) {}
    return ($i < @s2) ? undef : scalar @s2;
}

sub not_minus_1 {
    my($hits) = @_;

    my @poss = keys(%$hits);
    my $i;
    for ($i=0; ($i < @poss) && ($poss[$i] eq "-1"); $i++) {}
    return ($i < @poss);
}

sub active_subsystems
{
    my($self, $genome, $all) = @_;

    my $figv = $self->find_figv($genome);
    return $figv->active_subsystems($genome, $all);
}

sub seqs_with_role {
    my($self,$role,$who,$genome) = @_;

    my $figv = $self->find_figv($genome);
    return $figv->seqs_with_role($role,$who,$genome);
}

sub subsystems_for_peg_complete {
    my ($self, $peg) = @_;

    my $figv = $self->find_figv_for_fid($peg);
    return $figv->subsystems_for_peg_complete($peg);
}

sub subsystems_for_pegs_complete {
    my ($self, $peg, $include_aux) = @_;

    # divide the pegs into rast organisms and have one entry for all seed orgs
    my $orgs = { 'other' => [] };
    foreach my $p (@$peg) {
	my ($org) = $p =~ /^fig\|(\d+\.\d+)/;
	if (exists($self->{_figv_cache}->{$org})) {
	    if (exists($orgs->{$org})) {
		push(@{$orgs->{$org}}, $p);
	    } else {
		$orgs->{$org} = [ $p ];
	    }
	} else {
	    push(@{$orgs->{'other'}}, $p);
	}
    }
    
    # initialize the return variable
    my %results;
    
    # get the ss for each organism and push it 
    foreach my $key (keys(%$orgs)) {
	if ($key eq 'other') {
	    my %ret = $self->{_fig}->subsystems_for_pegs_complete($orgs->{$key}, $include_aux);
	    foreach my $key (keys(%ret)) {
		$results{$key} = $ret{$key};
	    }
	} else {
	    foreach my $p (@{$orgs->{$key}}) {
		my @ret = $self->{_figv_cache}->{$key}->subsystems_for_peg_complete($p);
		foreach my $entry (@ret) {
		    push(@{$results{$entry->[3]}}, [ $entry->[0], $entry->[1], $entry->[2] ]);
		}
	    }
	}
    }

    return %results;
}

sub protein_subsystem_to_roles {
    my ($self,$peg,$subsystem) = @_;

    my $figv = $self->find_figv_for_fid($peg);
    return $figv->protein_subsystem_to_roles($peg,$subsystem);
}

sub genus_species {
    my($self,$genome) = @_;

    my $figv = $self->find_figv($genome);
    return $figv->genus_species($genome);
}

sub get_genome_assignment_data {
    my($self,$genome) = @_;

    my $figv = $self->find_figv($genome);
    return $figv->get_genome_assignment_data($genome);
}

sub org_of {
    my($self,$peg) = @_;

    if ($peg =~ /^fig\|(\d+\.\d+)\.peg\.\d+/)
    {
	return $self->genus_species($1);
    }
    return "";
}

sub get_genome_subsystem_data {
    my($self,$genome) = @_;
    my $figv = $self->find_figv($genome);
    return $figv->get_genome_subsystem_data($genome);
}

sub get_genome_subsystem_count
{
    my($self,$genome) = @_;

    my $figv = $self->find_figv($genome);
    return $figv->get_genome_subsystem_count($genome);
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

#sub protein_subsystem_to_roles {
#die;
#}

sub contig_lengths {
    my ($self, $genome) = @_;
    my $figv = $self->find_figv($genome);
    return $figv->contig_lengths($genome);
}

sub contig_ln {
    my ($self, $genome, $contig) = @_;

    my $figv = $self->find_figv($genome);
    return $figv->contig_ln($genome, $contig);
}

sub contigs_of
{
    my ($self, $genome) = @_;
    my $figv = $self->find_figv($genome);
    return $figv->contigs_of($genome);
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

    my $figv = $self->find_figv($genome);
    return $figv->dna_seq($genome, @locations);
}

sub genome_szdna {
    my ($self, $genome) = @_;

    my $figv = $self->find_figv($genome);
    return $figv->genome_szdna($genome);
}

sub genome_version {
    my ($self, $genome) = @_;

    my $figv = $self->find_figv($genome);
    return $figv->genome_version($genome);
}

sub genome_pegs {
    my ($self, $genome) = @_;
    my $figv = $self->find_figv($genome);
    return $figv->genome_pegs($genome);
}

sub genome_rnas {
    my ($self, $genome) = @_;

    my $figv = $self->find_figv($genome);
    return $figv->genome_rnas($genome);
}

sub genome_domain {
    my ($self, $genome) = @_;
    my $figv = $self->find_figv($genome);
    return $figv->genome_domain($genome);
}

sub genes_in_region {
    my($self,$genome,$contig,$beg,$end) = @_;
    my $figv = $self->find_figv($genome);
    return $figv->genes_in_region($genome,$contig,$beg,$end);
}

sub overlaps {
    my($b1,$e1,$b2,$e2) = @_;

    if ($b1 > $e1) { ($b1,$e1) = ($e1,$b1) }
    if ($b2 > $e2) { ($b2,$e2) = ($e2,$b2) }
    return &FIG::between($b1,$b2,$e1) || &FIG::between($b2,$b1,$e2);
}

sub all_contigs {
    my($self,$genome) = @_;
    my $figv = $self->find_figv($genome);
    return $figv->all_contigs($genome);
}

sub all_features {
    my($self,$genome,$type) = @_;
    my $figv = $self->find_figv($genome);
    return $figv->all_features($genome,$type);
}

sub all_features_detailed_fast {
    my($self,$genome, $regmin, $regmax, $contig) = @_;
    my $figv = $self->find_figv($genome);
    return $figv->all_features_detailed_fast($genome, $regmin, $regmax, $contig);
}

sub compute_clusters {
    # Get the parameters.
    my ($self, $pegList, $subsystem, $distance) = @_;
    if (! defined $distance) {
        $distance = 5000;
    }

    my($peg,%by_contig);

    my @locs = $self->feature_location_bulk($pegList);

    for my $ent (@locs)
    {
	my($peg, $loc) = @$ent;
        if ($loc)
        {
            my ($contig,$beg,$end) = &FIG::boundaries_of($loc);
            my $genome = &FIG::genome_of($peg);
            push(@{$by_contig{"$genome\t$contig"}},[($beg+$end)/2,$peg]);
	}
    }
	
#     foreach $peg (@$pegList)
#     {
#         my $loc;
#         if ($loc = $self->feature_location($peg))
#         {
#             my ($contig,$beg,$end) = &FIG::boundaries_of($loc);
#             my $genome = &FIG::genome_of($peg);
#             push(@{$by_contig{"$genome\t$contig"}},[($beg+$end)/2,$peg]);
#         }
#     }

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
    return $fig->boundaries_of(@args);
}

sub feature_location {
    my($self,$fid) = @_;

    my $figv = $self->find_figv_for_fid($fid);
    return scalar $figv->feature_location($fid);
}

sub feature_location_bulk {
    my($self,$fids) = @_;

    my $fig     = $self->{_fig};

    my @ids;
    my @out;
    for my $fid (@$fids)
    {
	my $figv = $self->find_figv_for_fid($fid);
	if ($figv)
	{
	    push(@out, [$fid, scalar $figv->feature_location($fid)]);
	}
	else
	{
	    push(@ids, $fid);
	}
    }
    push(@out, $fig->feature_location_bulk(\@ids));
    return @out;
}

sub function_of {
    my($self,$fid) = @_;

    my $fig     = $self->{_fig};

    my $figv = $self->find_figv_for_fid($fid);
    return $figv->function_of($fid);
}

sub function_of_bulk
{
    my($self, $fid_list) = @_;

    my @for_fig;

    my $out = {};

    my $fallback_fig = $self->{_fig};
    
    for my $fid (@$fid_list)
    {
        my $fid_fig = $self->find_figv_for_fid($fid);
	if ($fid_fig == $fallback_fig)
	{
	    push(@for_fig, $fid);
	}
	else
	{
	    my $fn = $fid_fig->function_of($fid);
	    $out->{$fid} = $fn if defined($fn);
	}
    }

    my $others = $fallback_fig->function_of_bulk(\@for_fig);
    $out->{$_} = $others->{$_} for keys %$others;
    return $out;
}

sub assign_function
{
    my($self,$fid, $user, $function, $confidence) = @_;

    my $fig     = $self->{_fig};

    my $figv = $self->find_figv_for_fid($fid);
    return $figv->assign_function($fid, $user, $function, $confidence);
}

sub add_annotation
{
    my($self, $feature_id,$user,$annotation, $time_made) = @_;

    my $fig     = $self->{_fig};

    my $figv = $self->find_figv_for_fid($feature_id);
    return $figv->add_annotation($feature_id,$user,$annotation, $time_made);
}

sub feature_aliases {
    my($self,$fid) = @_;
    my $figv = $self->find_figv_for_fid($fid);
    return $figv->feature_aliases($fid);
}

sub feature_annotations {
    my($self,$fid,$rawtime) = @_;
    my $figv = $self->find_figv_for_fid($fid);
    return $figv->feature_annotations($fid);
}

sub get_translation {
    my($self,$peg) = @_;
    my $figv = $self->find_figv_for_fid($peg);
    return $figv->get_translation($peg);
}

sub translation_length
{
    my($self, $peg) = @_;
    my $figv = $self->find_figv_for_fid($peg);
    return $figv->translation_length($peg);
}

sub translatable
{
    my($self, $peg) = @_;
    my $figv = $self->find_figv_for_fid($peg);
    return $figv->translatable($peg);
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

    my $figv = $self->find_figv_for_fid($fid);
    return $figv->is_real_feature($fid);

}

sub is_deleted_fid
{
    my($self, $fid) = @_;

    my $figv = $self->find_figv_for_fid($fid);
    return $figv->is_deleted_fid($fid);

}

sub pegs_of
{
    my($self, $genome) = @_;
    my $figv = $self->find_figv($genome);
    return $figv->pegs_of($genome);
}


sub rnas_of
{
    my($self, $genome) = @_;
    my $figv = $self->find_figv($genome);
    return $figv->pegs_of($genome);
}

sub is_virtual_feature
{
    my($self, $peg) = @_;
    my $figv = $self->find_figv_for_fid($peg);

    return ref($figv) =~ /FIGV/ ? 1 : 0;
}

sub bbhs
{
    my($self,$peg,$cutoff) = @_;

    my $figv = $self->find_figv_for_fid($peg);

    return $figv->bbhs($peg, $cutoff);
}

=head3 sims

Retrieve sims. We partition the sims into sets based on their genome id, associating
them with the figv for the genome that they are part of. For the figv pegs
we use the sims_for_figm method that will retrieve the sims for the peer 
organims we've configured, as well as the usual figv sims.

=cut
    
sub sims
{
    my($self,$pegarg,$max,$maxP,$select, $max_expand, $filters) = @_;

    my %pegs_per_genome;
    my @other;

    for my $peg (ref($pegarg) ? @$pegarg : ($pegarg))
    {
	if ($peg =~ /^fig\|(\d+\.\d+)/)
	{
	    push(@{$pegs_per_genome{$1}}, $peg);
	}
	else
	{
	    push(@other, $peg);
	}
    }

    my @out;

    # print STDERR "SIMS: " . Dumper(\%pegs_per_genome);

    while (my($genome, $pegs) = each(%pegs_per_genome))
    {
	my $figv = $self->find_figv($genome);

	my @sims;

	if (ref($figv) eq "FIGV")
	{
	    @sims= $figv->sims_for_figm([keys %{$self->{_figv_cache}}], $pegs, $max, $maxP, $select, $max_expand, $filters);
	}
	else
	{
	    @sims= $figv->sims($pegs, $max, $maxP, $select, $max_expand, $filters);
	}

	push(@out, @sims);
    }
    push(@out, $self->{_fig}->sims(\@other, $max, $maxP, $select, $max_expand, $filters));
    # print STDERR "FINAL: " . Dumper(\@out);
    return @out;
}

sub coupled_to
{
    my($self,$peg) = @_;

    my $figv = $self->find_figv_for_fid($peg);
    return $figv->coupled_to($peg);
}

sub coupled_to_batch
{
    my($self, @pegs) = @_;

    return () unless @pegs;

    # divide the pegs into rast organisms and have one entry for all seed orgs
    my $orgs = { 'other' => [] };
    foreach my $peg (@pegs) {
	my ($org) = $peg =~ /^fig\|(\d+\.\d+)/;
	if (exists($self->{_figv_cache}->{$org})) {
	    if (exists($orgs->{$org})) {
		push(@{$orgs->{$org}}, $peg);
	    } else {
		$orgs->{$org} = [ $peg ];
	    }
	} else {
	    push(@{$orgs->{'other'}}, $peg);
	}
    }

    # initialize the return variable
    my $ret = [];
    
    # get the fc for each organism and push it 
    foreach my $key (keys(%$orgs)) {
	if ($key eq 'other') {
	    push(@$ret, $self->{_fig}->coupled_to_batch(@{$orgs->{$key}}));
	} else {
	    push(@$ret, $self->{_figv_cache}->{$key}->coupled_to_batch(@{$orgs->{$key}}));
	}
    }
    
    return @$ret;
}

sub coupling_evidence
{
    my($self,$peg1, $peg2) = @_;
    my $figv = $self->find_figv_for_fid($peg1);
    return $figv->coupling_evidence($peg1, $peg2);
}

sub coupling_and_evidence
{
    my($self,$peg1) = @_;
    my $figv = $self->find_figv_for_fid($peg1);
    return $figv->coupling_and_evidence($peg1);

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

    my $figv = $self->find_figv_for_fid($peg1);
    return $figv->in_pch_pin_with_and_evidence($peg1);
}

sub get_attributes {
    my($self, $id, $attr) = @_;

    my $fig     = $self->{_fig};

    if ($id =~ /^fig\|(\d+\.\d+)\.peg\.\d+$/ or $id =~ /^(\d+\.\d+)/)
    {
	my $figv = $self->find_figv($1);
	return $figv->get_attributes($id, $attr);
    }
    else
    {
	return $fig->get_attributes($id, $attr);
    }
}

sub taxonomy_of {
    my($self,$genome) = @_;
    my $figv = $self->find_figv($genome);
    return $figv->taxonomy_of($genome);
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

  my $figv = $self->find_figv($organism);
  my $directory = $figv->scenario_directory($organism);
  return $directory;
}

sub find_role_in_org {
    my ($self, $role, $organism, $user, $cutoff) = @_;

    my $figv = $self->find_figv($organism);
    return $figv->find_role_in_org($role, $organism, $user, $cutoff);
}

sub organism_directory {
    my ($self, $organism) = @_;

    my $figv = $self->find_figv($organism);
    return $figv->organism_directory($organism);
}

sub delete_feature {
    my($self,$user,$fid) = @_;

    my ($organism) = $fid =~ /fig\|(\d+\.\d+)/;
    my $figv = $self->find_figv($organism);
    return $figv->delete_feature($user, $fid);
}

sub add_feature {
    my( $self, $user, $organism, $type, $location, $aliases, $sequence) = @_;

    my $figv = $self->find_figv($organism);
    return $figv->add_feature($user, $organism, $type, $location, $aliases, $sequence);
}

sub genome_info {
    my ($self) = @_;

    my $info = $self->{_fig}->genome_info;
    foreach my $id (keys(%{$self->{_figv_cache}})) {
	my $f = $self->{_figv_cache}->{$id};
	push(@$info, [ $id, "Private: ".$f->genus_species($id), $f->genome_szdna($id), $f->genome_domain($id), $f->genome_pegs($id), $f->genome_rnas($id), $f->is_complete($id), $f->taxonomy_of($id) ]);
    }

    return $info;
}

sub check_db_peg_to_fams {
    my ($self, $peg_to_fams_hash) = @_;

    foreach my $id (keys(%{$self->{_figv_cache}})) {
	my $f = $self->{_figv_cache}->{$id};
	$f->check_db_peg_to_fams($peg_to_fams_hash);
    }

    return $peg_to_fams_hash;
}

sub check_db_genome_to_fams {
    my ($self, $genome_to_fams_hash) = @_;

    foreach my $id (keys(%{$self->{_figv_cache}})) {
	my $f = $self->{_figv_cache}->{$id};
	$f->check_db_peg_to_fams($genome_to_fams_hash);
    }

    return $genome_to_fams_hash;
}

1;
