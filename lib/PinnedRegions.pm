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

package PinnedRegions;

use strict;
use warnings;

use Data::Dumper;

use FIG;
use FIG_Config;
use BasicLocation;
use Tracer;

our $have_kmer_client;
eval {
    require KmerClient;
    $have_kmer_client = 1;
};

our $have_cds;
eval {
    require ConservedDomainSearch;
    $have_cds = 1;
};

use Time::HiRes qw( usleep gettimeofday tv_interval );

sub pinned_regions {
    my($fig, $pin_desc, $fast_color, $sims_from, $map_sz, $add_features, $extended, $color_by_function, $peg_functions, $fids_for_cds) = @_;

    my $cds;
    if ($have_cds)
    {
	$cds = ConservedDomainSearch->new($fig);
    }

    my $align = $pin_desc->{pin_alignment} || 'center';
    
    Trace("Pinned regions method called.") if T(3);
    # Get list of pegs required by the description in $pin_desc
    my $pinned_pegs = &expand_peg_list($fig, $pin_desc);
    Trace("Pinned pegs are " . join(", ", @$pinned_pegs) . ".") if T(3);
    # Filter out the pegs that don't exist.
    # Get regions around pinned pegs -- boundaries, features, etc.
    my $regions = &define_regions($fig, $map_sz, $pinned_pegs, $align);
    # Filter out overlapping regions caused by gene fusions, frame shifts etc. where multiple PEGs
    # are pinned (by similarity or PCH) to the input PEG
    $regions = &filter_regions_1($pin_desc, $regions);

    #
    # Add contig info to regions
    #
    my @cinfo = map { [ @$_{qw(genome_id contig)} ] } @$regions;
    my $contig_lens = $fig->contig_ln_bulk(\@cinfo);
    $_->{contig_length} = $contig_lens->{$_->{genome_id}}->{$_->{contig}} foreach @$regions;
    # print STDERR Dumper($regions);

    # Get information for each feature -- location, strand, function etc.
    my $feature_data = &feature_data($fig, $regions, $extended, $peg_functions);

    &add_functional_coupling($fig, $pin_desc, $regions, $feature_data);

#    &add_figfams($fig, $feature_data);
    &add_pattyfams($fig, $feature_data);
    &add_subsystem_data($fig, $pin_desc, $feature_data, $extended);

    my $cdd_trans = {};
    if ($cds && ref($fids_for_cds) eq 'ARRAY' && @$fids_for_cds)
    {
	$cdd_trans = &add_cdd($regions, $fig, $cds, $fids_for_cds, $feature_data);
    }

    Trace("Coloring pegs.") if T(3);
    # Assign a set number to some PEGs through transitive closure based on similarity, from blast scores
    if ($color_by_function)
    {
	&color_pegs_by_function($fig, $pin_desc, $pinned_pegs, $regions, $feature_data, $fast_color, $sims_from);
    }
    else
    {
	&color_pegs($fig, $pin_desc, $pinned_pegs, $regions, $feature_data, $fast_color, $sims_from, $cdd_trans);
    }

    # Filter out regions which have only a single PEG (the pinned one) colored.
    # $regions = &filter_regions_2($pin_desc, $regions, $feature_data);

    if ( defined( $add_features ) ) {
        &add_features_to_fdata( $fig, $feature_data, $add_features, $regions );
    }

    # Add feature data to the regions to make the final maps
    my $maps = &make_maps($fig, $regions, $feature_data);

    return $maps;
}

sub expand_peg_list {
    my($fig, $pin_desc) = @_;
    my @pegs = ();

    if ($pin_desc->{pegs}->[0] =~ /^fig/ && @{$pin_desc->{pegs}} > 1)
    {
	@pegs = @{$pin_desc->{pegs}};
	return \@pegs;
    }
    # Filter for legitimate genome IDS -- should handle deleted organisms and (in NMPDR) non-NMPDR orgs
    my %ok_genome_id = map {$_ => 1} $fig->genomes;

    # If the user has selected the genomes to be displayed, filter out all other genomes
    if ( @{ $pin_desc->{'show_genomes'} } )
    {
        # create new %ok_genome_id hash from intersection of selected genomes and legitimate genome IDs
        %ok_genome_id = map {$_ => 1} grep {$ok_genome_id{$_}} @{ $pin_desc->{'show_genomes'} };
    }

    if ( @{ $pin_desc->{'pegs'} } > 1 )
    {
        # Input already contains list of pegs, no need for expansion
	if ($pin_desc->{'pegs'}->[0] =~ /^\d/) {
	    @pegs = @{ $pin_desc->{'pegs'} };
	} else {
	    @pegs = grep {$ok_genome_id{$fig->genome_of($_)}} @{ $pin_desc->{'pegs'} };
	}
    }
    else
    {
        # Get PCH pinned pegs
        my $pch_pinned_pegs = &pch_pinned_pegs($fig, $pin_desc, \%ok_genome_id);

        $pch_pinned_pegs = &select_pegs($fig, $pin_desc, $pin_desc->{'n_pch_pins'}, $pch_pinned_pegs);

	# Get similarity pinned pegs
	my $sim_pegs = &sim_pegs($fig, $pin_desc, \%ok_genome_id);

        # Remove PCH PEGs (if any) from similarity pinned PEG list
        my %seen  = map {$_ => 1} @$pch_pinned_pegs;
        $sim_pegs = [grep {! $seen{$_}} @$sim_pegs];

	# print STDERR Dumper(simpegs1 => $sim_pegs);
	my $n_pegs = $pin_desc->{'n_sims'} + $pin_desc->{n_kmers} + $pin_desc->{n_plfams} + $pin_desc->{n_pgfams};
        $sim_pegs = &select_pegs($fig, $pin_desc, $n_pegs, $sim_pegs);
	# print STDERR Dumper(simpegs2 => $sim_pegs);

        @pegs = ();

        # Get input peg
        my $peg = $pin_desc->{'pegs'}[0];
        if ( $pin_desc->{'sort_by'} eq 'phylogeny' )
        {
            # First add input peg, then sort by phylogeny
            @pegs = $fig->sort_fids_by_taxonomy($peg, @$sim_pegs, @$pch_pinned_pegs);
        }
        elsif ( $pin_desc->{'sort_by'} eq 'phylogenetic_distance' )
        {
            # Sort by phylogenetic distance or organism name
            my $g1 = $fig->genome_of($peg);
            @pegs  = map {$_->[0] } 
                       sort {$a->[1] <=> $b->[1] or $a->[2] cmp $b->[2]} 
                         map {[$_, $fig->crude_estimate_of_distance($g1,$fig->genome_of($_)), $fig->org_of($_)]} @$sim_pegs, @$pch_pinned_pegs;
            # Add input peg at front of list
            unshift @pegs, $peg;
        }
        elsif ( $pin_desc->{'sort_by'} eq 'similarity' )
	{
	    # Sort by similarity to input peg (first peg in list of arguments)
            @pegs = &sort_pegs_by_similarity($fig, $peg, @$sim_pegs, @$pch_pinned_pegs);
	}
	
	if ( $n_pegs < @pegs )
	{
	    $#pegs = $n_pegs - 1;
	}
    }
    return \@pegs
}

sub pch_pinned_pegs {
    my($fig, $pin_desc, $ok_genome_id) = @_;

    my @pegs = ();

    if ( $pin_desc->{'n_pch_pins'} > 0 )
    {
        # Get input peg
        my $peg = $pin_desc->{'pegs'}[0];

        # Get pegs in pch pin with input peg, and filter out deleted or non-NMPDR organisms
        @pegs = grep {$ok_genome_id->{$fig->genome_of($_)}} 
                  grep {$_ ne $peg} 
                    map {$_->[0]} $fig->in_pch_pin_with_and_evidence($peg);
    }

    return \@pegs;
}

sub sim_pegs {
    my($fig, $pin_desc, $ok_genome_id) = @_;

    # Return a list of PEGS similar to the input PEG, with evalue less than or equal to 
    # the user defined similarity cutoff.
    # If the number of sims requested is 0, return the empty list 

    my $n_sims     = $pin_desc->{'n_sims'};
    my $n_kmers     = $pin_desc->{'n_kmers'};
    my $n_plfams     = $pin_desc->{'n_plfams'};
    my $n_pgfams     = $pin_desc->{'n_pgfams'};
    my $sim_cutoff = $pin_desc->{'sim_cutoff'};
    my @pegs = ();
    
    my $peg = $pin_desc->{'pegs'}[0];
    if ( $n_sims > 0 )
    {
        my %seen;

        @pegs = grep {$ok_genome_id->{$fig->genome_of($_)}} 
                  grep {! $seen{$_}++}
                    map {$_->[1]} $fig->sims($peg, 10000, $sim_cutoff, "fig");
    }
    elsif ($n_kmers > 0 && $have_kmer_client)
    {
        my %seen;

	my $kc;
	eval {
	    $kc = KmerClient->new($fig,
				  $FIG_Config::kmer_server_host,
				  $FIG_Config::kmer_server_port,
				  $FIG_Config::kmer_genome_map,
				  $FIG_Config::kmer_memcache_host,
				  $FIG_Config::kmer_memcache_port);
	};
	if ($@)
	{
	    warn "Cannot create kmer_client: $@";
	}
	else
	{
	    my @kp = $kc->get_close_pegs($peg, 10000);
	    
	    @pegs = grep {$ok_genome_id->{$fig->genome_of($_)}} grep {! $seen{$_}++} map {$_->[1]} @kp;

	    # open(O, ">","/homes/olson/tmp/kmer.out");
	    # print O Dumper($peg, \@kp, \@pegs);
	    # close(O);
	}
    }
    elsif ($n_pgfams > 0 || $n_plfams > 0)
    {
	my $re = ($n_pgfams > 0) ? "LIKE 'PGF%'" : "LIKE 'PLF%'";
	my $res = $fig->db_handle->SQL("SELECT family, score FROM family_membership WHERE fid = ? AND family $re", 0, $peg);

	# print STDERR Dumper('plfams', $n_pgfams, $n_plfams, $peg, $re, $res);
	if (@$res)
	{
	    my $fam = $res->[0]->[0];
	    my $score = $res->[0]->[1];
	    my $members = $fig->db_handle->SQL('SELECT fid, abs(score - ?) FROM family_membership WHERE family = ? ORDER BY score', undef,
					       $score, $fam);
#	    my $members = $fig->db_handle->SQL('SELECT fid, abs(score - ?) FROM family_membership WHERE family = ? ORDER BY abs(score - ?) DESC', undef,
#					       $score, $fam, $score);
	    @pegs = map { $_->[0] } @$members;
	    # print STDERR Dumper('plfams2', $fam, $score, $members, \@pegs);
	}
    }
    
    return \@pegs;
}

sub select_pegs {
    my($fig, $pin_desc, $n_pegs, $pegs) = @_;

    if ( $n_pegs == 0 )
    {
        return [];
    }

    if ( $pin_desc->{'collapse_close_genomes'} == 1 )
    {
        # To collapse the PEG list, we want to:
        # a. Return at most $n_pegs PEGs, 
        # b. include a representative PEG from every genus, subject to a, and 
        # c. include a representative PEG from each genus-species, subject to a and b.
        # Individual strains will get represented by a single PEG, chosen arbitrarily.

        # The PEG selection is defined by the order of the PEGs. This could be done beter.
        
        my(%seen, @unique_genus, @unique_genus_species);
        
        foreach my $peg ( @$pegs )
        {
            my $org = $fig->org_of($peg);
            # 'org_of' returns undef for deleted PEGs, these need to be omitted from the peg list returned
            
            if ( $org )
            {
                # Use only genus+species to drop strain information
                my($genus, $species) = split(/\s+/, $org);
                my $gs = "$genus $species";

                if ( not $seen{$genus}++ ) 
                {
                    # First PEG from this genus, add it to @unique_genus.
                    # Mark the genus+species as seen.
                    # A subsequent PEG with the same genus and species will be dropped.
                    # A subsequent PEG with the same genus but different species will be added to @unique_genus_species.
                    $seen{$gs}++;
                    push @unique_genus, $peg;
                } 
                elsif ( not $seen{$gs}++ ) 
                {
                    # First PEG from this genus+species, add it to @unique_genus_species.
                    push @unique_genus_species, $peg;
                }
            }
	}

        # Keep the unique_genus PEGS at the top, followed by the unique_genus_species
        $pegs = [@unique_genus, @unique_genus_species];
    }
    elsif ( $pin_desc->{'collapse_close_genomes'} == 2 )
    {
	#
	# Reduce the list to contain only one instance of any given tax id.
	#

	my %by_tax;
        
	my %keep;
	
        foreach my $peg ( @$pegs )
        {
            my ($org, $tax) = $peg =~ /^fig\|((\d+)\.\d+)/;
	    if ($tax =~ /^666666(6?)/)
	    {
		$keep{$peg}++;
		next;
	    }
            
            if ( $tax )
            {
		if (defined(my $last = $by_tax{$tax}))
		{
		    my($last_org, $last_tax) = @$last;
		    if (&FIG::by_genome_id($org, $last_org) > 0)
		    {
			$by_tax{$tax} = [$org, $tax, $peg];
		    }
		}
		else
		{
		    $by_tax{$tax} = [$org, $tax, $peg];
		}
            }
        }
	$keep{$_} = 1 foreach map { $_->[2] } values %by_tax;
	@$pegs = grep { $keep{$_} } @$pegs;
    }
    elsif ( $pin_desc->{'collapse_close_genomes'} == 3 )
    {
	#
	# Use the OTU reps to minimize the list.
	#
	
	my $gfile = "$FIG_Config::global/genome.sets";
	if (open(G, "<", $gfile))
	{
	    my %rep;
	    my %fam;
	    while (<G>)
	    {
		chomp;
		my($fam, $gid, $gname) = split(/\t/);
		if (!$rep{$fam})
		{
		    $rep{$fam} = $gid;
		    $fam{$gid} = $fam;
		}
	    }
	    close(G);

	    my @keep;
	    foreach my $peg ( @$pegs )
	    {
		my ($org, $tax) = $peg =~ /^fig\|((\d+)\.\d+)/;
		push(@keep, $peg) if defined($fam{$org});
	    }
	    @$pegs = @keep;
	}
    }

    # Truncate list if necessary
#    if ( $n_pegs < @$pegs )
#    {
#        $#{ $pegs } = $n_pegs - 1;
#    }
    
    return $pegs;
}

sub sort_pegs_by_similarity {
    my($fig, $input_peg, @other_pegs) = @_;

    # Set blast environment variables
    $ENV{"BLASTMAT"} ||= "$FIG_Config::blastmat";
    if ($ENV{"PATH"} !~ /fig\/bin/) { $ENV{"PATH"} = "$FIG_Config::bin:" . $ENV{"PATH"}; }

    my $temp = $FIG_Config::fast_temp // $FIG_Config::temp;

    # Get amino acid sequences
    my $sequences = &get_peg_sequences($fig, [$input_peg, @other_pegs]);

    #print STDERR Dumper(SORT => $input_peg, $sequences->{$input_peg}, $fig);

    # Write the input peg sequence to a fasta file
    my $input_seq = {$input_peg => $sequences->{$input_peg}};
    my $input_seq_file = "$temp/input_seq.$$.tmp.fasta";
    &write_seqs_to_file($input_seq, $input_seq_file);

    my @sorted = ($input_peg);

    # Write the other peg sequences to a fasta file
    delete $sequences->{$input_peg};
    my $other_seq_file = "$temp/other_seq.$$.tmp.fasta";
    if (%$sequences)
    {
	&write_seqs_to_file($sequences, $other_seq_file);
	
	# Run formatdb
	&formatdb_file($fig, $other_seq_file);
	
	# BLAST input peg sequence against other peg sequences
	# set high evalue cutoff so all hits get reported
	my $cutoff = 10;
	my $hits = &blast_files($fig, $input_seq_file, $other_seq_file, $cutoff);
	
	#
	# We might have multiple hits; just use the first.
	#
	my %seen;
	foreach my $hit ( @$hits )
	{
	    # BLAST output is already sorted by similarity
	    my($peg2) = (split(/\s+/, $hit))[1];
	    
	    push @sorted, $peg2 unless $seen{$peg2}++;
	}
	
	for my $file ($input_seq_file, $other_seq_file)
	{
	    for my $suffix ('', qw(.phr .psq .pin))
	    {
		my $path = "$file$suffix";
		unlink($path);
	    }
	}
    }

    return @sorted;
}

sub filter_regions_1 {
    my($pin_desc, $regions) = @_;
    my $new_regions;

    # Filter out some overlapping regions caused by gene fusions, frame shifts etc. where multiple PEGs
    # are pinned (by similarity or PCH) to the input PEG
    # Note that all overlaps are NOT eliminated

    if ( @{ $pin_desc->{'pegs'} } == 1 )
    {
        # Input pin description is for a single input peg
        my %seen;
        
        foreach my $region ( @$regions )
        {
            my $pinned_peg = $region->{'pinned_peg'};
            
            # If the pinned peg from a region has already been 'seen', don't
            # add the region into the @$new_regions array
            
            if ( not $seen{$pinned_peg} )
            {
                push @$new_regions, $region;
                
                my $fids = $region->{'features'};
                foreach my $fid ( @$fids )
                {
                    $seen{$fid} = 1;
                }
            }
        }
    }
    else
    {
        # input PEGs is a list -- keep all of them
        $new_regions = $regions;
    }

    return $new_regions;
}

sub filter_regions_2 {
    my($pin_desc, $regions, $feature_data) = @_;
    my $new_regions;

    if ( @{ $pin_desc->{'pegs'} } == 1 )
    {
        # Input is single peg
        my $input_peg = $pin_desc->{'pegs'}[0];

        foreach my $region ( @$regions )
        {
            my $fids = $region->{'features'};
            my $n_in_sets = grep {$feature_data->{$_}{'set_number'}} @$fids;
            
            if ( $n_in_sets > 1 or $region->{'pinned_peg'} eq $input_peg )
            {
                # Regions should be displayed only if:
                # a. more than one PEG is in a 'colored set' OR
                # b. the pinned PEG is the input PEG
                push @$new_regions, $region;
            }
        }
    }
    else
    {
        # Input is a list of pegs -- keep all of the regions
        $new_regions = $regions
    }

    return $new_regions;
}

sub make_maps {
    my($fig, $regions, $feature_data) = @_;

    for my $i (0..$#$regions)
    {
	my $region = $regions->[$i];
	
        my $features = [];
        my $fids = $region->{'features'};

        foreach my $fid ( @$fids )
        {
            push @$features, $feature_data->{$fid};
#           if ( !( $fid =~ /fig/ ) ) {
#             print STDERR $fid." FIDMAKEMAPS\n";
#           }
        }

        $region->{'features'} = $features;
    }

    return $regions;
}

sub define_regions {
    my($fig, $map_sz, $pinned_pegs, $align) = @_;

    # Define the 'region' around the pinned peg, get features, and store some region information
    my $regions = [];
    my $half_map_sz = int($map_sz/2);

    # print STDERR Dumper($pinned_pegs, $align);

    my $align_offset;
    foreach my $peg ( @$pinned_pegs )
    {
	my $genome;
	my $loc;
	if ($peg =~ /^(fig\|)*(\d+\.\d+)[\.\:]{1}(.+\_\d+\_\d+)$/) {
	    $loc = $3;
	    $genome = $2;
	} else {
	    $loc = $fig->feature_location($peg);
	    $genome = $fig->genome_of($peg);
	}
        # We only proceed if a location was found. Failure to find a location indicates
        # that the feature is deleted.
        if ($loc) {
            my($contig, $beg, $end) = $fig->boundaries_of($loc);

	    #
	    # Ensure that the focus peg is centered, even if we are aligning on start or stop.
	    #
	    my $align_point;
	    if ($align eq 'start') {
		$align_point = $beg;
	    } elsif ($align eq 'center') {
		$align_point = int(($beg + $end)/2);
	    } else {
		$align_point = $end;
	    };

	    my $center = int(($beg + $end) / 2);

	    my $region_mid;
	    if (!defined($align_offset))
	    {
		if ($beg < $end) {
		    $align_offset = $center - $align_point;
		} else {
		    $align_offset = $align_point - $center;
		}
		$region_mid = $center;
	    }
	    else
	    {
		if ($beg < $end) {
		    $region_mid = $align_point + $align_offset;
		} else {
		    $region_mid = $align_point - $align_offset;
		}
		
	    }
	    print STDERR "$peg AO=$align_offset  AP=$align_point C=$center M=$region_mid\n";
	    
            my $region_beg = $region_mid - $half_map_sz;
            my $region_end = $region_mid + $half_map_sz;
    
            my($fids) = $fig->genes_in_region($genome, $contig, $region_beg, $region_end);
            $fids = [ grep { $fig->is_real_feature($_) } @$fids];
	    
            my $region = {};
            $region->{'genome_id'} = $fig->genome_of($peg);
            $region->{'org_name'}  = $fig->genus_species($region->{'genome_id'});
            $region->{'contig'}    = $contig;
            $region->{'beg'}       = $region_beg;
            $region->{'mid'}       = $region_mid;
            $region->{'end'}       = $region_end;
            $region->{'features'}  = $fids;
            $region->{'pinned_peg_strand'} = ($beg <= $end)? '+' : '-';
            $region->{'pinned_peg'} = $peg;
    
            push @$regions, $region;
        }
    }

    return $regions;
}

sub add_features_to_fdata {
    my ( $fig, $feature_data, $add_features, $regions ) = @_;
    foreach my $region ( @$regions ) {
        my $new_feats = $add_features->{ $region->{ 'genome_id' } };
        foreach my $nf ( @$new_feats ) {
            if ( $nf->{ 'contig' } eq $region->{ 'contig' } &&
                 $nf->{ 'start' } < $region->{ 'end' } && 
                 $nf->{ 'start' } > $region->{ 'beg' } ) { 
                $feature_data->{ $nf->{ 'name' } } = &new_feature_entry( $fig, $nf->{ 'name' }, $region->{ 'mid' }, $region->{ 'contig' }, $nf->{ 'start' }, $nf->{ 'stop' }, $nf->{ 'type' }, $nf->{ 'function' } );
                push @{ $region->{ 'features' } }, $nf->{ 'name' };
            }
        }
    }
}

sub new_feature_entry {
    my ( $fig, $fid, $region_mid, $contig, $beg, $end, $type, $func ) = @_;

    if ( !defined( $type ) ) {
        $type = 'unknown';
    }

    if ( !defined( $func ) ) {
        $func = '';
    }

    my($left, $right)       = sort {$a <=> $b} ($beg, $end);
    my $size                = $right - $left + 1;
    my $strand              = ($beg <= $end)? '+' : '-';
    my $offset              = int(($left + $right)/2) - $region_mid;
    my $offset_beg          = $left  - $region_mid;
    my $offset_end          = $right - $region_mid;

    return {
             'fid'        => $fid,
             'type'       => $type,
             'contig'     => $contig,
             'beg'        => $beg,
             'end'        => $end,
             'size'       => $size,
             'strand'     => $strand,
             'offset'     => $offset,
             'offset_beg' => $offset_beg,
             'offset_end' => $offset_end,
             'function'   => $func,
	     'location'   => $contig."_".$beg."_".$end
    };
}

sub feature_data {
    my($fig, $regions, $extended, $peg_functions) = @_;
    my %feature_data;

    my %all_fids;
    my %all_genomes;

    foreach my $region ( @$regions )
    {
        my $fids       = $region->{'features'};
	$all_fids{$_}++, $all_genomes{&FIG::genome_of($_)}++ foreach @$fids;
    }

    my @all_fids = keys(%all_fids);
    my @locations = $fig->feature_location_bulk(\@all_fids);

    my %locations = map { $_->[0] => $_->[1] } @locations;

    my $aliases = $fig->feature_aliases_bulk(\@all_fids);

    my $functions;
    if (ref($peg_functions) eq 'HASH')
    {
	$functions = $peg_functions;
    }
    elsif (ref($peg_functions) eq 'PG')
    {
	$functions = $peg_functions->function_of_bulk(\@all_fids);
    }
    else
    {
	$functions = $fig->function_of_bulk(\@all_fids);
    }

#    print STDERR Dumper(\%locations, $aliases, $functions);
    foreach my $region ( @$regions )
    {
        my $region_mid = $region->{'mid'};
        my $fids       = $region->{'features'};

        foreach my $fid ( @$fids )
        {
            # Get feature data if this feature has not been seen before, this avoids repeating
            #  this step when pegs occur in multiple regions
            if ( not exists $feature_data{$fid} )
            {
                $feature_data{$fid} = &feature_entry($fig, $fid, $region_mid, $extended, $functions, \%locations, $aliases);
            }
        }
    }

    if ($extended) {
	# get the evidence codes
	my @all_fids = keys(%feature_data);
	my @codes = grep { $_->[1] =~ /^evidence_code/i } $fig->get_attributes(\@all_fids, 'evidence_code');
	my $pretty_codes = {};
	foreach my $code (@codes) {
	    my $pretty_code = $code->[2];
	    if ($pretty_code =~ /;/) {
		my ($cd, $ss) = split(";", $code->[2]);
		$pretty_code = $cd;
	    }
	    $pretty_code =~ s/lit\((\d+)\)/lit\(<a href='http:\/\/www\.ncbi\.nlm\.nih\.gov\/sites\/entrez\?cmd=Retrieve&db=PubMed&list_uids=$1&dopt=AbstractPlus' target=_blank>$1<\/a>\)/;
	    push(@{$pretty_codes->{$code->[0]}}, $pretty_code);
	}
	foreach my $entry (keys(%$pretty_codes)) {
	    $feature_data{$entry}{'evcodes'} = join(", ", @{$pretty_codes->{$entry}});
	}
    }

    return \%feature_data;
}

sub feature_entry {
    my($fig, $fid, $region_mid, $extended, $all_functions, $all_locations, $all_aliases) = @_;

    my $genome              = $fig->genome_of($fid);
    my $type                = $fig->ftype($fid);
    my $loc                 = $all_locations->{$fid};
    my($contig, $beg, $end) = $fig->boundaries_of($loc);
    my($left, $right)       = sort {$a <=> $b} ($beg, $end);
    my $strand              = ($beg <= $end)? '+' : '-';
    my $offset              = int(($left + $right)/2) - $region_mid;
    my $offset_beg          = $left  - $region_mid;
    my $offset_end          = $right - $region_mid;
    my $func                = $all_functions->{$fid};
    my $location            = FullLocation->new($fig, $genome, $loc);
    my $size = 0;
    map { $size += $_->Length } @{$location->Locs};

    my $retval = { 'fid'        => $fid,
		   'location'   => $loc,
		   'type'       => $type,
		   'contig'     => $contig,
		   'beg'        => $beg,
		   'end'        => $end,
		   'size'       => $size,
		   'strand'     => $strand,
		   'offset'     => $offset,
		   'offset_beg' => $offset_beg,
		   'offset_end' => $offset_end,
		   'function'   => $func };
    if ($extended) {
        my $a = $all_aliases->{$fid};
	my $aliases = ref($a) ? join("|", @$a) : '';

	$retval->{aliases} = $aliases;
    }

    return $retval;
}

sub add_functional_coupling {
    my($fig, $pin_desc, $regions, $feature_data) = @_;

    my $fc_cutoff = defined($pin_desc->{'fc_cutoff'})? $pin_desc->{'fc_cutoff'} : 4;

    foreach my $region ( @$regions )
    {
        my $pin  = $region->{'pinned_peg'};
        my @pegs = grep {$feature_data->{$_}{'type'} eq 'peg'} @{ $region->{'features'} };

        foreach my $couple ( $fig->coupled_to_batch(@pegs) )
        {
	    next unless $couple;
            my($peg1, $peg2, $sc) = @$couple;
            
            if ( $peg1 eq $pin and $sc >= $fc_cutoff )
            {
                $feature_data->{$peg2}{'fc_score'} = $sc;
		$feature_data->{$peg2}{'fc_pin'} = $peg1;
            }
        }
    }
}

sub add_figfams {
    return;
#     my($fig, $feature_data) = @_;

#     # Get FigFams directory from config file
#     my $figfam_dir = $fig->get_figfams_data();

#     # Check if FigFams directory is defined and exists on current machine
#     if ( defined($figfam_dir) and (-d $figfam_dir) )
#     {
#         # Get all PEG IDs
#         my @pegs            = grep {$feature_data->{$_}{'type'} eq 'peg'} keys %$feature_data;
#         # Get $figfams object
#         my $figfams         = new FigFams($fig, $figfam_dir);
#         # Get FigFam family ID for @pegs
#         my $figfam          = $figfams->families_containing_peg_bulk(\@pegs);
#         # Get hash of FigFam ID to family function
#         my $family_function = $figfams->family_functions();

#         # Some FigFams (the ones that are not subsystem related) have the FigFam ID in the family function.
#         # This results in the popup displaying the ID twice, so one should be removed.
        
#         my %figfam_text;
#         foreach my $fid  ( keys %$feature_data )
#         {
#             if ( $figfam->{$fid} )
#             {
#                 my $figfam_id = $figfam->{$fid};
#                 if ( ! exists $figfam_text{$figfam_id} ) {
#                     if ( $family_function->{$figfam_id} =~ /^FIG\d+/ ) {
#                         $figfam_text{$figfam_id} = $family_function->{$figfam_id};
#                     } else {
#                         $figfam_text{$figfam_id} = $figfam_id . ': ' . $family_function->{$figfam_id};
#                     }
#                 }

#                 # Add FigFam information to hash -- to go into the popup text
#                 $feature_data->{$fid}{'figfam'} = $figfam_text{$figfam_id};
#             }
#         }
#     }
}

sub add_pattyfams {
    my($fig, $feature_data) = @_;

    my @pegs = grep {$feature_data->{$_}{'type'} eq 'peg'} keys %$feature_data;
    my $fams = $fig->pattyfams_for_proteins(\@pegs);
    foreach my $fid  ( keys %$feature_data )
    {
	my $entlist = $fams->{$fid};
	next unless ref($entlist);
	my $fdata = $feature_data->{$fid};
	for my $ent (@$entlist)
	{
	    my($fam, $fun) = @$ent;
	    $fun //= "";
	    my $istart = my $iend = "";
#	    print STDERR Dumper(FUN => $fun, $fdata);
	    if (defined($fun) && defined($fdata->{function}) && ($fun ne $fdata->{function}))
	    {
		$istart = "<i>";
		$iend = "</i>";
	    }
	    if ($fam =~ /^PLF/)
	    {
		$fdata->{'plfam'} = "$fam: $istart$fun$iend";
	    }
	    elsif ($fam =~ /^PGF/)
	    {
		$fdata->{'pgfam'} = "$fam: $istart$fun$iend";
             }
	}
    }
}

sub add_subsystem_data {
    my($fig, $pin_desc, $feature_data, $extended) = @_;

    # Get subsystem_information for =all= pegs

    my @fids = keys(%$feature_data);
    my %ssdata;
    
    if ($FIG_Config::use_subsystem_estimates)
    {
	%ssdata = $fig->subsystems_for_pegs_complete_estimate(\@fids);
    }
    else
    {
	%ssdata = $fig->subsystems_for_pegs_complete(\@fids);
    }
    
    my @subsystems = ();
    foreach my $p (keys(%ssdata)) {
	foreach my $entry (@{$ssdata{$p}}) {
	    push(@subsystems, [ $entry->[0], $entry->[1], $entry->[2], $p ]);
	}
    }
    unless ($extended) {
	@subsystems = grep { $fig->usable_subsystem($_->[0]) } @subsystems;
	@subsystems = grep { $fig->is_exchangable_subsystem($_->[0]) } @subsystems;
    }

    my %peg_to_ss;
    foreach my $rec ( @subsystems ) {
	my($ss_name, $role, $variant, $fid) = @$rec;
	$ss_name =~ s/_/ /g;
	
	if ( $variant eq '0' ) {
	    # no subsystem
	    if ($extended) {
		my $ss_text = "$ss_name (not yet classified for this organism)";
		$peg_to_ss{$fid}{$ss_name} = $ss_text;
	    }
	} elsif ( $variant eq '-1' or $variant eq '*-1' ) {
	    # subsystem not functional in this organism		    
	    my $ss_text = "$ss_name (classified 'not active' in this organism)";
	    $peg_to_ss{$fid}{$ss_name} = $ss_text;
	} else {
	    $peg_to_ss{$fid}{$ss_name} = 1;
	}
    }

    # Count number of occurences of each subsystem
    my %ss_count;
    foreach my $fid ( keys %peg_to_ss )
    {
	foreach my $ss ( keys %{ $peg_to_ss{$fid} } )
	{
	    $ss_count{$ss}++;
	}
    }

    # Sort subsystems based on count and create hash with unique index for each subsystem where 
    # lower indices correspond to higher number of occurences of the subsystem.
    my %ss_index;
    my $index = 0;
    foreach my $ss ( sort {$ss_count{$b} <=> $ss_count{$a} or $a cmp $b} keys %ss_count )
    {
        $ss_index{$ss} = ++$index;
    }

    # Add subsystem information for pegs in subsystems
    foreach my $fid ( keys %peg_to_ss )
    {
        my @subsystems;
	foreach my $key (keys(%{$peg_to_ss{$fid}})) {
	    if ($peg_to_ss{$fid}{$key} ne "1") {
		push(@subsystems, [ $peg_to_ss{$fid}{$key}, $ss_index{$key} ]);
	    } else {
		push(@subsystems, [ $key, $ss_index{$key} ]);
	    }
	}

        if ( @subsystems )
        {
            $feature_data->{$fid}{'subsystems'} = [sort {$a->[1] <=> $b->[1]} @subsystems];
        }
    }
}

sub color_pegs {
    my($fig, $pin_desc, $pinned_pegs, $regions, $feature_data, $fast_color, $sims_from, $cdd_trans) = @_;

    # Run blastp and get connected pegs
    my $blast_hit = {};
    my $color_sim_cutoff = $pin_desc->{'color_sim_cutoff'}; 

    my $color_sets;
    if ($fast_color == 2 || $fast_color == 3)
    {
	my $like = ($fast_color == 2) ? 'PLF%' : 'PGF%';

	my @pegs = grep {$feature_data->{$_}{'type'} && $feature_data->{$_}{'type'} =~ /peg|domain/} keys %$feature_data;

	my %sets;
	if (@pegs)
	{
	    my $in = join(",", map { "?" } @pegs);
	    # print STDERR Dumper(query => \@pegs, $in, $feature_data);
	    my $res = $fig->db_handle->SQL("SELECT fid, family FROM family_membership WHERE fid IN ($in) AND family LIKE '$like'",
					  undef, @pegs);
	    # print STDERR Dumper($res);
	    for my $ent (@$res)
	    {
		my($fid, $fam) = @$ent;
		push(@{$sets{$fam}}, $fid);
	    }
	}
	$color_sets = [values(%sets)];
	# print STDERR Dumper(color_sets => \%sets, $color_sets);
    }
    else
    {
	if ( $sims_from eq 'blast' )
	{
	    $blast_hit = &blast_hits($fig, $regions, $feature_data, $fast_color, $color_sim_cutoff, $cdd_trans);
	}
	else
	{
	    $blast_hit = &blast_hits_from_sims($fig, $regions, $feature_data, $fast_color, $color_sim_cutoff, $cdd_trans);
	}
	
	# Assign set numbers to pegs based on blast scores
	$color_sets = &partition_pegs($blast_hit);
    }

    # Sort peg sets to that a) larger sets have smaller set numbers, and
    #                       b) sets closer to the pinned peg have smaller set numbers
    $color_sets = &sort_color_sets($pin_desc, $pinned_pegs, $color_sets, $feature_data);

    # Add color set number into $feature_data 
    for (my $i = 0; $i < @$color_sets; $i++)
    {
        # Use an index that starts at 1 (matches with coloring index)
        my $n = $i + 1;
        foreach my $peg ( @{ $color_sets->[$i] } )
        {
            $feature_data->{$peg}{'set_number'} = $n;
            
            # If this is the pinned set, set the 'color' to 'red'
            if ( $n == 1 )
            {
                $feature_data->{$peg}{'color'} = 'red';
            }
        }
    }
}

sub color_pegs_by_function {
    my($fig, $pin_desc, $pinned_pegs, $regions, $feature_data, $fast_color, $sims_from) = @_;

    my %by_function;
    my $next;
    while (my($fid, $data) = each %$feature_data)
    {
	my $set = $by_function{$data->{function}};
	if (!defined($set))
	{
	    $set = [];
	    $by_function{$data->{function}} = $set;
	}
	push(@$set, $fid);
    }

    my $color_sets = [values %by_function];

    # Sort peg sets to that a) larger sets have smaller set numbers, and
    #                       b) sets closer to the pinned peg have smaller set numbers
    $color_sets = &sort_color_sets($pin_desc, $pinned_pegs, $color_sets, $feature_data);

    # Add color set number into $feature_data 
    for (my $i = 0; $i < @$color_sets; $i++)
    {
        # Use an index that starts at 1 (matches with coloring index)
        my $n = $i + 1;
        foreach my $peg ( @{ $color_sets->[$i] } )
        {
            $feature_data->{$peg}{'set_number'} = $n;
            
            # If this is the pinned set, set the 'color' to 'red'
            if ( $n == 1 )
            {
                $feature_data->{$peg}{'color'} = 'red';
            }
        }
    }
}

sub sort_color_sets {
    my($pin_desc, $pinned_pegs, $color_sets, $feature_data) = @_;

    # Find the color set containing the set of pinned pegs, i.e. the set containing the input peg.
    # If the input peg is found (or if input is a list of pegs, returned value is the array index for the
    # set.
    # If the input peg is not found, the returned value is an reference to an array containing the input peg.
    my $set = &pinned_set($pin_desc, $pinned_pegs, $color_sets);
    
    my $pinned_color_set;
    if ( $set =~ /^\d+$/ )
    {
        # Splice out the color set containing the input peg
        ($pinned_color_set) = splice(@$color_sets,$set,1);
    }
    else
    {
        $pinned_color_set = $set;
    }

    # Add offset (summed distance from 
    my @color_sets = map {[$_, &offset($_, $feature_data)]}  @$color_sets;
    # Sort the remaining color sets based on:
    # a. size (number of pegs) and
    # b. offset from mid-point of region
    @$color_sets = map {$_->[0]} 
                     sort {@{$b->[0]} <=> @{$a->[0]} or $a->[1] <=> $b->[1]} 
#                     sort {$a->[1] <=> $b->[1] or @{$b->[0]} <=> @{$a->[0]}} 
                       map {[$_, &offset($_, $feature_data)]}  @$color_sets;

    # Add the set of pinned pegs at the front of the returned list so that it gets colored red
    return [$pinned_color_set, @$color_sets];
}

sub pinned_set {
    my($pin_desc, $pinned_pegs, $color_sets) = @_;

    # Returns an integer if the input is a peg-list or if the input peg is in a set.
    # Returns the input peg (as an array reference) if the input peg is not in a set.

    if ( @{ $pin_desc->{'pegs'} } == 1 )
    {
        # Get input peg if it exists -- the set containing this peg should be colored red
        my $peg = $pin_desc->{'pegs'}[0];

        # Iterate through the color sets until you find the set containing the input peg
        for (my $i = 0; $i < @$color_sets; $i++)
        {
            foreach my $peg2 ( @{ $color_sets->[$i] } )
            {
                if ( $peg2 eq $peg )
                {
                    # Return the set index
                    return $i;
                }
            }
        }

        # The input peg may not be in a set if there is only one region or the coloring cutoff is very stringent.
        return [$peg];
    }
    else
    {
        # Input is a peg-list, which may be split into multiple sets -- the largest will be called the 
        # "pinned" set.
        my $i = 0;
        my $max_count = 0;
        my %pinned = map {$_ => 1} @$pinned_pegs;

        for (my $j = 0; $j < @$color_sets; $j++)
        {
            # count how many pinned pegs are found in each set
            my $count = scalar(grep {$pinned{$_}} @{ $color_sets->[$j] }) || 0;

            if ( $max_count < $count )
            {
                # Keep track of the set having the largest number of pinned pegs
                $max_count = $count;
                $i = $j;
            }
        }

        return $i;
    }
}

sub offset {
    my($set, $feature_data) = @_;

    my $offset;
    foreach my $peg ( @$set )
    {
        $offset += abs($feature_data->{$peg}{'offset'});
    }
    
    return sprintf("%.2f", $offset/@$set);
}

sub partition_pegs {
    my($blast_hit) = @_;

    my %seen;
    my $sets = [];

    # Iterate through the pegs with blast hits in arbitrary order
    foreach my $peg1 ( keys %$blast_hit )
    {
        # If it has not been 'seen' (and placed in a set) use it as the seed for a set
        if ( not $seen{$peg1} )
        {
            my $set = [$peg1];

            # Mark $peg1 as seen
            $seen{$peg1} = 1;

            # Grow set through transitive closure -- note that @$set keeps growing
            for (my $i = 0; $i < @$set; $i++)
            {
                $peg1 = $set->[$i];

                # Iterate through blast hits
                foreach my $peg2 ( keys %{ $blast_hit->{$peg1} } )
                {
                    # If $peg2 has not been seen, put it in the set, and mark it as seen
                    if ( not exists $seen{$peg2} )
                    {
                        push @$set, $peg2;
                        $seen{$peg2} = 1;
                    }
                }
            }

            # Add the new set to the set of sets
            push @$sets, $set;
        }
    }

    return $sets;
}

sub blast_hits {
    my($fig, $regions, $feature_data, $fast_color, $color_sim_cutoff, $cdd_trans) = @_;

    # Set blast environment variables
    $ENV{"BLASTMAT"} ||= "$FIG_Config::blastmat";
    if ($ENV{"PATH"} !~ /fig\/bin/) { $ENV{"PATH"} = "$FIG_Config::bin:" . $ENV{"PATH"}; }

    my $temp = $FIG_Config::fast_temp // $FIG_Config::temp;

    # Get amino acid sequences
    my @pegs = grep {$feature_data->{$_}{'type'} && $feature_data->{$_}{'type'} =~ /peg|domain/} keys %$feature_data;
    my $sequences = &get_peg_sequences($fig, \@pegs, $cdd_trans);

    # Write the entire set of sequences to a fasta file
    my $all_seqs_file = "$temp/all_seqs.$$.tmp.fasta";
    &write_seqs_to_file($sequences, $all_seqs_file);

    # Run formatdb
    &formatdb_file($fig, $all_seqs_file);

    my @files_to_remove = ($all_seqs_file);

    # If $fast_color == 0, the complete blast of all vs. all sequences is performed
    # Otherwise, instead of all vs. all, the sequences from a single pinned peg region
    # is blasted against all sequences. If any of the hits are better than a given cutoff
    # ($cutoff_2), the pegs hit are deemed to be 'done', i.e. it is assumed that blasting this 
    # peg will not yield additional information for forming the peg sets. These 'done' pegs
    # are then omitted from the blast when the sequences from the region they belong to 
    # is blasted.
    my %blast_hit;
    if ( $fast_color == 1)
    {
        my $cutoff_2 = $color_sim_cutoff * 1e-20;
        my %done_with;

        foreach my $region ( @$regions )
        {
            # Iterate through each region
            my $fids = $region->{'features'};

            # Build up a hash of peg sequences which are not 'done_with'
            my %region_seqs;
            foreach my $peg ( grep {$feature_data->{$_}{'type'} eq 'peg'} @$fids )
            {
                if ( $sequences->{$peg} and not $done_with{$peg} )
                {
                    $region_seqs{$peg} = $sequences->{$peg};
                }
            }

            if ( scalar keys %region_seqs )
            {
                # Write the region sequences to a file
                my $region_seqs_file = "$temp/region_seqs.$$.tmp.fasta";
                &write_seqs_to_file(\%region_seqs, $region_seqs_file);
		push(@files_to_remove, $region_seqs_file);
                
                # BLAST the region sequences against all other sequences
                my $hits = &blast_files($fig, $region_seqs_file, $all_seqs_file, $color_sim_cutoff);
                
                # Build up hash of blast hits
                foreach my $hit ( @$hits )
                {
                    my($peg1, $peg2, $sc) = (split(/\s+/, $hit))[0,1,10];
                
                    if ( $peg1 ne $peg2 )
                    {
                        $blast_hit{$peg1}{$peg2} = 1;
                        $blast_hit{$peg2}{$peg1} = 1;

                        # If the blast score is less than the chosen cutoff, mark $peg2 as 'done_with' so 
                        # that it will not be blasted with it's region sequences
                        if ( $sc <= $cutoff_2 )
                        {
                            $done_with{$peg2} = 1;
                        }
                    }
                }
            }
        }
    }
    else
    {
        # BLAST sequence file against itself
        my $hits = &blast_files($fig, $all_seqs_file, $all_seqs_file, $color_sim_cutoff);

        # Build up hash of blast hits
        foreach my $hit ( @$hits )
        {
            my($peg1, $peg2, $sc) = (split(/\s+/, $hit))[0,1,10];

            if ( $peg1 ne $peg2 )
            {
                $blast_hit{$peg1}{$peg2} = 1;
                $blast_hit{$peg2}{$peg1} = 1;
            }
        }
    }

    for my $file (@files_to_remove)
    {
	for my $suffix ('', qw(.phr .psq .pin))
	{
	    my $path = "$file$suffix";
	    unlink($path);
	}
    }
    

    return \%blast_hit;
}

sub blast_hits_from_sims {
    my($fig, $regions, $feature_data, $fast_color, $color_sim_cutoff, $cdd_trans) = @_;

    my %blast_hit;

    my $maxN       = 2000;
    my $maxP       = $color_sim_cutoff;
    my $select     = 'fig';
    my $max_expand = '';
    my $filters    = '';    

    # Create a hash of all pegs
    my %region_peg = map {$_ => 1} grep {$feature_data->{$_}{'type'} eq 'peg'} keys %$feature_data;

    if ( $fast_color == 1 )
    {
        my $cutoff_2 = $color_sim_cutoff * 1e-20;
        my %done_with;

        # Iterate through each peg
        foreach my $peg1 ( keys %region_peg )
        {
            # Skip the 'done_with' pegs
            next if ($done_with{$peg1});

            my @sims = $fig->sims($peg1, $maxN, $maxP, $select, $max_expand, $filters);

            foreach my $sim ( grep {$region_peg{$_->[1]} and $_->[1] ne $peg1} @sims )
            {
                my($peg2, $sc) = @$sim[1,10];

                $blast_hit{$peg1}{$peg2} = 1;
                $blast_hit{$peg2}{$peg1} = 1;
                
                # If the blast score is less than the chosen cutoff, mark $peg2 as 'done_with' so 
                # that it will not be blasted with it's region sequences
                if ( $sc <= $cutoff_2 )
                {
                    $done_with{$peg2} = 1;
                }
            }
        }
    }
    else
    {
        # Iterate through each peg
        foreach my $peg1 ( keys %region_peg )
        {
            my @sims = $fig->sims($peg1, $maxN, $maxP, $select, $max_expand, $filters);

            foreach my $peg2 ( map {$_->[1]} grep {$region_peg{$_->[1]} and $_->[1] ne $peg1} @sims )
            {
                $blast_hit{$peg1}{$peg2} = 1;
                $blast_hit{$peg2}{$peg1} = 1;
            }
        }
    }

    return \%blast_hit;
}

sub blast_files {
    my($fig, $input, $database, $cutoff) = @_;

    my $trouble = 0;
    if (!-s $input) {
	$trouble = 1;
	warn "ERROR: query file \'$input\' is empty";
    }
    
    if (!-s $database) {
	$trouble = 1;
	warn "ERROR: database file \'$database\' is empty";
    }
    
    if ($trouble) { return []; }
    
    my $cmd  = "$FIG_Config::ext_bin/blastall";
    my @args = ('-p', 'blastp',  '-i', $input, '-d', $database, '-m', 8, '-e', $cutoff, '-F', 'F');
    my @blast_out = $fig->run_gathering_output($cmd, @args);

    return \@blast_out;
}

sub formatdb_file {
    my($fig, $file) = @_;

    my $cmd = "$FIG_Config::ext_bin/formatdb -i $file -p T";
    $fig->run($cmd);
}

sub write_seqs_to_file {
    my($seq, $fasta_file) = @_;

    open(FASTA, ">$fasta_file") or die "could not create FASTA file '$fasta_file': $!";
    foreach my $peg ( keys %$seq )
    {
        print FASTA ">$peg\n$seq->{$peg}\n";
    }
    close(FASTA) or die "could not close file FASTA file '$fasta_file': $!";
}

sub get_peg_sequences {
    my($fig, $pegs, $cdd_trans) = @_;
    my %sequences;

    # Iterate over each peg
    my %lookup;
    foreach my $peg ( @$pegs )
    {
	my $seq = $cdd_trans->{$peg};

        if ( $seq )
        {
            $sequences{$peg} = $seq;
        }
        else
        {
	    $lookup{$peg} = 1;
	}
    }

    my $res = $fig->get_translation_bulk([keys %lookup]);
    #print STDERR Dumper(GFS => \@lookup, $res, $fig);
    for my $id (keys %$res)
    {
	$sequences{$id} = $res->{$id};
	delete $lookup{$id};
    }

    #
    # For FIGM, we need this.
    #
    for my $peg (keys %lookup)
    {
	$sequences{$peg} = $fig->get_translation($peg);
    }

    return \%sequences;
}

sub add_cdd
{
    my($regions, $fig, $cds, $fids_for_cds, $feature_data) = @_;

    return unless ref($fids_for_cds) eq 'ARRAY';
    my %genomes = map { $fig->genome_of($_) => 1 } @$fids_for_cds;
    my %fids = map { $_ => 1 } @$fids_for_cds;

    my @new;
    my $trans = {};
    
    for my $row (@$regions)
    {
	push(@new, $row);
	next unless $genomes{$row->{genome_id}};

	#
	# new row with cdd data.
	#

	my $feats = [];
	my $new = {
	    beg => $row->{beg},
	    mid => $row->{mid},
	    end => $row->{end},
	    contig => $row->{contig},
	    contig_length => $row->{contig_length},
	    pinned_peg_strand => $row->{pinned_peg_strand},
	    org_name => "$row->{org_name} CDD",
	    genome_id => $row->{genome_id},
	    features => $feats,
	};
	for my $fid (grep { $fids{$_} } @{$row->{features}})
	{
	    my @x = $cds->create_cdd_features($fid, { data_mode => 'rep', cached_only => 0 });
	    for my $f (@x)
	    {
		my($cfid, $type, $canno, $cloc, $ctrans) = @$f;
		push(@$feats, $cfid);
		$trans->{$cfid} = $ctrans;

		my $genome              = $row->{genome_id};
		my $loc                 = $cloc;
		my($contig, $beg, $end) = $fig->boundaries_of($loc);
		my($left, $right)       = sort {$a <=> $b} ($beg, $end);
		my $strand              = ($beg <= $end)? '+' : '-';
		my $offset              = int(($left + $right)/2) - $row->{mid};
		my $offset_beg          = $left  - $row->{mid};
		my $offset_end          = $right - $row->{mid};
		my $func                = $canno;
		my $location            = FullLocation->new($fig, $genome, $loc);
		my $size = 0;
		map { $size += $_->Length } @{$location->Locs};
		
		$feature_data->{$cfid} = {
		    'parent'     => $fid,
		    'fid'        => $cfid,
		    'location'   => $loc,
		    'type'       => $type,
		    'contig'     => $contig,
		    'beg'        => $beg,
		    'end'        => $end,
		    'size'       => $size,
		    'strand'     => $strand,
		    'offset'     => $offset,
		    'offset_beg' => $offset_beg,
		    'offset_end' => $offset_end,
		    'function'   => $func
		    };		    
	    }
	}
	push(@new, $new);
    }
    @$regions = @new;
    return $trans;
}

1;
