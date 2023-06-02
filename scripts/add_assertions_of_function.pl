# -*- perl -*-
#
# Copyright (c) 2003-2013 University of Chicago and Fellowship
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


# usage: add_assertions_of_function [G1 G2 G3 ...]

use strict;
use FIG;
use Tracer;
my $fig = new FIG;
my $dbf = $fig->db_handle;

my($temp_dir,$user,$prot_id,$quality,$made_by,$file,$line,$assigned_function,$org,@pieces,$piece,@terms,$term);
my($remove_dups,%assignments,$must_correct);

# Get the genome list. If all genomes are requested, $mode will be set to 'all'.
my ($mode, @genomes) = FIG::parse_genome_args(@ARGV);

Trace("Clearing old data.") if T(2);
if ($mode eq 'all') {
    $dbf->drop_table( tbl => "assigned_functions" );
    $dbf->create_table( tbl  => "assigned_functions",
            flds => qq(prot varchar(64),
		       made_by varchar(32),
		       assigned_function text,
		       quality char,
		       org varchar(64),
		       annotation_written char(1),
		       last_modification timestamp
		       )
            );

    $dbf->drop_table( tbl => "assigned_functions_log" ) if $dbf->table_exists( "assigned_functions_log" );
    # We now recreate it below, outside of the condition

    $dbf->drop_table( tbl => "roles" );
    $dbf->create_table( tbl  => "roles",
                        flds => "prot varchar(64), role varchar(1024),"  .
                                "made_by varchar(32), org varchar(64)"
            );
}
else
{
    my $genome;
    foreach $genome (@ARGV)
    {
        $dbf->SQL("DELETE FROM assigned_functions WHERE ( org = ? )", undef, $genome);
        $dbf->SQL("DELETE FROM roles WHERE ( org = ? )", undef, $genome);
    }
}

#  This table did not exist on my system, so let's test and create if necessary -- GJO
$dbf->table_exists( "assigned_functions_log" )
    or $dbf->create_table( tbl  => "assigned_functions_log",
            flds => qq(prot varchar(64),
		       made_by varchar(32),
		       assigned_function text,
		       quality char,
		       org varchar(64),
		       created timestamp default current_timestamp
		       )
            );

$temp_dir = $FIG_Config::temp;
Trace("Finding assignments.") if T(2);
foreach $_ (&files_with_assignments($mode, @genomes))
{
    #   The third value (remove_dups) indicates whether or not updates may have been
    #   appended to the end of the file.
    ($file,$made_by,$remove_dups) = @$_;
    Trace("Processing $file for user $made_by: remove_dups = $remove_dups") if T(3);
    $made_by =~ s/\s/_/g;

    if ($remove_dups) {
        undef %assignments;
        $must_correct = 0;
        Open(\*TMP, "<$file");
        while (defined($line = <TMP>))
        {
            if ($line =~ /^(\S+)\t(\S.*\S)\s*$/)
            {
                if ($assignments{$1})
                {
                    $must_correct = 1;
                }
                $assignments{$1} = $2;
            }
        }
        close(TMP);
    
        if ($must_correct)
        {
            Trace("Removing duplicates for $file.") if T(3);
            unlink("$file~");
            rename($file,"$file~") || Confess("could not rename $file");
            Open(\*TMP, ">$file");
            foreach $prot_id (sort { &FIG::by_fig_id($a,$b) } keys(%assignments))
            {
                print TMP "$prot_id\t$assignments{$prot_id}\n";
            }
            close(TMP);
            chmod(0777,$file);
        }
    }
    # Now $file contains the assignments with all the duplicates removed.
    my  $aN = 0;
    my  $aLN = 0;
    my  $rN = 0;

    if (open(TMP, "<$file"))
    {
        Open(\*ASSIGNMENTS, ">$temp_dir/tmp$$");
        Open(\*ASSIGNMENT_LOGS, ">$temp_dir/tmpl$$");
        Open(\*INDEX, "| sort -u -T $temp_dir > $temp_dir/tmpgen$$");
    
        while (defined($line = <TMP>))
        {
            chomp $line;
            ($prot_id,$assigned_function,$quality) = split(/\t/,$line);
            $assigned_function =~ s/^\s+//;
            $assigned_function =~ s/(\t\S)?\s*$//;
    
            if (($prot_id !~ /^fig\|/) && $quality) { ($assigned_function,$quality) = ("$assigned_function $quality","") }
            if ($prot_id =~ /^pir/) { $assigned_function =~ s/\s+\[imported\]\s*$// }
    
            next if (! $prot_id);
    
            $org = $fig->genome_of($prot_id);
            $org = $org ? $org : "unknown";
    
            $assigned_function =~ s/\\$//;    #...Backslashes appear to cause problems for PostGres...
            $assigned_function =~ s/\\/ /g;   #...Backslashes appear to cause problems for PostGres...
    
            if (defined($prot_id) && defined($assigned_function) && 
                ((! $quality) || (length($quality) == 1)))
            {
                $quality = $quality ? $quality : "";
                if (&verify_row($prot_id,$made_by,$assigned_function,$quality,$org))
                {
                    &add_assignment_to_db($dbf,"$temp_dir/tmp$$",\$aN,
                        "$prot_id\t$made_by\t$assigned_function\t$quality\t$org\tL\t\\N\n");
                    &add_assignment_log_to_db($dbf,"$temp_dir/tmpl$$",\$aLN,
                        "$prot_id\t$made_by\t$assigned_function\t$quality\t$org\t\\N\n");
                    if ($prot_id =~ /^fig\|/) {
                        #  The original code also had a minimum length on roles, but ...
                        my @roles = grep { /\S/ && length $_ <= 1024 }
                                    FIG::roles_of_function( $assigned_function );
                        foreach my $role (@roles)
                        {
                            &add_role_to_db($dbf,"$temp_dir/tmpgen$$",\$rN,"$prot_id\t$role\t$made_by\t$org\n");
                        }
                    }
                }
                elsif ( length($quality) > 1 )
                {
                    print STDERR "Invalid quality score in data: $quality.\n";
                }
            }
        }
        close(TMP);
        close(ASSIGNMENTS);
        close(ASSIGNMENT_LOGS);
        close(INDEX);
    
        if ($aN > 0)
        {
            $dbf->load_table( tbl => "assigned_functions",
                      file => "$temp_dir/tmp$$" );
	}
        if ($aLN > 0)
        {
	    eval {
		$dbf->load_table( tbl => "assigned_functions_log",
			  file => "$temp_dir/tmpl$$" );
		}; 
	    if ($@)
	    {
		    die "assigned_functions_log load failed, file $temp_dir/tmpl$$\n$@";
		}
        }
        if ($rN > 0)
        {
            $dbf->load_table( tbl => "roles",
                      file => "$temp_dir/tmpgen$$" );
        }
    }
    else
    {
        Trace("Could not open $file.") if T(0);
    }
}

if ($mode eq 'all')
{
    Trace("Creating assignments index.") if T(2);
    $dbf->create_index( idx  => "assignments_ix",
            tbl  => "assigned_functions",
            type => "btree",
            flds => "prot,made_by" );
    $dbf->create_index( idx  => "assignments_ix2",
            tbl  => "assigned_functions",
            type => "btree",
            flds => "org" );
    $dbf->create_index( idx  => "assignments_log_ix",
            tbl  => "assigned_functions_log",
            type => "btree",
            flds => "prot,made_by" );
    Trace("Creating roles index.") if T(2);
    $dbf->create_index( idx  => "roles_ix",
            tbl  => "roles",
            type => "btree",
            flds => "role,org" );
    Trace("Creating protein index.") if T(2);
    $dbf->create_index( idx  => "roles_ix2",
            tbl  => "roles",
            type => "btree",
            flds => "prot" );
    Trace("Creating org index.") if T(2);
    $dbf->create_index( idx  => "roles_ix3",
            tbl  => "roles",
            type => "btree",
            flds => "org" );
    Trace("Vaccuuming tables.") if T(2);
    $dbf->vacuum_it("assigned_functions");
    $dbf->vacuum_it("assigned_functions_log");
    $dbf->vacuum_it("roles");
}

unlink("$temp_dir/tmp$$","$temp_dir/tmpgen$$");
undef $fig;
Trace("Function assertions added.") if T(2);

sub files_with_assignments {
    my($mode,@genomes) = @_;
    my(@files,$genome,@users,$user) ;

    @files = ();
    if ($mode eq "all")
    {
        @files = (["$FIG_Config::global/ext_func.table","master",0]);
        opendir(ORG,"$FIG_Config::organisms") || die "Where are the organisms?";
        @genomes = sort {$a <=> $b} grep { $_ =~ /^\d+\.\d+$/ } readdir(ORG);
        closedir(ORG);
    }

    foreach $genome (@genomes)
    {
        if (-s "$FIG_Config::organisms/$genome/assigned_functions")
        {
            push(@files,["$FIG_Config::organisms/$genome/assigned_functions","master",1]);
        }
    
        if (0 &&   ##### I am turning off user models as of June 1, 2006 [RAO]
	    (-d "$FIG_Config::organisms/$genome/UserModels") &&
            opendir(USERS,"$FIG_Config::organisms/$genome/UserModels"))
        {
            @users = grep { $_ !~ /^\./ } readdir(USERS);
            closedir(USERS);
            foreach $user (@users)
            {
                if (-s "$FIG_Config::organisms/$genome/UserModels/$user/assigned_functions")
                {
                    push(@files,["$FIG_Config::organisms/$genome/UserModels/$user/assigned_functions",$user,1]);
                }
            }
        }
    }
    return @files;
}

# The idea here is we build up the assignment file, and if it exceeds a certain
# size, we copy it into the database.
sub add_assignment_to_db {
    my($dbf, $file, $aNP, $row) = @_;
    if ($$aNP > 50000)
    {
        my @tmp = `date`;
        Trace("Copying assignments: $tmp[0]") if T(4);
        close(ASSIGNMENTS);
    
        $dbf->load_table( tbl => "assigned_functions",
                  file => "$temp_dir/tmp$$" );
    
        $$aNP = 0;
        Open(\*ASSIGNMENTS, ">$file");
    }
    print ASSIGNMENTS $row;
    $$aNP++;
}

sub add_assignment_log_to_db {
    my($dbf, $file, $aNP, $row) = @_;
    if ($$aNP > 50000)
    {
        my @tmp = `date`;
        Trace("Copying assignments: $tmp[0]") if T(4);
        close(ASSIGNMENT_LOGS);
    
        $dbf->load_table( tbl => "assigned_functions_log",
                  file => "$temp_dir/tmp$$" );
    
        $$aNP = 0;
        Open(\*ASSIGNMENT_LOGS, ">$file");
    }
    print ASSIGNMENT_LOGS $row;
    $$aNP++;
}

sub add_role_to_db {
    my($dbf,$file,$rNP,$row) = @_;

    if ($$rNP > 50000)
    {
        Trace("Copying role index for $file.") if T(3);
        close(INDEX);
    
        $dbf->load_table( tbl => "roles",
                  file => "$temp_dir/tmpgen$$" );
        $$rNP = 0;
        Open(\*INDEX, "| sort -u -T $temp_dir > $file");
    }
    print INDEX $row;
    $$rNP++;
}

sub verify_row {
    my($prot_id,$made_by,$assigned_function,$quality,$org) = @_;

    if ((length($prot_id) <= 64) &&
        (length($made_by) <= 32) &&
        (length($quality) <= 1) &&
        (length($org) <= 64))
    {
        return 1;
    }

    Trace("Field-width overflow in entry: \"$prot_id\t$made_by\t$assigned_function\t$quality\t$org\"") if T(0);
    return 0;
}
