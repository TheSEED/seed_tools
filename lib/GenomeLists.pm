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

package GenomeLists;

use FIG;
use FileHandle;
use strict;

my $setdir = $FIG_Config::data . '/Sets/';
&FIG::verify_dir( $setdir );
my $listdir = $FIG_Config::data . '/Sets/Genome_Sets/';
&FIG::verify_dir( $listdir );

1;

sub new {
  
  my ( $class, $self ) = @_;
  bless( $self, $class );

  return $self;
}

sub save {

  my ( $self ) = @_;
  my $file = $listdir . $self->{ 'name' };

  open ( FILE, ">", $file ) or die "Cannot write $file: $!";

  print FILE "Name:\t";
  print FILE $self->{ 'name' };
  print FILE "\n\nDescription:\t";
  print FILE $self->{ 'description' };
  my $genomes = $self->{ 'genomes' };

  print FILE "\n\n";
  foreach my $g ( @$genomes ) {
    print FILE "Genome:\t$g\n";
  }

  close(FILE);
  return 1;
}

sub create
{
    my($name, $description, $genomes) = @_;

    my $file = $listdir . $name;
    
    open ( FILE, ">", $file ) or die "Cannot write $file: $!";

    print FILE "Name:\t";
    print FILE $name;
    print FILE "\n\nDescription:\t";
    print FILE $description;

    print FILE "\n\n";
    foreach my $g ( @$genomes ) {
	print FILE "Genome:\t$g\n";
    }
    close(FILE);
}

sub load {
  my ( $name ) = @_;

  my $file = $listdir . $name;
  open ( FILE, $file ) or return -1;
  my $self;

  my @genomes = ();
  while ( <FILE> ) {
    chomp;
    next if ( $_ eq '' );
    # get name
    if ( $_ =~ /Name:\t(.*)/ ) {
      $self->{ 'name' } = $1;
    }
    # get description
    if ( $_ =~ /Description:\t(.*)/ ) {
      $self->{ 'description' } = $1;
    }
    # get genomes
    if ( $_ =~ /Genome:\t(.*)/ ) {
      my $g = $1;
      if ( $g ne '' ) {
	push @genomes, $g;
      }
    }
  }
  $self->{ 'genomes' } = \@genomes;

  bless( $self, 'GenomeLists' );
  return $self;
}

sub getListsForUser {
  my ( $class, $user ) = @_;
  
  opendir( D, "$listdir" );
  my @filesR = grep { not /^\./ } readdir( D );
  my @files;

 return @filesR if ( !defined( $user ) );

  foreach my $f ( @filesR ) {
    if ( $f =~ /$user\_(.*)/ ) {
      push @files, $1;
    }
  }

  return @files;

}
