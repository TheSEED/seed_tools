#!/usr/bin/perl -w

package Dlits;

    use strict;
    use Tracer;
    use FIG_Config;
    use XML::LibXML;
    use LWP::UserAgent;

=head1 DLIT Manager

=head2 Introduction

This module contains replacement methods for the FIG.pm routines that manipulate
literature references. Literature references were formerly stored as evidence
code attributes for features in a special database table; now they will be stored
as evidence code attributes for proteins in the L<CustomAttributes> database.

These methods are called directly from the L<FIG> module with the object itself
as the first parameter. This is not technically an object-oriented package;
however, its methods are functionally equivalent to object-oriented methods
for the FIG package itself. They are separated out in order to reduce the code
bloat in FIG.pm, which is sufficient to crash the syntax highlighting in some
editors.

The publication data is still stored in a database table and has not changed.

=cut

################################# Dlit Stuff       ###################################

=head3 add_dlit

    $rc = $fig->add_dlit(
                          -status   => 'D',       # required
                          -peg      => $peg,      # or -md5 => $md5,  # one is required
                          -pubmed   => $pubmed,   # required
                          -curator  => 'RossO',   # required
                          -override => 1);        # default = 0

This adds a dlit tuple.

=over 4

=item -status

Must be C<D>: anything else is ignored.

=item -md5

ID of the protein to which the literature reference is attached.

=item -peg

ID of the PEG feature to which the literature reference is attached. This parameter
is mutually exclusive with C<-md5>, and is provided as a convenience. The MD5 hash
of the PEG's protein will be computed as used instead.

=item -pubmed

PUBMED ID for the literature in question (all numeric, but stored as string)

=item -curator

Name of the curator making the assertion (30 char max).

=item -override (optional)

If TRUE, then an existing literature reference will be replaced with this new
information. Otherwise, an existing literature reference will cause this
add to be ignored. An I<existing literature reference> in this case is a
connection between the same protein and PUBMED ID. Replacement means changing
the curator.

=item RETURN

Returns TRUE if successful, FALSE if the literature reference was not added.

=cut

sub add_dlit {
    my( $self, @parms ) = @_;
    my %parms = @parms;        #  Previous code clobbered the defaults
    $parms{-override} ||=  0;  #  Moved default here

    #  Check for required parameters
    return 0 if ! $parms{-status} || $parms{-status} ne 'D';
    return 0 if ! ( $parms{-peg} || $parms{-md5} );
    return 0 if ! $parms{-pubmed};
    return 0 if ! $parms{-curator};

    my $status   = $parms{-status};
    my $peg      = $parms{-peg};
    my $md5      = $peg ? $self->md5_of_peg($peg) : lc $parms{-md5};
    my $pubmed   = $parms{-pubmed};
    my $curator  = $parms{-curator};
       $curator  =~ s/^master://i;     # Strip master from the recorded curator
    my $override = $parms{-override};  # Moved here to collect initializations
    # Compute the attribute ID for this protein.
    my $protID = "Protein:$md5";
    # Assume we've failed unless we decide otherwise.
    my ($rc, $action) = (0, "insert");

    # Get all of the evidence codes for this protein.
    my @codes = $self->get_attributes($protID, "evidence_code");
    # Look for a duplicate.
    my @dups = grep { $_ =~ /dlit\($pubmed\)/ } @codes;
    if (@dups) {
        # Here a duplicate code was found. Check the override mode.
        if ($override) {
            # Overriding, so we delete the old codes and set the action to
            # replace.
            $action = "replace";
            for my $dup (@dups) {
                $self->delete_attribute($protID, "evidence_code", $dup);
            }
        } else {
            # Not overriding, so kill the whole operation.
            $action = "ignore";
        }
    }
    # Only proceed if we're NOT ignoring the update.
    if ($action ne "ignore") {
        # Insert the new DLIT code.
        $self->add_attribute($protID, "evidence_code", "dlit($pubmed);$curator");
        # Denote we've been successful.
        $rc = 1;
        # Add logging
        &FIG::verify_dir( "$FIG_Config::data/Dlits" );
        if ( open LOG, ">>$FIG_Config::data/Dlits/dlits.log" ) {
            print LOG join( "\t", $action, $status, $md5, $pubmed, $curator ), "\n";
            close LOG;
        }
    }

    return $rc;
}


=head3 add_title

    $rc = $fig->add_title( $pubmed_id, $title )

Add a pubmed title to the database.  If the pubmed_id is not already
present, the id and title are added.  The return code reflects that success
or failure of the add.  If the pubmed_id is already defined,
and the titles match, there is no change, and the return code is 2.
If the id exists and the title is different, no change is made, and
the return code is 0.  To change an existing title, use:

    $rc = $fig->update_title( $pubmed_id, $title )

The returned values are:

    0  attempting to change a title, or failure;
    1  successful addition of a new title; or
    2  existing and new titles are the same

=cut

sub add_title {
    my( $self, $pubmed, $title ) = @_;

    my $rdbH  = $self->db_handle;

    #  If there is already a title, do not duplicate it.

    my $db_resp = $rdbH->SQL( "SELECT title
                               FROM pubmed_titles
                               WHERE ( pubmed = '$pubmed' )"
                            );

    #  Same title is success; different title is failure.

    if ( $db_resp && @$db_resp )
    {
        return $db_resp->[0]->[0] eq $title ? 2 : 0;  # Same title is success
    }

    #  If it does not exist, add it
    $title = quotemeta $title;
    return $rdbH->SQL( "INSERT
                        INTO pubmed_titles ( pubmed, title )
                        VALUES ( ?, ? )",
                        undef, $pubmed, $title
                     );
}


=head3 update_title

    $rc = $fig->update_title( $pubmed_id, $title )

Add or change a pubmed title to the database.  If the pubmed_id is not already
present, the id and title are added.  The return code reflects that success
or failure of the add.  If the pubmed_id is already defined,
and the titles match, there is no change, and the return code is 2.
If the id exists and the title is different, no change is made, and
the return code is 0.  To change an existing title, use:

    $rc = $fig->update_title( $pubmed_id, $title )

The returned values are:

    0  on failure;
    1  successful addition or change of a title; or
    2  existing and new titles are the same

=cut

sub update_title {
    my( $self, $pubmed, $title ) = @_;

    my $rdbH  = $self->db_handle;

    #  If there is already a title, do not duplicate it.

    my $rc = $rdbH->SQL( "SELECT title
                          FROM pubmed_titles
                          WHERE ( pubmed = '$pubmed' )"
                       );

    #  If there is a title, make sure that it is what we want.
    $title = quotemeta $title;
    if ( $rc && @$rc )
    {
        return 2 if $rc->[0]->[0] eq $title;     # Same title is success

        # title is different, change it

        return $rdbH->SQL( "UPDATE pubmed_titles
                            SET title = '$title'
                            WHERE ( pubmed = ? )",
                            undef, $pubmed
                         );
    }

    return $rdbH->SQL( "INSERT
                        INTO pubmed_titles ( pubmed, title )
                        VALUES ( ?, ? )",
                        undef, $pubmed, $title
                     );
}


=head3 get_title

    $title = $fig->get_title( $pubmed_id )

Get a title for a literature id

Returned value:

    $title   upon success
    undef    upon failure

=cut

sub get_title {
    my( $self, $pubmed ) = @_;

    my $rdbH  = $self->db_handle;

    #  If there is already a title, do not duplicate it.

    my $db_resp = $rdbH->SQL( "SELECT title
                               FROM pubmed_titles
                               WHERE ( pubmed = '$pubmed' )"
                            );

    return ( $db_resp && @$db_resp ) ? $db_resp->[0]->[0] : undef;
}


=head3 all_titles

    [ [ id, title ], ... ] = $fig->all_titles()

Get all pubmed_id, title pairs

Returned value:

    [ [ id, title ], ... ]   upon success
    []                       upon failure

=cut

sub all_titles {
    my( $self ) = @_;
    my $rdbH  = $self->db_handle;

    my $db_resp = $rdbH->SQL( "SELECT DISTINCT * FROM pubmed_titles" );
    return $db_resp ? [ sort { $a->[0] <=> $b->[0] }  @$db_resp ] : [];
}


=head3 all_dlits

    $dlits = $fig->all_dlits();
    @dlits = $fig->all_dlits();

Returns a reference to an array, or a list, of all current dlit data.

The returned value is

    [ [ status, md5_hash, pubmed, curator ], ... ]

or

    ( [ status, md5_hash, pubmed, curator ], ... )
=cut

# This function was called by FIG.pm, but was missing in this module

sub all_dlits {
    my $self  = shift or return wantarray ? () : [];
    my @dlits = map { $_->[1] =~ s/^Protein://; $_ }
                map { $_->[2] =~ m/^dlit\((\d+)\)(;(\S+))?/ ? [ 'D', $_->[0], $1, $3 || '' ] : () }
                $self->get_attributes( undef, 'evidence_code', 'dlit%' );
    wantarray ? @dlits : \@dlits;
}


=head3 get_dlits_for_peg

    $dlits = $fig->get_dlits_for_peg( $fid );

Returns a reference to an array of current dlit data for a peg.

The returned value is

    [ [ status, md5_hash, pubmed, curator ], ... ]
=cut

sub get_dlits_for_peg {
    my ( $self, $peg ) = @_;
    return wantarray ? () : [] if ! $peg;  # Slightly less antisocial than undef

    my $md5 = $self->md5_of_peg( $peg );
    $md5 = $peg if ! $md5 && $peg =~ /^[0-9A-Fa-f]{32}$/;
    return wantarray ? () : [] if ! $md5;

    my @codes = map { $_->[2] =~ /^dlit\((\d+)\)(;(.+))?/ ? ['D', $md5, $1, $3] : () }
                $self->get_attributes( "Protein:$md5", 'evidence_code', 'dlit%' );

    wantarray ? @codes : \@codes;
}

=head3 get_dlits_for_pegs

    $dlits = $fig->get_dlits_for_pegs( \@fids );
    @dlits = $fig->get_dlits_for_pegs( \@fids );

Returns a reference to an array of current dlit data for a list of pegs.
This interface consolidates queries for proteins of identical sequence.
For very complex sets, it might be worth requesting all dlits, and filtering
the output.

The returned value is

    [ [ status, fid, pubmed, curator ], ... ]

or

    ( [ status, fid, pubmed, curator ], ... )

Note that this differs from the output of get_dlits_for_peg(), which
returns the md5 value as the protein identifier.
=cut

sub get_dlits_for_pegs {
    my ( $self, $ids ) = @_;

    return wantarray ? () : [] unless $self && $ids;
    $ids = [ $ids ] unless ref $ids eq 'ARRAY';

    my %md5s;
    foreach my $id ( @$ids )
    {
        my $md5 = $self->md5_of_peg( $id ) || ( $id =~ /^[0-9A-Fa-f]{32}$/ ? $id : '' );
        push @{ $md5s{ $md5 } }, $id if $md5;
    }
    return wantarray ? () : [] unless keys %md5s;
    my @keys = map { "Protein:$_" } keys %md5s;

    my @codes;
    foreach ( $self->get_attributes( \@keys, 'evidence_code', 'dlit%' ) )
    {
        my ( $pmid, undef, $curator ) = $_->[2] =~ /^dlit\((\d+)\)(;(.+))?/;
        next unless $pmid;
        my ( $md5 ) = $_->[0] =~ /^Protein:(.*)$/;
        # Convert the md5 values back to the requested peg fids.
        foreach my $peg ( @{ $md5s{ $md5 } } )
        {
            push @codes, ['D', $peg, $pmid, $curator ];
        }
    }

    wantarray ? @codes : \@codes;
}

=head3 protein_evidence_codes

    my @attrs = Dlits::protein_evidence_codes($fig, $peg);

Return a list of the evidence codes for the specified PEG that are associated with
its protein sequence. This method is called from "get_attributes" to create the
illusion that the evidence codes for the proteins are still associated with features.

=over 4

=item fig

L<FIG> object used to get the protein sequence and its attributes.

=item peg

ID of the PEG feature whose evidence codes are desired.

=item RETURN

Returns a list of 3-tuples, each consisting of (0) the incoming feature ID, (1) a
constant string C<evidence_code>, and (2) an evidence code inherited from the
feature's protein sequence.

=back

=cut

sub protein_evidence_codes {
    # Get the parameters.
    my ($fig, $peg) = @_;
    # Compute the protein ID for the feature.
    my $md5 = $fig->md5_of_peg($peg);
    return () unless $md5;
    my $protID = "Protein:$md5";
    # Get its list of evidence codes.
    my @retVal = $fig->get_attributes($protID, "evidence_code");
    # Replace the protein IDs with the feature ID.
    for my $tuple (@retVal) {
        $tuple->[0] = $peg;
    }
    # Return the result.
    return @retVal;
}
#
# Routines for manipulating the SeedGlobal document mapping code.
#

sub get_mapped_documents_for_feature
{
    my($fig, $fid, $include_document_details) = @_;

    #
    # For protein features the identifier used is the MD5 of the translation.#
    #

    if (FIG::ftype($fid) eq 'peg')
    {
	$fid = $fig->md5_of_peg($fid);
    }

    my $res = $fig->seed_global_dbh->SQL_returning_hash(qq(SELECT pmid, added_by, added_on, status
							   FROM mapped_document
							   WHERE feature_id = ?), 'pmid', undef, $fid);
    my @out = values %$res;
    if ($include_document_details)
    {
	for my $ent (@out)
	{
	    my $x = get_pubmed_document_details($fig, $ent->{pmid});
	    $ent->{$_} = $x->{$_} foreach keys %$x;
	}
    }
    return @out;
}

sub get_mapped_status_types
{
    my($fig) = @_;
    my $res = $fig->seed_global_dbh->SQL(qq(SELECT status, description
					    FROM document_mapping_status));
    return $res;
}

=head3 get_pubmed_document_details

    my $docHash = Dlits::get_pubmed_document_details($figm, $pmid, $force);

Retrieve a hash of the document information for the specified pubmed
document ID.

=over 4

=item fig

A L<FIG> object or a C<DBKernel> object for the seed global database.

=item pmid

The pubmed ID of the document whose details are desired.

=item force

If TRUE, then the data will always be retrieved from PUBMED; otherwise,
it will only be retrieved if it has not been cached in the SEED global
database.

=item RETURN

Returns a reference to a hash containing the following keys and
values.

=over 8

=item title

The title of the document.

=item authors

A comma-separated list of the authors.

=item source

The source journal for the document.

=item pubdate

The publication date of the document.

=back

=back

=cut

sub get_pubmed_document_details
{
    my($fig, $pmid, $force) = @_;

    if ($pmid !~ /^\d{1,10}$/ || $pmid > 2147483647)
    {
	warn "invalid pmid $pmid\n";
	return;
    }

    my $dbh = (ref $fig eq 'FIG' ? $fig->seed_global_dbh : $fig);
    my $res = $dbh->SQL(qq(SELECT authors, title, source, ts
					    FROM pubmed_document
					    WHERE pmid = ?), undef, $pmid);
    my $exists = @$res;
    if ($exists && !$force)
    {
	my($authors, $title, $source, $pubdate) = @{$res->[0]};
	return {
	    authors => $authors,
	    title => $title,
	    source => $source,
	    pubdate => $pubdate
	};
    }

    #
    # Not in database, look up, cache results, and return.
    #

    my $entrez_base = "https://eutils.ncbi.nlm.nih.gov/entrez/";
    my $journal_url = "$entrez_base"."eutils/esummary.fcgi?db=pubmed&id=";
    my $url_format = "&retmode=xml";

    my $url = join("", $journal_url, $pmid, $url_format);

    my $ua = LWP::UserAgent->new();
    my $result = $ua->get($url);

    if ($result->is_success)
    {
	my $txt = $result->content;
	my $doc = XML::LibXML->new->parse_string($txt);

        if ($doc)
        {
            my($source) = get_entrez_item($doc, "Source");
            my($title) = get_entrez_item($doc, "Title");
            my($ref) = get_entrez_item($doc, "SO");
	    my($pubdate) = get_entrez_item($doc, "PubDate");
	    my @authors = get_entrez_item($doc, "Author");
	    my $authors = join(", ", @authors);

	    $dbh->SQL(qq(DELETE FROM pubmed_document WHERE pmid = ?), undef, $pmid) if $exists;
	    $dbh->SQL(qq(INSERT INTO pubmed_document (pmid, authors, title, source)
					  VALUES (?, ?, ?, ?)), undef, $pmid, $authors, $title, $source);
	    return {
		authors => $authors,
		title => $title,
		source => $source,
		pubdate => $pubdate,
		ref => $ref,
	    };
        }
	else
	{
	    warn "Error parsing pubmed document: $txt\n";
	}
    }
    else
    {
	warn "Failure looking up document: " . $result->content;
    }

    return undef;
}


sub get_entrez_item
{
    my($doc, $str) = @_;
    my @list = $doc->findnodes('//Item[@Name="' . $str . '"]');
    return map { $_->textContent } @list;
}

1;
