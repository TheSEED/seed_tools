# -*- perl -*-
########################################################################
# Copyright (c) 2003-2014 University of Chicago and Fellowship
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

package FIG;

use strict;
no warnings 'redefine';  ## prevents spurious warnings due to use recursion
use FIG_Config;

#
# See if we need to use fcntl-based file locking. If so, import
# the package and override the global definition of flock.
# This is in place at least initially for the GPFS-based install on
# the NMPDR cluster.
#

use FileLocking;
use DB_File;
use FF;

use Fcntl qw/:flock/;  # import LOCK_* constants

use POSIX;
use Devel::Backtrace;
use Scalar::Util qw(weaken);
use Errno;
use IPC::Open2;
use MIME::Base64;
use File::Basename;
use File::Copy;
use File::Path;
use File::Spec;
use FileHandle;
use DirHandle;
use IO::Socket;
use SOAP::Lite;
use LWP::UserAgent;
use LWP::Simple; # for ncbi connection - get genetic code
use Digest::MD5;
use URI::Escape;
use Time::Local;

use Carp qw(confess croak carp cluck);
use Data::Dumper;

use BasicLocation;
use gjoseqlib;
use DBrtns;
use Sim;
use Annotation;
use Blast;
use FullLocation;
use tree_utilities;
use Subsystem;
use SeedDas;
use Construct;
use FIGRules;
use Tracer;
use GenomeIDMap;
use RemoteCustomAttributes;
use FIGV;
use Dlits;
use SAPserver;
use SeedUtils;
use GenomeTypeObject;

our $haveDateParse;
eval {
    require Date::Parse;
    import Date::Parse;
    $haveDateParse = 1;
    require CustomAttributes;
    import CustomAttributes;
};

eval {
    require Net::RabbitMQ;
    require JSON::XS;
    JSON::XS->import();
};

our $haveUUID;
eval {
    require UUID;
    $haveUUID = 1;
};

eval { require FigGFF; };
if ($@ and T(1)) {
    warn $@;
}

#
# Conditionally evaluate this in case its prerequisites are not available.
#

our $ClearinghouseOK;
eval {
    require Clearinghouse;
    $ClearinghouseOK = 1;
};

#
# Try to load the RPC stuff; it might fail on older versions of the software.
#
eval {
    require FIGrpc;
};

my $xmlrpc_available = 1;
if ($@ ne "") {
    $xmlrpc_available = 0;
}

use FIGAttributes;
use base 'FIGAttributes';

use vars qw(%_FunctionAttributes);

#
# Force all new files to be all-writable.
#

umask 0;

=head1 FIG Genome Annotation System

=head2 Introduction

This is the main object for access to the SEED data store. The data store
itself is a combination of flat files and a database. The flat files can
be moved easily between systems and the database rebuilt as needed.

A reduced set of this object's functions are available via the B<SFXlate>
object. The SFXlate object uses a single database to represent all its
genomic information. It provides a much smaller capability for updating
the data, and eliminates all similarities except for bidirectional best
hits.

The key to making the FIG system work is proper configuration of the
C<FIG_Config.pm> file. This file contains names and URLs for the key
directories as well as the type and login information for the database.

FIG was designed to operate as a series of peer instances. Each instance is
updated independently by its owner, and the instances can be synchronized
using a process called a I<peer-to-peer update>. The terms
I<SEED instance> and I<peer> are used more-or-less interchangeably.

The POD documentation for this module is still in progress, and is provided
on an AS IS basis without warranty. If you have a correction and you're
not a developer, EMAIL the details to B<bruce@gigabarb.com> and I'll fold
it in.

B<NOTE>: The usage example for each method specifies whether it is static

    FIG::something

or dynamic

    $fig->something

If the method is static and has no parameters (C<FIG::something()>) it can
also be invoked dynamically. This is a general artifact of the
way PERL implements object-oriented programming.

=head2 Hiding/Caching in a FIG object

We save the DB handle, cache taxonomies, and put a few other odds and ends in the
FIG object.  We expect users to invoke these services using the object $fig constructed
using:

    use FIG;
    my $fig = new FIG;

$fig is then used as the basic mechanism for accessing FIG services.  It is, of course,
just a hash that is used to retain/cache data.  The most commonly accessed item is the
DB filehandle, which is accessed via $self->db_handle.

We cache genus/species expansions, taxonomies, distances (very crudely estimated) estimated
between genomes, and a variety of other things.

=cut


#: Constructor FIG->new();

=head2 Public Methods

=head3 new

    my $fig = FIG->new();

This is the constructor for a FIG object. It uses no parameters. If tracing
has not yet been turned on, it will be turned on here. The tracing type and
level are specified by the configuration variables C<$FIG_Config::trace_levels>
and C<$FIG_Config::trace_type>. These defaults can be overridden using the
environment variables C<Trace> and C<TraceType>, respectively.

=cut

sub new
{
    my($class) = @_;

    #
    # For the FCGI and other persistent web configurations:
    #
    # If $FIG_Config::use_fig_singleton is set,
    # we store a singleton FIG object in
    # $FIG_Config::fig_singleton and return that
    # if no FIG has yet been created. Otherwise
    # create a new FIG each time.
    #

    print STDERR "\n\nNEW FIG\n";
    print STDERR Devel::Backtrace->new;
    print STDERR "****  $FIG_Config::use_fig_singleton\n";

    if ($FIG_Config::use_fig_singleton)
    {
	if (!ref($FIG_Config::fig_singleton))
	{
	    my $fig = &_new($class);
	    $FIG_Config::fig_singleton = $fig;
	}

	return $FIG_Config::fig_singleton;
    }
    else
    {
	return &_new($class);
    }
}

sub _new {
    my($class) = @_;



    #
    # Check to see if we have a FIG_URL environment variable set.
    # If we do, don't actually create a FIG object, but rather
    # create a FIGrpc and return that as the return from this constructor.
    #
    if ($ENV{FIG_URL} && $xmlrpc_available) {
        my $figrpc = new FIGrpc($ENV{FIG_URL});
        return $figrpc;
    }
    Trace("Connecting to the database.") if T(2);
    # Connect to the database, then return ourselves.
    my $rdbH = new DBrtns;

    my $self = {
        _dbf  => $rdbH,
    };

    print STDERR "FIG Create $self $$\n";
    $self->{gdata} = GenomeDataCache->new($self);
    $self->{sdata} = SubsystemDataCache->new($self);

    if ($FIG_Config::attrOld) {
        # Use the old attribute system. This is normally only done if we
        # need to reload.
        Trace("Legacy attribute system chosen using the override feature.") if T(3);
    } elsif ($FIG_Config::attrURL) {
        Trace("Remote attribute server $FIG_Config::attrURL chosen.") if T(3);
        $self->{_ca} = RemoteCustomAttributes->new($FIG_Config::attrURL);
    } elsif ($FIG_Config::attrHost) {
        eval {
            Trace("Attribute database on $FIG_Config::attrHost chosen.") if T(3);
            my $user = ($FIG_Config::arch eq 'win' ? 'self' : "seed"); # scalar(getpwent()));
            $self->{_ca} = CustomAttributes->new(user => $user);
        };
        if ($@) {
            Tracer::Warn("Attribute server error: $@");
        }
    }
    Trace("Attribute connection complete.") if T(3);
    #
    # If we have a readonly-database defined in the config,
    # create a handle for that as well.
    #

    if (defined($FIG_Config::readonly_dbhost)) {
        my $ro = new DBrtns($FIG_Config::dbms, $FIG_Config::readonly_db, $FIG_Config::readonly_dbuser,
                            $FIG_Config::readonly_dbpass, $FIG_Config::readonly_dbport, $FIG_Config::readonly_dbhost,
                            $FIG_Config::readonly_dbsock);
        $self->{_ro_dbf} = $ro;

        #
        # Oh ick. All of the queries made go through the one dbf that a FIG holds. We want
        # to redirect the select queries through this readonly object. We'll need
        # to tell the main handle about the readonly one, and let it decide.
        #

        $rdbH->set_readonly_handle($ro);
    }

    #
    # Check for memcached.
    #

    if (ref($FIG_Config::memcached_config))
    {
	eval {
	    require Cache::Memcached::Fast;
	    $self->{memcache} = new Cache::Memcached::Fast($FIG_Config::memcached_config);
	    # print STDERR "Configured memcached\n";
	};
    }

    #
    # Check for table support for assign_functions auditing.
    #
    eval {
	my $dbh = $rdbH->{_dbh};
	local $dbh->{RaiseError} = 1;
        local $dbh->{PrintError} = 0;

	my $res = $dbh->selectall_arrayref(qq(SELECT annotation_written FROM assigned_functions LIMIT 1));
	$res = $dbh->selectall_arrayref(qq(SELECT prot FROM assigned_functions_log LIMIT 1));

	$self->{have_assignment_auditing} = 1;
    };

    return bless $self, $class;
}

=head3 SaplingCheck

    my $sapLoader = $fig->SaplingCheck();

Return the L<SaplingFunctionLoader> object for updating the Sapling
when there are assignment or annotation changes. If Sapling updates
are turned off, this method will return an undefined value.

=cut

sub SaplingCheck {
    # Get the parameters.
    my ($self) = @_;
    # Declare the return variable.
    my $retVal;
    # Are Sapling updates turned on?
    if ($FIG_Config::update_sapling) {
        # Return the Sapling updater already in place if there is on.
        $retVal = $self->{sapUpdater};
        if (! $retVal) {
            # Here we need to create the updater. Connect to the Sapling.
            require Sapling;
            my $sap = Sapling->new();
            # Create the updater.
            require SaplingFunctionLoader;
            $retVal = SaplingFunctionLoader->new($sap);
            # Cache the updater for future calls.
            $self->{sapUpdater} = $retVal;
        }
    }
    # Return the updater.
    return $retVal;
}


=head3 CacheTrick

    my $value = $fig->CacheTrick($self, $field => $evalString);

This is a helper method used to create simple field caching in another object. If the
named field is found in $self, then it will be returned directly. Otherwise, the eval
string will be executed to compute the value. The value is then cahced in the $self
object so it can be retrieved easily when needed. Use this method to make a FIG
data-access object more like an object created by PPO or ERDB.

=over 4

=item self

Hash or blessed object containing the cached fields.

=item field

Name of the field desired.

=item evalString

String that can be evaluated to compute the field value.

=item RETURN

Returns the value of the desired field.

=back

=cut

sub CacheTrick {
    # Get the parameters. Note that we get this object under the name "$fig" rather than
    # "$self", because $self represents the caller's object.
    my ($fig, $self, $field, $evalString) = @_;
    # Declare the return variable.
    my $retVal;
    # Check the cache.
    if (exists $self->{$field}) {
        # Return the cached data.
        $retVal = $self->{$field};
    } else {
        # Compute the field value.
        Trace("Retrieving data for $field using formula: $evalString") if T(4);
        $retVal = eval($evalString);
        # Cache it for future use.
        $self->{$field} = $retVal;
    }
    # Return the field value.
    return $retVal;
}

=head3 go_number_to_term

Returns GO term for GO number from go_number_to_term table in database

=cut

sub go_number_to_term {
    my($self,$id) = @_;
    my $rdbH = $self->db_handle;
    my $relational_db_response = $rdbH->SQL("SELECT go_desc  FROM go_terms where go_id = \'$id\'");
    return (@$relational_db_response == 1) ? $relational_db_response->[0]->[0] : "";
    return "";
}

sub go_number_to_info {
    my($self,$id) = @_;
    my $rdbH = $self->db_handle;
    my $relational_db_response = $rdbH->SQL("SELECT go_desc,go_type,obsolete  FROM go_terms where go_id = \'$id\'");
    return (@$relational_db_response == 1) ? $relational_db_response->[0] : "";
    return "";
}


=head3 db_handle

    my $dbh = $fig->db_handle;

Return the handle to the internal B<DBrtns> object. This allows direct access to
the database methods.

=cut

sub db_handle {
    my($self) = @_;
    $self->{_dbf}->{_dbh}->ping();
    return $self->{_dbf};
}

=head3 seed_global_dbh

    my $dbh = $fig->seed_global_dbh

Return the handle to the Seed Global database. This is a DBKernel object.

=cut

sub seed_global_dbh
{
    my($self) = @_;


    return undef;

    my $dbh = $self->{_seed_global_dbh};
    if (!$dbh)
    {
	$dbh = DBKernel->new('mysql', 'seed_global', 'seed', undef, undef,
			     'seed-db-write.mcs.anl.gov');
	$self->{_seed_global_dbh} = $dbh;
    }
    return $dbh;
}


sub table_exists {
    my($self,$table) = @_;

    my $rdbH = $self->db_handle;
    return $rdbH->table_exists($table);
}

=head3 cached

    my $x = $fig->cached($name);

Return a reference to a hash containing transient data. If no hash exists with the
specified name, create an empty one under that name and return it.

The idea behind this method is to allow clients to cache data in the FIG object for
later use. (For example, a method might cache feature data so that it can be
retrieved later without using the database.) This facility should be used sparingly,
since different clients may destroy each other's data if they use the same name.

=over 4

=item name

Name assigned to the cached data.

=item RETURN

Returns a reference to a hash that is permanently associated with the specified name.
If no such hash exists, an empty one will be created for the purpose.

=back

=cut

sub cached {
    my($self,$what) = @_;

    my $x = $self->{$what};
    if (! $x) {
        $x = $self->{$what} = {};
    }
    return $x;
}

=head3 get_system_name

    my $name = $fig->get_system_name;

Returns C<seed>, indicating that this is object is using the SEED
database. The same method on an SFXlate object will return C<sprout>.

=cut
#: Return Type $;
sub get_system_name {
    return "seed";
}

=head3 DESTROY

The destructor releases the database handle.

=cut

 sub DESTROY {
    my($self) = @_;
    my($rdbH);

    print STDERR "FIG destroy $$\n";
    if ($rdbH = $self->db_handle) {
        $rdbH->DESTROY;
    }
}

=head3 same_seqs

    my $sameFlag = FIG::same_seqs($s1, $s2);

Return TRUE if the specified protein sequences are considered equivalent and FALSE
otherwise. The sequences should be presented in I<nr-analysis> form, which is in
reverse order and upper case with the stop codon omitted.

The sequences will be considered equivalent if the shorter matches the initial
portion of the long one and is no more than 30% smaller. Since the sequences are
in nr-analysis form, the equivalent start potions means that the sequences
have the same tail. The importance of the tail is that the stop point of a PEG
is easier to find than the start point, so a same tail means that the two
sequences are equivalent except for the choice of start point.

=over 4

=item s1

First protein sequence, reversed and with the stop codon removed.

=item s2

Second protein sequence, reversed and with the stop codon removed.

=item RETURN

Returns TRUE if the two protein sequences are equivalent, else FALSE.

=back

=cut

sub same_seqs {
    my ($s1,$s2) = @_;

    my $ln1 = length($s1);
    my $ln2 = length($s2);

    return ((abs($ln1-$ln2) < (0.3 * (($ln1 < $ln2) ? $ln1 : $ln2))) &&
            ((($ln1 <= $ln2) && (index($s2,$s1) == 0)) ||
             (($ln1 > $ln2)  && (index($s1,$s2) == 0))));
}

=head3 is_locked_fid

    $fig->is_locked_fid($fid);

returns 1 iff $fid is locked

=cut

sub is_locked_fid {
    my($self,$fid) = @_;

    return 0; #### turns off recognition of locks (RAO, 4/27/2010)

    if (! $self->table_exists('fid_locks')) { return 0 }
    my $rdbH = $self->db_handle;
    my $relational_db_response = $rdbH->SQL("SELECT fid FROM fid_locks WHERE fid = \'$fid\' ");
    return (@$relational_db_response > 0) ? 1 : 0;
}

=head3 lock_fid

    $fig->lock_fid($user,$fid);

Sets a lock on annotations for $fid.

=cut

sub lock_fid {
    my($self,$user,$fid) = @_;

    if (! $self->table_exists('fid_locks'))       { return 0 }
    if ((! $user) || ($fid !~ /^fig\|\d+\.\d+/))  { return 0 }
    if ($self->is_locked_fid($fid))               { return 1 }

    my $func = $self->function_of($fid);
    $self->add_annotation($fid,$user,"locked assignments to '$func'");

    my $rdbH = $self->db_handle;
    my $relational_db_response = $rdbH->SQL("SELECT fid FROM fid_locks WHERE fid = \'$fid\' ");
    if (! (@$relational_db_response > 0))
    {
        $rdbH->SQL("INSERT INTO fid_locks ( fid ) VALUES ( '$fid' )");
	if ($fid =~ /^fig\|(\d+\.\d+)\.([^\.]+)/)
	{
	    my $genome = $1;
	    my $type = $2;
	    if (open(TMP,">>$FIG_Config::organisms/$genome/Features/$type/locks"))
	    {
		print TMP "$fid\t1\n";
	    }
	    close(TMP);
	}
    }
}

=head3 unlock_fid

    $fig->unlock_fid($user,$fid);

Sets a unlock on annotations for $fid.

=cut

sub unlock_fid {
    my($self,$user,$fid) = @_;

    if (! $self->table_exists('fid_locks'))       { return 0 }
    if ((! $user) || ($fid !~ /^fig\|\d+\.\d+/))  { return 0 }
    if (! $self->is_locked_fid($fid))             { return 1 }

    $self->add_annotation($fid,$user,"unlocked assignments");
    my $rdbH = $self->db_handle;
    my $relational_db_response = $rdbH->SQL("SELECT fid FROM fid_locks WHERE fid = '$fid' ");
    $rdbH->SQL("DELETE FROM fid_locks WHERE ( fid = '$fid' )");
    if ($fid =~ /^fig\|(\d+\.\d+)\.([^\.]+)/)
    {
	my $genome = $1;
	my $type = $2;
	if (open(TMP,">>$FIG_Config::organisms/$genome/Features/$type/locks"))
	{
	    print TMP "$fid\t0\n";
	}
	close(TMP);
    }
}

=head3 is_propagation_locked_fid

    $fig->is_propagation_locked_fid($fid);

returns 1 iff $fid is propagation_locked

=cut

sub is_propagation_locked_fid {
    my($self,$fid) = @_;

    if (! $self->table_exists('fid_propagation_locks')) { return 0 }
    my $rdbH = $self->db_handle;
    my $relational_db_response = $rdbH->SQL("SELECT fid FROM fid_propagation_locks WHERE fid = \'$fid\' ");
    return (@$relational_db_response > 0) ? 1 : 0;
}

=head3 propagation_lock_fid

    $fig->propagation_lock_fid($user,$fid);

Sets a propagation_lock on annotations for $fid.

=cut

sub propagation_lock_fid {
    my($self,$user,$fid) = @_;

    if (! $self->table_exists('fid_propagation_locks'))       { return 0 }
    if ((! $user) || ($fid !~ /^fig\|\d+\.\d+/))  { return 0 }
    if ($self->is_propagation_locked_fid($fid))               { return 1 }

    my $func = $self->function_of($fid);
    $self->add_annotation($fid,$user,"propagation_locked assignments to '$func'");

    my $rdbH = $self->db_handle;
    my $relational_db_response = $rdbH->SQL("SELECT fid FROM fid_propagation_locks WHERE fid = \'$fid\' ");
    if (! (@$relational_db_response > 0))
    {
        $rdbH->SQL("INSERT INTO fid_propagation_locks ( fid ) VALUES ( '$fid' )");
	if ($fid =~ /^fig\|(\d+\.\d+)\.([^\.]+)/)
	{
	    my $genome = $1;
	    my $type = $2;
	    if (open(TMP,">>$FIG_Config::organisms/$genome/Features/$type/propagation_locks"))
	    {
		print TMP "$fid\t1\n";
	    }
	    close(TMP);
	}
    }
}

=head3 propagation_unlock_fid

    $fig->propagation_unlock_fid($user,$fid);

Unsets a propagation_lock on annotations for $fid.

=cut

sub propagation_unlock_fid {
    my($self,$user,$fid) = @_;

    if (! $self->table_exists('fid_propagation_locks'))       { return 0 }
    if ((! $user) || ($fid !~ /^fig\|\d+\.\d+/))  { return 0 }
    if (! $self->is_propagation_locked_fid($fid))             { return 1 }

    $self->add_annotation($fid,$user,"propagation_unlocked assignments");
    my $rdbH = $self->db_handle;
    my $relational_db_response = $rdbH->SQL("SELECT fid FROM fid_propagation_locks WHERE fid = '$fid' ");
    $rdbH->SQL("DELETE FROM fid_propagation_locks WHERE ( fid = '$fid' )");
    if ($fid =~ /^fig\|(\d+\.\d+)\.([^\.]+)/)
    {
	my $genome = $1;
	my $type = $2;
	if (open(TMP,">>$FIG_Config::organisms/$genome/Features/$type/propagation_locks"))
	{
	    print TMP "$fid\t0\n";
	}
	close(TMP);
    }
}

##################
use SOAP::Lite;

sub get_all_assertions {
    my($pegs) = @_;

    my $response = SOAP::Lite
        -> uri('http://www.nmpdr.org/AnnoClearinghouse_SOAP')
        -> proxy('http://clearinghouse.nmpdr.org/aclh-soap.cgi')
        -> get_all_annotations( $pegs );

    if (! $response) { return () }
    my $result = $response->result;
    if (! $result)   { return () }

    my @assertions = ();
    foreach my $peg (@$pegs)
    {
	push @assertions, $result->{$peg};
    }
    return @assertions;
}

sub get_expert_assertions {
    my($pegs) = (@_ == 1) ? $_[0] : $_[1];

    my $response = SOAP::Lite
        -> uri('http://www.nmpdr.org/AnnoClearinghouse_SOAP')
        -> proxy('http://clearinghouse.nmpdr.org/aclh-soap.cgi')
        -> get_user_annotations( $pegs );

    if (! $response) { return () }
    my $result = $response->result;
    if (! $result)   { return () }
    my @assertions = ();
    foreach my $peg (keys(%$result))
    {
	my $x = $result->{$peg};
	push(@assertions,map { [$peg,@$_] } @$x);
    }

    return sort { &FIG::by_fig_id($a->[0],$b->[0]) } @assertions;
}
###############


=head3 delete_genomes

    $fig->delete_genomes(\@genomes);

Delete the specified genomes from the data store. This requires making
system calls to move and delete files.

=cut
#: Return Type ;
################################# make damn sure that you have enough disk ######################
### The following code represents a serious, major update.  Normally, one simply "marks" deleted
### genomes, which is quick and does not require halting the system.
sub delete_genomes {
    my($self,$genomes) = @_;

    die "Do not attempt to use delete_genomes";

    my $tmpD     = "$FIG_Config::temp/tmp.deleted.$$";
    my $tmp_Data = "$FIG_Config::temp/Data.$$";

    my %to_del = map { $_ => 1 } @$genomes;
    open(TMP,">$tmpD") || die "could not open $tmpD";

    my $genome;
    foreach $genome ($self->genomes) {
        if (! $to_del{$genome}) {
            print TMP "$genome\n";
        }
    }
    close(TMP);

    &run("extract_genomes $tmpD $FIG_Config::data $tmp_Data");
    print STDERR "Please bring the system down for a bit\n";
    system "echo \"System down due to update of genomes\n\" >> $tmp_Data/Global/why_down";
    &run("mv $FIG_Config::data $FIG_Config::data.deleted");
    &run("mv $tmp_Data $FIG_Config::data");
    &run("fig load_all");
    print STDERR "Now, you should think about deleting $FIG_Config::data.deleted\n";
    unlink("$FIG_Config::global/why_down");            ### start allowing CGIs to run
#   &run("rm -rf $FIG_Config::data.deleted");
}

### Mark a genome as deleted, but do not actually clean up anything.  That whole event
### requires "delete_genomes"
###
sub mark_deleted_genomes {
    my($self,$user,$genomes) = @_;
    my($genome);

    foreach $genome (@$genomes)
    {
	$self->broker_log("mark_genome_deleted", {
	    user => $user,
	    genome => $genome,
	});
        $self->log_update($user,$genome,$self->genus_species($genome),"Marked Deleted Genome $genome");
    }
    return $self->mark_deleted_genomes_body($user,$genomes);
}

sub mark_deleted_genomes_body {
    my($self,$user,$genomes) = @_;
    my($genome);

    my $rdbH = $self->db_handle;

    my $n = 0;
    foreach $genome (@$genomes)
    {
        if ($self->is_genome($genome) && open(DEL,">$FIG_Config::organisms/$genome/DELETED"))
        {
            print DEL "deleted\n";
            $rdbH->SQL("DELETE FROM genome WHERE ( genome = '$genome' )");
            $n++;
        }
        close(DEL);
    }
    $self->{_is_genome} = {};
    return $n;
}

sub unmark_deleted_genomes {
    my($self,$user,$genomes) = @_;
    my($genome);

    foreach $genome (@$genomes)
    {
        $self->log_update($user,$genome,$self->genus_species($genome),"Unmarked Deleted Genome $genome");
    }

    my $rdbH = $self->db_handle;

    my $n = 0;
    foreach $genome (@$genomes)
    {
        if (-s "$FIG_Config::organisms/$genome/DELETED")
        {
            unlink("$FIG_Config::organisms/$genome/DELETED");
            &run("compute_genome_counts $genome");
            $n++;
        }
    }
    $self->{_is_genome} = {};
    return $n;
}

sub log_corr {
    my($self,$user,$genome, $mapping,$msg) = @_;

    my $gs = $self->genus_species($genome);
    $self->log_update($user,$genome,$gs,"Logged correspondence for $genome [$msg]",$mapping);
}

=head3 replaces

	my $old_genome_id = $fig->replaces($new_genome_id);

If the new genome replaces any old ones, as denoted by the REPLACES file contents, that will be returned. Else undef is returned

=cut

sub replaces {
	my ($self, $genome) =@_;
	my $ret;
	if (-e "$FIG_Config::organisms/$genome/REPLACES") {
		$ret=`head -n 1 $FIG_Config::organisms/$genome/REPLACES`;
		chomp($ret);
	}
	return $ret;
}

sub replace_genome {
    my($self,$user,$old_genome,$genomeF, $mapping, $force, $skipnr) = @_;

    ($genomeF =~ /(\d+\.\d+)$/)
        || die "$genomeF must have a valid genome ID as the last part of the path";
    my $genome = $1;

    open(TMP,"<$genomeF/GENOME") || die "could not open $genome/GENOME";
    my $gs = <TMP>;
    chomp $gs;
    close(TMP);

    $self->log_update($user,$genome,$gs,"Replaced genome $old_genome with $genome\n$genomeF $force $skipnr",$genomeF,$mapping);

    $self->mark_deleted_genomes($user,[$old_genome]);
    return $self->add_genome_body($user,$genomeF,$force,$skipnr);
}

=head3 add_genome

    my $ok = $fig->add_genome($genomeF, $force, $skipnr);

Add a new genome to the data store. A genome's data is kept in a directory
by itself, underneath the main organism directory. This method essentially
moves genome data from an external directory to the main directory and
performs some indexing tasks to integrate it.

=over 4

=item genomeF

Name of the directory containing the genome files. This should be a
fully-qualified directory name. The last segment of the directory
name should be the genome ID.

=item force

This will ignore errors thrown by verify_genome_directory. This is bad, and you should
never do it, but I am in the situation where I need to move a genome from one machine
to another, and although I trust the genome I can't.

=item skipnr

We don't always want to add the proteins into the nr database. For example wih a metagnome that has been called by blastx. This will just skip appending the proteins into the NR file.

=item RETURN

Returns TRUE if successful, else FALSE.

=back

=cut
#: Return Type $;
sub add_genome {
    my($self,$user,$genomeF, $force, $skipnr, $dont_mark_complete) = @_;

    ($genomeF =~ /(\d+\.\d+)$/)
        || die "$genomeF must have a valid genome ID as the last part of the path";
    my $genome = $1;

    open(TMP,"<$genomeF/GENOME") || die "could not open $genome/GENOME";
    my $gs = <TMP>;
    chomp $gs;
    close(TMP);

    my $rc = $self->add_genome_body($user,$genomeF,$force,$skipnr,$dont_mark_complete);

    if ($rc)
    {
	$self->broker_log("add_genome", {
	    user => $user,
	    genome => $genome,
	    genome_file => $genomeF,
	    force => ($force ? 1 : 0),
	});
        $self->log_update($user,$genome,$gs,"Added genome $genome\n$genomeF $force $skipnr",$genomeF);
    }

    return $rc;
}

sub add_genome_body {
    my($self,$user,$genomeF, $force, $skipnr,$dont_mark_complete) = @_;

    my $rc = 0;

    my(undef, $path, $genome) = File::Spec->splitpath($genomeF);

    if ($genome !~ /^\d+\.\d+$/) {
        warn "Invalid genome filename $genomeF\n";
        return $rc;
    }

    if (-d $FIG_Config::organisms/$genome) {
        warn "Organism already exists for $genome\n";
        return $rc;
    }


    #
    # We're okay, it doesn't exist.
    #

    my @errors = `$FIG_Config::bin/verify_genome_directory $genomeF`;

    if (@errors) {
        print STDERR "Errors found while verifying genome directory $genomeF:\n";
        print STDERR join("", @errors);

	#
	# Special case check: If the only errors returned are peg_tbl_stop_missing, we're
	# probably hitting a possibly_truncated bug. Let the process continue.
	#

	my @corrupt = grep { /corrupt/ } @errors;
	if (@corrupt == 1 and $corrupt[0] =~ /is corrupt \(peg_tbl_stop_missing=(\d+)\)/)
	{
	    my $count = $1;
	    my $s = $count > 1 ? "s" : "";
	    print "Only error is $count peg_tbl_stop_missing error$s, continuing\n";
	}
        elsif (!$force)
        {
            return $rc;
        }
        else
        {
            warn "Skipped these errors and continued. You should not do this";
        }
    }

    my $sysrc = system("cp -r $genomeF $FIG_Config::organisms");
    if ($sysrc != 0)
    {
        warn "Failure copying $genomeF to $FIG_Config::organisms\n";
        return $rc;
    }

    my $genome_dir = "$FIG_Config::organisms/$genome";

    $sysrc = system("chmod -R 777 $genome_dir");
    if ($sysrc != 0)
    {
        warn "Command failed: chmod -R 777 $genome_dir\n";
        return $rc;
    }

    if (-s "$genome_dir/COMPLETE")
    {
        if ($dont_mark_complete)
        {
            print STDERR "$genome was marked as \"complete\", but moving to WAS_MARKED_COMPLETE\n";
            rename("$genome_dir/COMPLETE", "$genome_dir/WAS_MARKED_COMPLETE");
        }
        else
        {
            print STDERR "$genome was marked as \"complete\"\n";
        }
    }
    else
    {
        #
        # Not marked complete; assess completeness.
        #

        my $sysrc = system("$FIG_Config::bin/assess_completeness $genome_dir > $genome_dir/assess_completeness.out 2>&1");
        if ($sysrc != 0)
        {
            warn "assess_completeness $genome_dir failed; continuing with installation.\n";
        }
        else
        {
            if (-s "$genome_dir/PROBABLY_COMPLETE")
            {
                print STDERR "Assessed $genome to be probably complete\n";
                if ($dont_mark_complete)
                {
                    print STDERR "Not copying PROBABLY_COMPLETE to COMPLETE; this will need to be done later\n";
                }
                else
                {
                    my $cp = "cp -p $genome_dir/PROBABLY_COMPLETE $genome_dir/COMPLETE";
                    $sysrc = system($cp);
                    $sysrc == 0 or warn "Command failed, continuing: $cp\n";
                }
            }
            else
            {
                print STDERR "Assessed $genome to not be probably complete\n";
            }
        }
    }

    #
    # If this is an NMPDR organism and wasn't marked COMPLETE, mark it anyway so that it
    # get imported into the NMPDR. This will go away at some point.
    #

    my $nmpdr_group = &FIG::file_head("$genome_dir/NMPDR");
    chomp $nmpdr_group;
    if (! -s "$genome_dir/COMPLETE" and $nmpdr_group ne '')
    {
	open(P, ">$genome_dir/COMPLETE");
	print P "Marked complete due to NMPDR membership in $nmpdr_group\n";
	close(P);
    }

    #
    # If this was a RAST genome that has imp_annotations and imp_assigned_functions files,
    # rename any existing annotations/assigned_functions files to rast_XX and copy
    # imp_XX to XX.
    #
    #
    # HOWEVER, do not do this if there is a PSEED_RAST file since this
    # is a PSEED->SEED import.
    #

    if (-f "$genome_dir/RAST" && ! -f "$genome_dir/PSEED_RAST")
    {
	for my $base ('annotations', 'assigned_functions')
	{
	    my $imp = "$genome_dir/imp_$base";
	    my $file = "$genome_dir/$base";
	    my $rast = "$genome_dir/rast_$base";

	    if (-f $file)
	    {
		print "Rename $file to $rast\n";
		rename($file, $rast);
	    }
	    if (-f $imp)
	    {
		print "Copy $imp to $file\n";
		copy($imp, $file);
	    }
	}
    }

    $rc = $self->index_genome($genome, $user);

    return $rc;
}

sub index_genome
{
    my($self, $genome, $user) = @_;
    my $sysrc;
    my $rc;

    my $genome_dir = $self->organism_directory($genome);

    print "index_contigs $genome\n";
    $sysrc = system("index_contigs $genome");
    $sysrc == 0 or
        warn "index_contigs $genome failed; continuing with installation\n";

    print "compute_genome_counts $genome\n";
    $sysrc = system("compute_genome_counts $genome");
    $sysrc == 0 or
        warn "compute_genome_counts $genome failed; continuing with installation\n";

    print "load_features $genome\n";
    $sysrc = system("load_features $genome");
    $sysrc == 0 or
        warn "load_features $genome failed; continuing with installation\n";

    $rc = 1;
    if (-s "$genome_dir/Features/peg/fasta")
    {
	print "index_translations $genome\n";
        $sysrc = system("index_translations $genome");
        $sysrc == 0 or
            warn "index_translations $genome failed; continuing with installation\n";
    }

    if ((-s "$genome_dir/assigned_functions") ||
        (-d "$genome_dir/UserModels"))
    {
	print "add_assertions_of_function $genome\n";
        $sysrc = system("add_assertions_of_function $genome");
        $sysrc == 0 or warn "add_assertions_of_function $genome failed; continuing with installation\n";
    }

    if (-s "$genome_dir/annotations")
    {
	print "index_annotations $genome\n";
        $sysrc = system("index_annotations $genome");
        $sysrc == 0 or warn "index_annoations $genome failed; continuing with installation\n";
    }

    #
    # New support for installing precomputed data coming out of the RAST runs.
    #
    # PCHs are installed with install_new_coupling_data.
    #

    my $pchs = "$genome_dir/pchs";
    my $pch_scores = "$genome_dir/pchs.scored";

    if (-f $pchs and  -f $pch_scores)
    {
	print "install_new_coupling_data $genome $pchs $pch_scores\n";
	$sysrc = system("$FIG_Config::bin/install_new_coupling_data",
			$genome,
			$pchs,
			$pch_scores);
	if ($sysrc == 0)
	{
	    print "PCHs installed, indexing.\n";
	    $sysrc = system("$FIG_Config::bin/load_coupling", $genome);
	    if ($sysrc != 0)
	    {
		warn "load_coupling $genome failed with rc=$sysrc\n";
	    }
	}
	else
	{
	    warn "Error $sysrc installing coupling data";
	}
    }

    #
    # Make sure that the features are registered for this genome. We assume here that
    # the genome is already registered (as it should be if we came from RAST).
    #

    my $dh = new DirHandle("$genome_dir/Features");
    for my $ftype ($dh->read())
    {
	my $path = "$genome_dir/Features/$ftype";
	next if $ftype =~ /^\./ or ! -d $path;

	my $fh = new FileHandle("<$path/tbl");
	if (!$fh)
	{
	    warn "Cannot open tbl file in feature directory $path: $!";
	    next;
	}
	#
	# Find the largest feature in use.
	#
	my $max = -1;
	while (<$fh>)
	{
	    chomp;
	    my($fid) = split(/\t/);
	    if ($fid =~ /^fig\|\d+\.\d+\.[^.]+\.(\d+)/)
	    {
		$max = $1 > $max ? $1 : $max;
	    }
	}
	close($fh);

	print "Done registering features\n";
	#
	# See what the clearinghouse has, and register features if they are not there.
	#
	if ($ClearinghouseOK)
	{
	    my $clnext = $self->clearinghouse_next_feature_id($genome, $ftype);
	    if ($clnext <= $max)
	    {
		#
		# Not enough features are registered in the clearinghouse. ($clnext needs to be $max + 1)
		# Register some more.
		#
		
		my $missing = $max - $clnext + 1;
		my $start = $self->clearinghouse_register_features($genome, $ftype, $missing);
		if (defined($start))
		{
		    print "Registered $missing new features of type $ftype on $genome (start=$start)\n";
		}
	    }
	}
    }

    #
    # Walk the functions we have just assigned and apply any renames from the funcrole rename log.
    #
    my $rename = $self->read_role_rename_log();
    my @feats = $self->all_features($genome);
    my $funcs = $self->function_of_bulk(\@feats);
    for my $fid (@feats)
    {
	my $func = $funcs->{$fid};
	my $orig = $func;
	my $new;
        my $last_n;
	while (my $new_ent = $rename->{$func})
	{
	    my($new, $n) = @$new_ent;
	    if ($n < $last_n)
	    {
		warn "Breaking off renames for $fid $orig due to $n < $last_n\n";
		last;
	    }
	    $last_n = $n;
	    $func = $new;
	}
	if ($func ne $orig)
	{
	    print "Rename $fid: $orig => $func\n";
	    $self->add_annotation($fid, $user, "Changing assignment from $orig to $func based on role rename log");
	    $self->assign_function($fid, $user, $func);
	}
    }

    return $rc;
}

=head3 assess_completeness

    my $lengths = [[contigid1 => $length1], [contigid2 => $lengt2], ...];
    my($complete, $fraction_in_large_contigs, $total_len) = &FIG::assess_completeness($lengths);

=cut

sub assess_completeness
{
    my($contig_lengths) = @_;

    my $minfrac =    0.7;
    my $minlen  =  20000;
    my $minsize = 300000;

    my $ttlen = 0;
    my $inbig = 0;
    foreach my $ent (@$contig_lengths)
    {
	my($id, $len) = @$ent;
	$ttlen += $len;
	if ($len >= $minlen)
	{
	    $inbig += $len;
	}
    }

    my $frac = 100 * $inbig / $ttlen;

    my $complete = (($ttlen >= $minsize) && ($inbig >= $minfrac * $ttlen) ) ? 1 : 0;
    return ($complete, $frac, $ttlen);
}

sub read_role_rename_log
{
    my($self) = @_;
    my $logfile = "$FIG_Config::data/Logs/functionalroles.rewrite";

    my $log = {};

    my $lf = new FileHandle("<$logfile");
    if (!$lf)
    {
	warn "Cannot read role rename log $logfile: $!";
	return $log;
    }

    my $n = 0;

    while (<$lf>)
    {
        if (/^Role\s+(.*)\s+was\s+replaced\s+by\s+(.*)/)
        {
	    my $from = $1;
	    my $to = $2;
            $n++;
            if ($from eq $to)
            {
                warn "CYCLE $from->$to in rename table\n";
                next;
            }
	    elsif ($to eq '')
	    {
		warn "EMPTY rename $from->EMPTY in rename table\n";
		next;
	    }
            $log->{$1} = [$2, $n];
        }
    }
    close($lf);
    return $log;
}

sub get_index {
    my($self,$gs) = @_;

    my($index,$max);
    $gs || confess "MISSING GS";

    my $indexF = "$FIG_Config::data/Logs/GenomeLog/index";
    if (open(INDEX,"<$indexF"))
    {
        while ((! $index) && ($_ = <INDEX>))
        {
            if ($_ =~ /^(\d+)/)
            {
                $max = $1;
                if (($_ =~ /^(\d+)\t(\S.*\S)/) && ($2 eq $gs))
                {
                    $index = $1;
                }
            }
        }
        close(INDEX);
    }

    if (! $index)
    {
        open(INDEX,">>$indexF") || die "could not open $indexF";
        $index = defined($max) ? $max+1 : 1;
        print INDEX "$index\t$gs\n";
        close(INDEX);
        &verify_dir("$FIG_Config::data/Logs/GenomeLog/Entries/$index");
    }
    return $index;
}

sub log_update {
    my($self,$user,$genome,$gs,$msg,@data) = @_;

    my $time_made = time;
    &verify_dir("$FIG_Config::data/Logs/GenomeLog");
    my $index_id = $self->get_index($gs);
    $index_id || die "could not make an index entry for $gs";
    my $gs_dir = "$FIG_Config::data/Logs/GenomeLog/Entries/$index_id";
    &verify_dir($gs_dir);
    my($i,$file_or_dir,@tars);
    for ($i=0; ($i < @data); $i++)
    {
        $file_or_dir = $data[$i];
        my($dir,$file);
        if ($file_or_dir =~ /^(.*)\/([^\/]+)$/)
        {
            ($dir,$file) = ($1,$2);
        }
        else
        {
            ($dir,$file) = (".",$file_or_dir);
        }
        my $tar = "$gs_dir/$time_made.$i.tgz";
        &run("cd $dir; tar czf $tar $file");
        push(@tars,$tar);
    }
    open(LOG,">>$gs_dir/log")
        || die "could not open $gs_dir/log";
    print LOG "$time_made\n$user\n$genome\n$msg\n";
    if (@tars > 0)
    {
        print LOG join(",",@tars),"\n";
    }
    print LOG "//\n";
    close(LOG);
}

=head3 parse_genome_args

    my ($mode, @genomes) = FIG::parse_genome_args(@args);

Extract a list of genome IDs from an argument list. If the argument list is empty,
return all the genomes in the data store.

This is a function that is performed by many of the FIG command-line utilities. The
user has the option of specifying a list of specific genome IDs or specifying none
in order to get all of them. If your command requires additional arguments in the
command line, you can still use this method if you shift them out of the argument list
before calling. The $mode return value will be C<all> if the user asked for all of
the genomes or C<some> if he specified a list of IDs. This is useful to know if,
for example, we are loading a table. If we're loading everything, we can delete the
entire table; if we're only loading some genomes, we must delete them individually.

This method uses the genome directory rather than the database because it may be used
before the database is ready.

=over 4

=item args1, args2, ... argsN

List of genome IDs. If all genome IDs are to be processed, then this list should be
empty.

=item RETURN

Returns a list. The first element of the list is C<all> if the user is asking for all
the genome IDs and C<some> otherwise. The remaining elements of the list are the
desired genome IDs.

=back

=cut

sub parse_genome_args {
    # Get the parameters.
    my @args = @_;
    # Check the mode.
    my $mode = (@args > 0 ? 'some' : 'all');
    # Build the return list.
    my @retVal = ($mode);
    # Process according to the mode.
    if ($mode eq 'all') {
        # We want all the genomes, so we get them from the organism directory.
        my $orgdir = "$FIG_Config::organisms";
        opendir( GENOMES, $orgdir ) || Confess("Could not open directory $orgdir");
        push @retVal, grep { $_ =~ /^\d/ } readdir( GENOMES );
        closedir( GENOMES );
    } else {
        # We want only the genomes specified by the user.
        push @retVal, @args;
    }
    # Return the result.
    return @retVal;
}

=head3 reload_table

    $fig->reload_table($mode, $table, $flds, $xflds, $fileName, $keyList, $keyName);

Reload a database table from a sequential file. If I<$mode> is C<all>, the table
will be dropped and re-created. If I<$mode> is C<some>, the data for the individual
items in I<$keyList> will be deleted before the table is loaded. Thus, the load
process is optimized for the type of reload.

=over 4

=item mode

C<all> if we are reloading the entire table, C<some> if we are only reloading
specific entries.

=item table

Name of the table to reload.

=item flds

String defining the table columns, in SQL format. In general, this is a
comma-delimited set of field specifiers, each specifier consisting of the
field name followed by the field type and any optional qualifiers (such as
C<NOT NULL> or C<DEFAULT>); however, it can be anything that would appear
between the parentheses in a C<CREATE TABLE> statement. The order in which
the fields are specified is important, since it is presumed that is the
order in which they are appearing in the load file.

=item xflds

Reference to a hash that describes the indexes. The hash is keyed by index name.
The value is the index's field list. This is a comma-delimited list of field names
in order from most significant to least significant. If a field is to be indexed
in descending order, its name should be followed by the qualifier C<DESC>. For
example, the following I<$xflds> value will create two indexes, one for name followed
by creation date in reverse chronological order, and one for ID.

    { name_index => "name, createDate DESC", id_index => "id" }

=item fileName

Fully-qualified name of the file containing the data to load. Each line of the
file must correspond to a record, and the fields must be arranged in order and
tab-delimited. If the file name is omitted, the table is dropped and re-created
but not loaded.

=item keyList

Reference to a list of the IDs for the objects being reloaded. This parameter is
only used if I<$mode> is C<some>.

=item keyName (optional)

Name of the key field containing the IDs in the keylist. If omitted, C<genome> is
assumed.

=back

=cut

sub reload_table {
    # Get the parameters.
    my ($self, $mode, $table, $flds, $xflds, $fileName, $keyList, $keyName) = @_;
    if (!defined $keyName) {
        $keyName = 'genome';
    }
    # Get the database handler.
    my $dbf = $self->{_dbf};
    # Call the DBKernel method.
    $dbf->reload_table($mode, $table, $flds, $xflds, $fileName, $keyList, $keyName);
}

=head3 enqueue_similarities

    FIG::enqueue_similarities(\@fids);

Queue the passed Feature IDs for similarity computation. The actual
computation is performed by L</create_sim_askfor_pool>. The queue is a
persistent text file in the global data directory, and this method
essentially writes new IDs on the end of it.

=over 4

=item fids

Reference to a list of feature IDs.

=back

=cut
#: Return Type ;
sub enqueue_similarities {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($fids) = @_;
    my $fid;

    my $sim_q = "$FIG_Config::global/queued_similarities";

    open(TMP,">>$sim_q")
        || die "could not open $sim_q";

    #
    # We need to lock here so that if a computation is creating a snapshot of the
    # queue, we block until it's done.
    #

    flock(TMP, LOCK_EX) or die "Cannot lock $sim_q\n";
    seek(TMP, 0, 2);

    foreach $fid (@$fids) {
        print TMP "$fid\n";
    }
    close(TMP);
}

=head3 export_similarity_request

Creates a similarity computation request from the queued similarities and
the current NR.

We keep track of the exported requests in case one gets lost.

=cut

sub export_similarity_request {
    my($self, $user_req_dir) = @_;

    my $nr_file = "$user_req_dir/nr";
    my $fasta_file = "$user_req_dir/fasta";
    my $peg_syn_file = "$user_req_dir/peg.synonyms";

    my $req_dir = "$FIG_Config::fig/var/sim_requests";
    &verify_dir("$FIG_Config::fig/var");
    &verify_dir($req_dir);

    $req_dir = "$req_dir/" . time;
    &verify_dir($req_dir);

    #
    # Open all of our output files before zeroing out the sim queue, in case
    # there is a problem.
    #

    open(my $user_fasta_fh, ">$fasta_file") or confess "Cannot open $fasta_file for writing: $!";
    open(my $fasta_fh, ">$req_dir/fasta.in");

    open(my $user_nr_fh, ">$nr_file") or confess "Cannot open $nr_file for writing: $!";
    open(my $nr_fh, ">$req_dir/nr") or confess "Cannot open $req_dir/nr for writing: $!";

    open(my $user_peg_syn_fh, ">$peg_syn_file") or confess "Cannot open $peg_syn_file for writing: $!";
    open(my $peg_syn_fh, ">$req_dir/peg.synonyms") or confess "Cannot open $req_dir/peg.synonyms for writing: $!";

    open(my $nr_read_fh, "<$FIG_Config::data/Global/nr") or die "Cannot open $FIG_Config::data/Global/nr for reading: $!";
    open(my $peg_syn_read_fh, "<$FIG_Config::data/Global/peg.synonyms") or die "Cannot open $FIG_Config::data/Global/peg.synonyms for reading: $!";

    my $sim_q = "$FIG_Config::global/queued_similarities";

    #
    # We need to lock here so that if a computation is creating a snapshot of the
    # queue, we block until it's done.
    #

    open(my $sim_q_lock, ">>$sim_q") or confess "could not open $sim_q";
    flock($sim_q_lock, LOCK_EX) or confess "Cannot lock $sim_q\n";

    #
    # Everything open & locked, start copying.
    #

    copy("$sim_q", "$req_dir/q") or confess "Copy $sim_q $req_dir/q failed: $!";
    copy("$sim_q", "$user_req_dir/q") or confess "Copy $sim_q $user_req_dir/q failed: $!";

    #
    # Copy the contents of the sim queue to the "expected import" queue;
    # this is a list of pegs for which we expect sims to be computed and installed
    # at some point.
    #
    # We also lock on the pending queue file.
    #

    if (not(open(SQ, "<$sim_q")))
    {
        warn "Cannot open $sim_q for reading: $!\n";
    }
    else
    {
        if (open(AW, ">>$FIG_Config::global/pending_similarities"))
        {
            flock(AW, LOCK_EX);
            seek(AW, 0, 2);

            while (<SQ>)
            {
                print AW @_;
            }
            close(AW);
        }
        else
        {
            warn "Could not open $FIG_Config::global/pending_similarities: $!\n";
        }
        close(SQ);
    }

    my($buf);
    while (1) {
        my $n = read($nr_read_fh, $buf, 4096);
        defined($n) or confess "Error reading nr: $!";
        last unless $n;
        syswrite($user_nr_fh, $buf) or confess "Error writing $nr_file: $!";
        syswrite($nr_fh, $buf) or confess "Error writing $req_dir/nr: $!";
    }

    close($nr_read_fh);
    close($nr_fh);
    close($user_nr_fh);

    while (1) {
        my $n = read($peg_syn_read_fh, $buf, 4096);
        defined($n) or confess "Error reading peg.synonyms: $!";
        last unless $n;
        syswrite($user_peg_syn_fh, $buf) or confess "Error writing $peg_syn_file: $!";
        syswrite($peg_syn_fh, $buf) or confess "Error writing $req_dir/peg.synonyms: $!";
    }

    close($peg_syn_read_fh);
    close($peg_syn_fh);
    close($user_peg_syn_fh);

    #
    # We can zero out the queue and unlock now.
    #

    open(F, ">$sim_q") or die "Cannot open $sim_q to truncate it: $!\n";
    close(F);

    close($sim_q_lock);

    #
    # Generate the fasta input from the queued ids.
    #

    open(my $q_fh, "<$req_dir/q");
    while (my $id = <$q_fh>) {
        chomp $id;

        my $seq = $self->get_translation($id);

        display_id_and_seq($id, \$seq, $user_fasta_fh);
        display_id_and_seq($id, \$seq, $fasta_fh);
    }
    close($q_fh);

    close($user_fasta_fh);
    close($fasta_fh);
}

=head3 create_sim_askfor_pool

    $fig->create_sim_askfor_pool($chunk_size);

Creates an askfor pool, which a snapshot of the current NR and similarity
queue. This process clears the old queue.

The askfor pool needs to keep track of which sequences need to be
calculated, which have been handed out, etc. To simplify this task we
chunk the sequences into fairly small numbers (20k characters) and
allocate work on a per-chunk basis. We make use of the relational
database to keep track of chunk status as well as the seek locations
into the file of sequence data. The initial creation of the pool
involves indexing the sequence data with seek offsets and lengths and
populating the sim_askfor_index table with this information and with
initial status information.

=over 4

=item chunk_size

Number of features to put into a processing chunk. The default is 15.

=back

=cut
#: Return Type $;
sub create_sim_askfor_pool {
    my($self, $chunk_size) = @_;

    $chunk_size = 20000 unless $chunk_size =~ /^\d+$/;

    my $pool_dir = "$FIG_Config::fig/var/sim_pools";
    &verify_dir($pool_dir);

    #
    # Lock the pool directory.
    #
    open(my $lock, ">$pool_dir/lockfile");

    flock($lock, LOCK_EX);

    my $num = 0;
    if (open(my $toc, "<$pool_dir/TOC")) {
        while (<$toc>) {
            chomp;
            # print STDERR "Have toc entry  $_\n";
            my ($idx, $time, $str) = split(/\s+/, $_, 3);

            $num = max($num, $idx);
        }
        close($toc);
    }
    $num++;
    open(my $toc, ">>$pool_dir/TOC") or die "Cannot write $pool_dir/TOC: $!\n";

    print $toc "$num ", time(), " New toc entry\n";
    close($toc);

    my $cpool_id = sprintf "%04d", $num;
    my $cpool_dir = "$pool_dir/$cpool_id";

    #
    # All set, create the directory for this pool.
    #

    &verify_dir($cpool_dir);

    #
    # Now we can copy the nr and sim queue here.
    # Do this stuff inside an eval so we can clean up
    # the lockfile.
    #

    eval {
        my $sim_q = "$FIG_Config::global/queued_similarities";

        copy("$sim_q", "$cpool_dir/q");
        copy("$FIG_Config::data/Global/nr", "$cpool_dir/nr");

        open(F, ">$sim_q") or die "Cannot open $sim_q to truncate it: $!\n";
        close(F);
    };

    unlink("$pool_dir/lockfile");
    close($lock);

    #
    # We've created our pool; we can now run the formatdb and
    # extract the sequences for the blast run.
    #
    my $child_pid = $self->run_in_background(
        sub {
            #
            # Need to close db or there's all sorts of trouble.
            #

            my $cmd = "$FIG_Config::ext_bin/formatdb -i $cpool_dir/nr -p T -l $cpool_dir/formatdb.log";
            print "Will run '$cmd'\n";
            &run($cmd);
            print "finished. Logfile:\n";
            print &FIG::file_read("$cpool_dir/formatdb.log");
            unlink("$cpool_dir/formatdb.pid");
        });
    warn "Running formatdb in background job $child_pid\n";
    open(FPID, ">$cpool_dir/formatdb.pid");
    print FPID "$child_pid\n";
    close(FPID);

    my $db = $self->db_handle();
    if (!$db->table_exists("sim_queue")) {
        $db->create_table(tbl => "sim_queue",
                  flds => "qid varchar(32), chunk_id INTEGER, seek INTEGER, len INTEGER, " .
                  "assigned BOOL, finished BOOL, output_file varchar(255), " .
                  "worker_pid INTEGER, start_time timestamp, " .
                  "assignment_expires INTEGER, worker_info varchar(255)"
                 );
    }

    #
    # Write the fasta input file. Keep track of how many have been written,
    # and write seek info into the database as appropriate.
    #

    open(my $seq_fh, ">$cpool_dir/fasta.in");

    my($chunk_idx, $chunk_begin, $seq_idx);

    my $cur_size = 0;

    $chunk_idx = 0;
    $chunk_begin = 0;
    $seq_idx = 0;

    my $tmpfile = "$FIG_Config::temp/simseek.$$";
    open(my $tmpfh, ">$tmpfile") or confess "Cannot open tmpfile $tmpfile: $!";

    open(my $q_fh, "<$cpool_dir/q");
    while (my $id = <$q_fh>) {
        chomp $id;

        my $seq = $self->get_translation($id);

        #
        # check if we're at the beginning of a chunk
        #

        print $seq_fh ">$id\n$seq\n";

        #
        # Check if we're at the end of a chunk
        #

        $cur_size += length($seq);
        if ($cur_size >= $chunk_size) {
            my $chunk_end = tell($seq_fh);
            my $chunk_len = $chunk_end - $chunk_begin;

            print $tmpfh join("\t", $cpool_id, $chunk_idx, $chunk_begin, $chunk_len, 'FALSE', 'FALSE',
                              '\N', '\N', '\N', '\N', '\N'), "\n";
            $chunk_idx++;
            $chunk_begin = $chunk_end;
            $cur_size = 0;
        }
        $seq_idx++;
    }

    if ($cur_size > 0) {
        my $chunk_end = tell($seq_fh);
        my $chunk_len = $chunk_end - $chunk_begin;

        print $tmpfh join("\t", $cpool_id, $chunk_idx, $chunk_begin, $chunk_len, 'FALSE', 'FALSE',
                          '\N', '\N', '\N', '\N', '\N'), "\n";
    }

    close($q_fh);
    close($seq_fh);
    close($tmpfh);

    warn "Write seqs from $tmpfile\n";

    $self->db_handle->load_table(tbl => 'sim_queue',
                                 file => $tmpfile);

#    unlink($tmpfile);

#     for my $seek (@seeks)
#     {
#       my($cpool_id, $chunk_idx, $chunk_begin, $chunk_len) = @$seek;

#       $db->SQL("insert into sim_queue (qid, chunk_id, seek, len, assigned, finished) " .
#                "values('$cpool_id', $chunk_idx, $chunk_begin, $chunk_len, FALSE, FALSE)");
#     }

    return $cpool_id;
}

#=head3 get_sim_queue
#
#usage: get_sim_queue($pool_id, $all_sims)
#
#Returns the sims in the given pool. If $all_sims is true, return the entire queue. Otherwise,
#just return the sims awaiting processing.
#
#=cut

sub get_sim_queue  {
    my($self, $pool_id, $all_sims) = @_;
}

=head3 get_sim_work

    my ($nrPath, $fasta) = $fig->get_sim_work();

Get the next piece of sim computation work to be performed. Returned are
the path to the NR and a string containing the fasta data.

=cut

sub get_sim_work {

    my ($self) = @_;

    #
    # For now, just don't care about order of data that we get back.
    #

    my $db = $self->db_handle();
    my $lock = FIG::SimLock->new;

    my $work = $db->SQL(qq(SELECT qid, chunk_id, seek, len
                           FROM sim_queue
                           WHERE not finished AND not assigned
                           LIMIT 1));
    print "Got work ", Dumper($work), "\n";

    if (not $work or @$work == 0) {
        return undef;
    }

    my($cpool_id, $chunk_id, $seek, $len) = @{$work->[0]};

    my $pool_dir = "$FIG_Config::fig/var/sim_pools";
    my $cpool_dir = "$pool_dir/$cpool_id";

    my $nr = "$cpool_dir/nr";
    open(my $fh, "<$cpool_dir/fasta.in");
    seek($fh, $seek, 0);
    my $fasta;
    read($fh, $fasta, $len);

    $db->SQL(qq(UPDATE sim_queue
                SET assigned = true
                WHERE qid = ? AND chunk_id = ?), undef,
             $cpool_id, $chunk_id);

    return($cpool_id, $chunk_id, $nr, $fasta, "$cpool_dir/out.$chunk_id");
}

sub sim_work_working
{
    my($self, $pool, $chunk, $host, $pid) = @_;

    my $db = $self->db_handle();
    my $lock = FIG::SimLock->new;

    my $res = $db->SQL(qq(UPDATE sim_queue
                          SET worker_pid = ?, start_time = NOW(), worker_info = ?
                          WHERE qid = ? AND chunk_id = ?),
                       undef,
                       $pid, $host, $pool, $chunk);
}

=head3 sim_work_done

    $fig->sim_work_done($pool_id, $chunk_id, $out_file);

Declare that the work in pool_id/chunk_id has been completed, and output written
to the pool directory (get_sim_work gave it the path).

=over 4

=item pool_id

The ID number of the pool containing the work that just completed.

=item chunk_id

The ID number of the chunk completed.

=item out_file

The file into which the work was placed.

=back

=cut

sub sim_work_done {
    my ($self, $pool_id, $chunk_id, $out_file) = @_;

    if (! -f $out_file) {
        Confess("sim_work_done: output file $out_file does not exist");
    }

    my $db = $self->db_handle();
    my $lock = FIG::SimLock->new;

    my $dbh = $db->{_dbh};

    my $rows = $dbh->do(qq(UPDATE sim_queue
                           SET finished = TRUE, output_file = ?
                           WHERE qid = ? and chunk_id = ?), undef, $out_file, $pool_id, $chunk_id);
    if ($rows != 1) {
        if ($dbh->errstr) {
            Confess("Update not able to set finished=TRUE: ", $dbh->errstr);
        } else {
            Confess("Update not able to set finished=TRUE");
        }
    }
    #
    # Determine if this was the last piece of work for this pool. If so, we can
    # schedule the postprocessing work.
    #
    # Note we're still holding the lock.
    #

    my $out = $db->SQL(qq(SELECT chunk_id
                          FROM sim_queue
                          WHERE qid = ? AND not finished), undef, $pool_id);
    if (@$out == 0) {
        #
        # Pool is done.
        #
        $self->schedule_sim_pool_postprocessing($pool_id);
    }
}

=head3 schedule_sim_pool_postprocessing

    $fig->schedule_sim_pool_postprocessing($pool_id);

Schedule a job to do the similarity postprocessing for the specified pool.

=over 4

=item pool_id

ID of the pool whose similarity postprocessing needs to be scheduled.

=back

=cut

sub schedule_sim_pool_postprocessing {

    my($self, $pool_id) = @_;

    my $pool_dir = "$FIG_Config::fig/var/sim_pools";
    my $cpool_dir = "$pool_dir/$pool_id";

    my $js = JobScheduler->new();
    my $job = $js->job_create();

    my $spath = $job->get_script_path();
    open(my $sfh, ">$spath");
    print $sfh <<END;
    #!/bin/sh
    . $FIG_Config::fig_disk/config/fig-user-env.sh
    $FIG_Config::bin/postprocess_computed_sims $pool_id
END

    close($sfh);
    chmod(0775, $spath);

    #
    # Write the job ID to the subsystem queue dir.
    #

    open(J, ">$cpool_dir/postprocess_jobid");
    print J $job->get_id(), "\n";
    close(J);

    $job->enqueue();
}

=head3 postprocess_computed_sims

    $fig->postprocess_computed_sims($pool_id);

Set up to reduce, reformat, and split the similarities in a given pool. We build
a pipe to this pipeline:

    reduce_sims peg.synonyms 300 | reformat_sims nr | split_sims dest prefix

Then we put the new sims in the pool directory, and then copy to NewSims.

=over 4

=item pool_id

ID of the pool whose similarities are to be post-processed.

=back

=cut

sub postprocess_computed_sims {
    my($self, $pool_id) = @_;

    #
    # We don't lock here because the job is already done, and we
    # shouldn't (ha, ha) ever postprocess twice.
    #

    my $pool_dir = "$FIG_Config::fig/var/sim_pools";
    my $cpool_dir = "$pool_dir/$pool_id";

    my $sim_dir = "$cpool_dir/NewSims";
    &verify_dir($sim_dir);

    #
    # Open the processing pipeline.
    #

    my $reduce = "$FIG_Config::bin/reduce_sims $FIG_Config::global/peg.synonyms 300";
    my $reformat = "$FIG_Config::bin/reformat_sims $cpool_dir/nr";
    my $split = "$FIG_Config::bin/split_sims $sim_dir sims.$pool_id";
    open(my $process, "| $reduce | $reformat | $split");

    #
    # Iterate over all the sims files, taken from the database.
    #

    my $dbh = $self->db_handle()->{_dbh};
    my $files = $dbh->selectcol_arrayref(qq(SELECT output_file
                                            FROM sim_queue
                                            WHERE qid = ? and output_file IS NOT NULL
                                            ORDER BY chunk_id), undef, $pool_id);
    for my $file (@$files) {
        my $buf;
        open(my $fh, "<$file") or confess "Cannot sim input file $file: $!";
        while (read($fh, $buf, 4096)) {
            print $process $buf;
        }
        close($fh);
    }
    my $res = close($process);
    if (!$res) {
        if ($!) {
            confess "Error closing process pipeline: $!";
        } else {
            confess "Process pipeline exited with status $?";
        }
    }

    #
    # If we got here, it worked.  Copy the new sims files over to NewSims.
    #

    opendir(my $simdh, $sim_dir) or confess "Cannot open $sim_dir: $!";
    my @new_sims = grep { $_ !~ /^\./ } readdir($simdh);
    closedir($simdh);

    &verify_dir("$FIG_Config::data/NewSims");

    for my $sim_file (@new_sims) {
        my $target = "$FIG_Config::data/NewSims/$sim_file";
        if (-s $target) {
            Confess("$target already exists");
        }
        print "copying sim file $sim_file\n";
        &FIG::run("cp $sim_dir/$sim_file $target");
        &FIG::run("$FIG_Config::bin/index_sims $target");
    }
}

=head3 get_active_sim_pools

    @pools = $fig->get_active_sim_pools();

Return a list of the pool IDs for the sim processing queues that have
entries awaiting computation.

=cut
#: Return Type @;
sub get_active_sim_pools {
    my($self) = @_;

    my $dbh = $self->db_handle();

    my $res = $dbh->SQL("select distinct qid from sim_queue where not finished");
    return undef unless $res;

    return map { $_->[0] } @$res;
}

=head3 compute_clusters

    my @clusterList = $fig->compute_clusters(\@pegList, $subsystem, $distance);

Partition a list of PEGs into sections that are clustered close together on
the genome. The basic algorithm used builds a graph connecting PEGs to
other PEGs close by them on the genome. Each connected subsection of the graph
is then separated into a cluster. Singleton clusters are thrown away, and
the remaining ones are sorted by length. All PEGs in the incoming list
should belong to the same genome, but this is not a requirement. PEGs on
different genomes will simply find themselves in different clusters.

=over 4

=item pegList

Reference to a list of PEG IDs.

=item subsystem

Subsystem object for the relevant subsystem. This parameter is not used, but is
required for compatability with Sprout.

=item distance (optional)

The maximum distance between PEGs that makes them considered close. If omitted,
the distance is 5000 bases.

=item RETURN

Returns a list of lists. Each sub-list is a cluster of PEGs.

=back

=cut

sub compute_clusters {
    # Get the parameters.
    my ($self, $pegList, $subsystem, $distance, $location_hash) = @_;
    if (! defined $distance) {
        $distance = 5000;
    }

    my($peg,%by_contig);

    my $locsH;
    if (ref($location_hash))
    {
	$locsH = $location_hash;
    }
    else
    {
	my @locs = $self->feature_location_bulk($pegList);
	$locsH = {};
	$locsH->{$_->[0]} = $_->[1] for @locs;
    }

    foreach $peg (@$pegList)
    {
        my $loc;
        if ($loc = $locsH->{$peg})
        {
            my ($contig,$beg,$end) = $self->boundaries_of($loc);
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

=head3 get_sim_pool_info

    my ($total_entries, $n_finished, $n_assigned, $n_unassigned) = $fig->get_sim_pool_info($pool_id);

Return information about the given sim pool.

=over 4

=item pool_id

Pool ID of the similarity processing queue whose information is desired.

=item RETURN

Returns a four-element list. The first is the number of features in the
queue; the second is the number of features that have been processed; the
third is the number of features that have been assigned to a
processor, and the fourth is the number of features left over.

=back

=cut
#: Return Type @;
sub get_sim_pool_info {

    my($self, $pool_id) = @_;
    my($dbh, $res, $total_entries, $n_finished, $n_assigned, $n_unassigned);

    $dbh = $self->db_handle();

    $res = $dbh->SQL("select count(chunk_id) from sim_queue where qid = '$pool_id'");
    $total_entries = $res->[0]->[0];

    $res = $dbh->SQL("select count(chunk_id) from sim_queue where qid = '$pool_id' and finished");
    $n_finished = $res->[0]->[0];

    $res = $dbh->SQL("select count(chunk_id) from sim_queue where qid = '$pool_id' and assigned and not finished");
    $n_assigned = $res->[0]->[0];

    $res = $dbh->SQL("select count(chunk_id) from sim_queue where qid = '$pool_id' and not finished and not assigned");
    $n_unassigned = $res->[0]->[0];

    return ($total_entries, $n_finished, $n_assigned, $n_unassigned);
}

#=head3 get_sim_chunk
#
#usage: get_sim_chunk($n_seqs, $worker_id)
#
#Returns a chunk of $n_seqs of work.
#
#From Ross, about how sims are processed:
#
#Here is how I process them:
#
#
#    bash$ cd /Volumes/seed/olson/Sims/June22.out
#    bash$ for i in really*
#    > do
#    > cat  < $i >> /Volumes/laptop/new.sims
#    > done
#
#
#Then, I need to "reformat" them by adding to columns to each one
# and split the result into files of about 3M each This I do using
#
#reduce_sims /Volumes/laptop/NR/NewNR/peg.synonyms.june21  300 < /Volumes/laptop/new.sims |
#    reformat_sims /Volumes/laptop/NR/NewNR/checked.nr.june21   > /Volumes/laptop/reformated.sims
#rm /Volumes/laptop/new.sims
#split_sims /Volumes/laptop/NewSims sims.june24  reformated.sims
#rm reformatted.sims
#
#=cut

sub get_sim_chunk  {
    my($self, $n_seqs, $worker_id) = @_;
}

=head3 get_local_hostname

    my $result = FIG::get_local_hostname();

Return the local host name for the current processor. The name may be
stored in a configuration file, or we may have to get it from the
operating system.

=cut
#: Return Type $;
sub get_local_hostname {

    #
    # See if there is a FIGdisk/config/hostname file. If there
    # is, force the hostname to be that.
    #

    my $hostfile = "$FIG_Config::fig_disk/config/hostname";
    if (-f $hostfile) {
        my $fh;
        if (open($fh, $hostfile)) {
            my $hostname = <$fh>;
            chomp($hostname);
            return $hostname;
        }
    }

    #
    # First check to see if we our hostname is correct.
    #
    # Map it to an IP address, and try to bind to that ip.
    #

    local $/ = "\n";

    my $tcp = getprotobyname('tcp');

    my $hostname = `hostname`;
    chomp $hostname;

    my @hostent = gethostbyname($hostname);

    if (@hostent > 0) {
        my $sock;
        my $ip = $hostent[4];

        socket($sock, PF_INET, SOCK_STREAM, $tcp);
        if (bind($sock, sockaddr_in(0, $ip))) {
            #
            # It worked. Reverse-map back to a hopefully fqdn.
            #

            my @rev = gethostbyaddr($ip, AF_INET);
            if (@rev > 0) {
                my $host = $rev[0];
                #
                # Check to see if we have a FQDN.
                #

                if ($host =~ /\./) {
                    #
                    # Good.
                    #
                    return $host;
                } else {
                    #
                    # We didn't get a fqdn; bail and return the IP address.
                    #
                    return get_hostname_by_adapter()
                }
            } else {
                return inet_ntoa($ip);
            }
         } else {
            #
            # Our hostname must be wrong; we can't bind to the IP
            # address it maps to.
            # Return the name associated with the adapter.
            #
            return get_hostname_by_adapter()
        }
    } else {
        #
        # Our hostname isn't known to DNS. This isn't good.
        # Return the name associated with the adapter.
        #
        return get_hostname_by_adapter()
    }
}

=head3 get_hostname_by_adapter

    my $name = FIG::get_hostname_by_adapter();

Return the local host name for the current network environment.

=cut
#: Return Type $;
sub get_hostname_by_adapter {
    #
    # Attempt to determine our local hostname based on the
    # network environment.
    #
    # This implementation reads the routing table for the default route.
    # We then look at the interface config for the interface that holds the default.
    #
    #
    # Linux routing table:
    # [olson@yips 0.0.0]$ netstat -rn
    #     Kernel IP routing table
    #     Destination     Gateway         Genmask         Flags   MSS Window  irtt Iface
    #     140.221.34.32   0.0.0.0         255.255.255.224 U         0 0          0 eth0
    #     169.254.0.0     0.0.0.0         255.255.0.0     U         0 0          0 eth0
    #     127.0.0.0       0.0.0.0         255.0.0.0       U         0 0          0 lo
    #     0.0.0.0         140.221.34.61   0.0.0.0         UG        0 0          0 eth0
    #
    #     Mac routing table:
    #
    #     bash-2.05a$ netstat -rn
    #     Routing tables
    #
    #  Internet:
    #     Destination        Gateway            Flags    Refs      Use  Netif Expire
    #     default            140.221.11.253     UGSc       12      120    en0
    #     127.0.0.1          127.0.0.1          UH         16  8415486    lo0
    #     140.221.8/22       link#4             UCS        12        0    en0
    #     140.221.8.78       0:6:5b:f:51:c4     UHLW        0      183    en0    408
    #     140.221.8.191      0:3:93:84:ab:e8    UHLW        0       92    en0    622
    #     140.221.8.198      0:e0:98:8e:36:e2   UHLW        0        5    en0    691
    #     140.221.9.6        0:6:5b:f:51:d6     UHLW        1       63    en0   1197
    #     140.221.10.135     0:d0:59:34:26:34   UHLW        2     2134    en0   1199
    #     140.221.10.152     0:30:1b:b0:ec:dd   UHLW        1      137    en0   1122
    #     140.221.10.153     127.0.0.1          UHS         0        0    lo0
    #     140.221.11.37      0:9:6b:53:4e:4b    UHLW        1      624    en0   1136
    #     140.221.11.103     0:30:48:22:59:e6   UHLW        3      973    en0   1016
    #     140.221.11.224     0:a:95:6f:7:10     UHLW        1        1    en0    605
    #     140.221.11.237     0:1:30:b8:80:c0    UHLW        0        0    en0   1158
    #     140.221.11.250     0:1:30:3:1:0       UHLW        0        0    en0   1141
    #     140.221.11.253     0:d0:3:e:70:a      UHLW       13        0    en0   1199
    #     169.254            link#4             UCS         0        0    en0
    #
    #     Internet6:
    #     Destination                       Gateway                       Flags      Netif Expire
    #                                                                     UH          lo0
    #     fe80::%lo0/64                                                   Uc          lo0
    #                                       link#1                        UHL         lo0
    #     fe80::%en0/64                     link#4                        UC          en0
    #     0:a:95:a8:26:68               UHL         lo0
    #     ff01::/32                                                       U           lo0
    #     ff02::%lo0/32                                                   UC          lo0
    #     ff02::%en0/32                     link#4                        UC          en0

    my($fh);

    if (!open($fh, "netstat -rn |")) {
        warn "Cannot run netstat to determine local IP address\n";
        return "localhost";
    }

    my $interface_name;

    while (<$fh>) {
        my @cols = split();

        if ($cols[0] eq "default" || $cols[0] eq "0.0.0.0") {
            $interface_name = $cols[$#cols];
        }
    }
    close($fh);

    # print "Default route on $interface_name\n";

    #
    # Find ifconfig.
    #

    my $ifconfig;

    for my $dir ((split(":", $ENV{PATH}), "/sbin", "/usr/sbin")) {
        if (-x "$dir/ifconfig") {
            $ifconfig = "$dir/ifconfig";
            last;
        }
    }

    if ($ifconfig eq "") {
        warn "Ifconfig not found\n";
        return "localhost";
    }
    # print "Foudn $ifconfig\n";

    if (!open($fh, "$ifconfig $interface_name |")) {
        warn "Could not run $ifconfig: $!\n";
        return "localhost";
    }

    my $ip;
    while (<$fh>) {
        #
        # Mac:
        #         inet 140.221.10.153 netmask 0xfffffc00 broadcast 140.221.11.255
        # Linux:
        #           inet addr:140.221.34.37  Bcast:140.221.34.63  Mask:255.255.255.224
        #

        chomp;
        s/^\s*//;

        # print "Have '$_'\n";
        if (/inet\s+addr:(\d+\.\d+\.\d+\.\d+)\s+/) {
            #
            # Linux hit.
            #
            $ip = $1;
            # print "Got linux $ip\n";
            last;
        } elsif (/inet\s+(\d+\.\d+\.\d+\.\d+)\s+/) {
            #
            # Mac hit.
            #
            $ip = $1;
            # print "Got mac $ip\n";
            last;
        }
    }
    close($fh);

    if ($ip eq "") {
        warn "Didn't find an IP\n";
        return "localhost";
    }

    return $ip;
}

=head3 get_seed_id

    my $id = FIG::get_seed_id();

Return the Universally Unique ID for this SEED instance. If one
does not exist, it will be created.

=cut
#: Return type $;
sub get_seed_id {
    #
    # Retrieve the seed identifer from FIGdisk/config/seed_id.
    #
    # If it's not there, create one, and make it readonly.
    #
    my $id;
    my $id_file = "$FIG_Config::fig_disk/config/seed_id";
    if (! -f $id_file) {
        my $newid = `uuidgen`;
        if (!$newid) {
            die "Cannot run uuidgen: $!";
        }

        chomp($newid);
        my $fh = new FileHandle(">$id_file");
        if (!$fh) {
            die "error creating $id_file: $!";
        }
        print $fh "$newid\n";
        $fh->close();
        chmod(0444, $id_file);
    }
    my $fh = new FileHandle("<$id_file");
    $id = <$fh>;
    chomp($id);
    return $id;
}

=head3 get_release_info

    my ($name, $id, $inst, $email, $parent_id, $description) = FIG::get_release_info();

Return the current data release information..

The release info comes from the file FIG/Data/RELEASE. It is formatted as:

 <release-name>
 <unique id>
 <institution>
 <contact email>
 <unique id of data release this release derived from>
 <description>

For instance:

 -----
 SEED Data Release, 09/15/2004.
 4148208C-1DF2-11D9-8417-000A95D52EF6
 ANL/FIG
 olson@mcs.anl.gov

 Test release.
 -----

If no RELEASE file exists, this routine will create one with a new unique ID. This
lets a peer optimize the data transfer by being able to cache ID translations
from this instance.

=cut
#: Return Type @;
sub get_release_info {
    my($fig, $no_create) = @_;

    my $rel_file = "$FIG_Config::data/RELEASE";

    if (! -f $rel_file and !$no_create) {
        #
        # Create a new one.
        #

        my $newid = `uuidgen`;
        if (!$newid) {
            die "Cannot run uuidgen: $!";
        }

        chomp($newid);

        my $relinfo = "Automatically generated release info " . localtime();
        my $inst = "Unknown";
        my $contact = "Unknown";
        my $parent = "";
        my( $a, $b, $e, $v, $env ) = $fig->genome_counts;
        my $description = "Automatically generated release info\n";
        $description .= "Contains $a archaeal, $b bacterial, $e eukaryal, $v viral and $env environmental genomes.\n";

        my $fh = new FileHandle(">$rel_file");
        if (!$fh) {
            warn "error creating $rel_file: $!";
            return undef;
        }
        print $fh "$relinfo\n";
        print $fh "$newid\n";
        print $fh "$inst\n";
        print $fh "$contact\n";
        print $fh "$parent\n";
        print $fh $description;
        $fh->close();
        chmod(0444, $rel_file);
    }

    if (open(my $fh, $rel_file)) {
        my(@lines) = <$fh>;
        close($fh);

        chomp(@lines);

        my($info, $id, $inst, $contact, $parent, @desc) = @lines;

        return ($info, $id, $inst, $contact, $parent, join("\n", @desc));
    }

    return undef;
}

=head3 Title

    my $title = $fig->Title();

Return the title of this database. For SEED, this will return SEED, for Sprout
it will return NMPDR, and so forth.

=cut

sub Title {
    return "SEED";
}

=head3 FIG

    my $realFig = $fig->FIG();

Return this object. This method is provided for compatability with SFXlate.

=cut

sub FIG {
    my ($self) = @_;
    return $self;
}

=head3 get_peer_last_update

    my $date = $fig->get_peer_last_update($peer_id);

Return the timestamp from the last successful peer-to-peer update with
the given peer. If the specified peer has made updates, comparing this
timestamp to the timestamp of the updates can tell you whether or not
the updates have been integrated into your SEED data store.

We store this information in FIG/Data/Global/Peers/<peer-id>.

=over 4

=item peer_id

Universally Unique ID for the desired peer.

=item RETURN

Returns the date/time stamp for the last peer-to-peer updated performed
with the identified SEED instance.

=back

=cut
#: Return Type $;
sub get_peer_last_update  {
    my($self, $peer_id) = @_;

    my $dir = "$FIG_Config::data/Global/Peers";
    &verify_dir($dir);
    $dir .= "/$peer_id";
    &verify_dir($dir);

    my $update_file = "$dir/last_update";
    if (-f $update_file)  {
        my $time = file_head($update_file, 1);
        chomp $time;
        return $time;
    } else {
        return undef;
    }
}

=head3 set_peer_last_update

    $fig->set_peer_last_update($peer_id, $time);

Manually set the update timestamp for a specified peer. This informs
the SEED that you have all of the assignments and updates from a
particular SEED instance as of a certain date.

=cut
#: Return Type ;

sub set_peer_last_update {
    my($self, $peer_id, $time) = @_;

    my $dir = "$FIG_Config::data/Global/Peers";
    &verify_dir($dir);
    $dir .= "/$peer_id";
    &verify_dir($dir);

    my $update_file = "$dir/last_update";
    open(F, ">$update_file");
    print F "$time\n";
    close(F);
}

=head3 clean_spaces

Remove any extra spaces from input fields. This will (currently) remove ^\s, \s$, and concatenate multiple spaces into one.

my $input=$fig->clean_spaces($cgi->param('input'));

=cut

sub clean_spaces
{
 my ($self, $s)=@_;
 # note at the moment I do not use \s because that recognizes \t and \n too. This should only remove multiple spaces.
 $s =~ s/^ +//;
 $s =~ s/ +$//;
 $s =~ s/ +/ /g;
 return $s;
}



=head3 cgi_url

    my $url = FIG::$fig->cgi_url();

Return the URL for the CGI script directory.

=cut
#: Return Type $;
sub cgi_url {
#    return &plug_url($FIG_Config::cgi_url);

    #
    # In order to globally make relative references work properly, return ".".
    # This might break some stuff in p2p, but this will get us most of the way there.
    # The things that break we can repair by inspecting the value of $ENV{SCRIPT_NAME}
    #
    return ".";
}

=head3 top_link

    my $url = FIG::top_link();

Return the relative URL for the top of the CGI script directory.

We determine this based on the SCRIPT_NAME environment variable, falling
back to FIG_Config::cgi_base if necessary.

=cut

sub top_link
{

    #
    # Determine if this is a toplevel cgi or one in one of the subdirs (currently
    # just /p2p).
    #

    my @parts = split(/\//, $ENV{SCRIPT_NAME});
    my $top;
    if ($parts[-2] eq 'FIG')
    {
        $top = '.';
#       warn "toplevel @parts\n";
    }
    elsif ($parts[-3] eq 'FIG')
    {
        $top = '..';
#       warn "subdir @parts\n";
    }
    else
    {
        $top = $FIG_Config::cgi_base;
#       warn "other @parts\n";
    }

    return $top;
}

=head3 temp_url

    my $url = FIG::temp_url();

Return the URL of the temporary file directory.

=cut
#: Return Type $;
sub temp_url {
#    return &plug_url($FIG_Config::temp_url);

    #
    # Similarly, make this relative.
    #
    return "../FIG-Tmp";
}

=head3 plug_url

    my $url2 = $fig->plug_url($url);

or

    my $url2 = $fig->plug_url($url);

Change the domain portion of a URL to point to the current domain. This essentially
relocates URLs into the current environment.

=over 4

=item url

URL to relocate.

=item RETURN

Returns a new URL with the base portion converted to the current operating host.
If the URL does not begin with C<http://>, the URL will be returned unmodified.

=back

=cut
#: Return Type $;
sub plug_url {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($url) = @_;

    my $name;

    #  Revised by GJO
    #  First try to get url from the current http request

    if (      defined( $ENV{ 'HTTP_HOST' } )   # This is where $cgi->url gets its value
         && ( $name =  $ENV{ 'HTTP_HOST' } )
         && ( $url  =~ s~^http://[^/]*~http://$name~ )  # ~ is delimiter
       ) {}

    #  Otherwise resort to alternative sources

    elsif ( ( $name = &get_local_hostname )
         && ( $url  =~ s~^http://[^/]*~http://$name~ )  # ~ is delimiter
       ) {}

    return $url;
}

=head3 file_read

    my $text = $fig->file_read($fileName);

or

    my @lines = $fig->file_read($fileName);

or

    my $text = FIG::file_read($fileName);

or

    my @lines = FIG::file_read($fileName);

Read an entire file into memory. In a scalar context, the file is returned
as a single text string with line delimiters included. In a list context, the
file is returned as a list of lines, each line terminated by a line
delimiter. (For a method that automatically strips the line delimiters,
use C<Tracer::GetFile>.)

=over 4

=item fileName

Fully-qualified name of the file to read.

=item RETURN

In a list context, returns a list of the file lines. In a scalar context, returns
a string containing all the lines of the file with delimiters included.

=back

=cut
#: Return Type $;
#: Return Type @;
sub file_read {

    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($fileName) = @_;
    return file_head($fileName, '*');

}


=head3 file_head

    my $text = $fig->file_head($fileName, $count);

or

    my @lines = $fig->file_head($fileName, $count);

or

    my $text = FIG::file_head($fileName, $count);

or

    my @lines = FIG::file_head($fileName, $count);

Read a portion of a file into memory. In a scalar context, the file portion is
returned as a single text string with line delimiters included. In a list
context, the file portion is returned as a list of lines, each line terminated
by a line delimiter.

=over 4

=item fileName

Fully-qualified name of the file to read.

=item count (optional)

Number of lines to read from the file. If omitted, C<1> is assumed. If the
non-numeric string C<*> is specified, the entire file will be read.

=item RETURN

In a list context, returns a list of the desired file lines. In a scalar context, returns
a string containing the desired lines of the file with delimiters included.

=back

=cut
#: Return Type $;
#: Return Type @;
sub file_head {

    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($file, $count) = @_;

    my ($n, $allFlag);
    if ($count eq '*') {
        Trace("Full file read for \"$file\".") if T(3);
        $allFlag = 1;
        $n = 0;
    } else {
        $allFlag = 0;
        $n = (!$count ? 1 : $count);
        Trace("Reading $n record(s) from \"$file\".") if T(3);
    }

    if (open(my $fh, "<$file")) {
        my(@ret, $i);
        $i = 0;
        while (<$fh>) {
            push(@ret, $_);
            $i++;
            last if !$allFlag && $i >= $n;
        }
        close($fh);
        if (wantarray) {
            return @ret;
        } else {
            return join("", @ret);
        }
    }
}

################ Basic Routines [ existed since WIT ] ##########################

=head3 min

    my $min = FIG::min(@x);

or

    my $min = $fig->min(@x);

Return the minimum numeric value from a list.

=over 4

=item x1, x2, ... xN

List of numbers to process.

=item RETURN

Returns the numeric value of the list entry possessing the lowest value. Returns
C<undef> if the list is empty, or has no defined values.

=back

=cut
#: Return Type $;
sub min {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);

    shift while ( @_ && ! defined( $_[0] ) );
    return undef unless @_;

    my $min = shift;
    foreach ( @_ ) { $min = $_ if defined($_) && ( $_ < $min ) }

    return $min;
}

=head3 max

    my $max = FIG::max(@x);

or

    my $max = $fig->max(@x);

Return the maximum numeric value from a list.

=over 4

=item x1, x2, ... xN

List of numbers to process.

=item RETURN

Returns the numeric value of t/he list entry possessing the highest value. Returns
C<undef> if the list is empty, or has no defined values.

=back

=cut
#: Return Type $;
sub max {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);

    shift while ( @_ && ! defined( $_[0] ) );
    return undef unless @_;

    my $max = shift;
    foreach ( @_ ) { $max = $_ if defined($_) && ( $_ > $max ) }

    return $max;
}

=head3 between

    my $flag = FIG::between($x, $y, $z);

or

    my $flag = $fig->between($x, $y, $z);

Determine whether or not $y is between $x and $z.

=over 4

=item x

First edge number.

=item y

Number to examine.

=item z

Second edge number.

=item RETURN

Return TRUE if the number I<$y> is between the numbers I<$x> and I<$z>. The check
is inclusive (that is, if I<$y> is equal to I<$x> or I<$z> the function returns
TRUE), and the order of I<$x> and I<$z> does not matter. If I<$x> is lower than
I<$z>, then the return is TRUE if I<$x> <= I<$y> <= I<$z>. If I<$z> is lower,
then the return is TRUE if I<$x> >= I$<$y> >= I<$z>.

=back

=cut
#: Return Type $;
sub between {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($x,$y,$z) = @_;

    if ($x < $z) {
        return (($x <= $y) && ($y <= $z));
    } else {
        return (($x >= $y) && ($y >= $z));
    }
}


=head3 get_organism_info_from_ncbi

C<< my $code = FIG::get_organism_info_from_ncbi( $taxonomyID ); >>

For a given taxonomy ID returns a hash containing scientific name , genetic code , synonyms and lineage

=cut

# originally by Andreas, rewritten to use XML::Simple by /gdp
sub get_organism_info_from_ncbi{
    my ($self , $tax_id) = @_;
    my $overview = {};
    
    #query url
    my $url="http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&report=xml&id=". $tax_id;
    my $content = get($url);
    
    use XML::Simple;
    my $xml = XML::Simple->new();
    my $parsed =  $xml->XMLin($content);
    print STDERR Dumper($parsed) if $ENV{FIG_VERBOSE};
    
    # get genetic code
    my $genetic_code = "";
    if (defined($genetic_code = $parsed->{Taxon}->{GeneticCode}->{GCId})) {
	print STDERR "genetic_code = $genetic_code\n" if $ENV{FIG_VERBOSE};
	$overview->{ genetic_code } = $genetic_code;
    }
    
    # set scientific name
    my $scientific_name = "";
    if (defined($parsed->{Taxon}->{ScientificName})) {
	$scientific_name =  $parsed->{Taxon}->{ScientificName};
        $scientific_name =~ s/Candidatus\s+//go;
	$scientific_name =~ s/^\s*//o;
        $scientific_name =~ s/\s+/ /go;
        $scientific_name =~ s/\s+$//o;
	
	print STDERR "scientific_name = $scientific_name\n" if $ENV{FIG_VERBOSE};
	$overview->{ scientific_name } = $scientific_name;
    }
    
    # set synonyms
    my @alternate_names = ();
    if (defined($parsed->{Taxon}->{OtherNames})) {
	if (defined($parsed->{Taxon}->{OtherNames}->{EquivalentName})) {
	    my $alt_names = $parsed->{Taxon}->{OtherNames}->{EquivalentName};
	    my @alt_names = (ref($alt_names) eq q(ARRAY)) ? @$alt_names : ($alt_names);
	    push @alternate_names, @alt_names;
	}
	
	if (defined($parsed->{Taxon}->{OtherNames}->{Synonym})) {
	    my $synonyms = $parsed->{Taxon}->{OtherNames}->{Synonym};
	    my @synonyms = (ref($synonyms) eq q(ARRAY)) ? @$synonyms : ($synonyms);
	    push @alternate_names, @synonyms;
	}
	
	map { s/^\s+//o;
	      s/^Candidatus\s+//o;
	      s/\s+/ /go;
	      s/\s+$//o;
	      $overview->{synonyms}->{$_} = 1 
	    } @alternate_names;
    }
    
    
    #get lineage
    my $lineage = "";
	if (defined($lineage = $parsed->{Taxon}->{Lineage})) {
	$lineage =~ s/^\s+//o;
	$lineage =~ s/^cellular\s+organisms\;\s+//o;
	$lineage =~ s/Candidatus\s+//go;
	$lineage =~ s/\s+/ /go;
	$lineage =~ s/\s+$//o;
	
	print STDERR "lineage = \'$lineage\'\n" if $ENV{FIG_VERBOSE};
	$overview->{ lineage } = $lineage;
    }
    
#...Extract Genus, species, strain from LineageEx
    my $genus   = "";
    my $species = "";
    my $lineageEx = $parsed->{Taxon}->{LineageEx}->{Taxon};
    foreach my $field (@$lineageEx) {
	if ($field->{Rank} eq q(genus)) {
	    # set genus
	    $genus =  $field->{ScientificName};
	    
	    $genus =~ s/^\s+//o;
	    $genus =~ s/^Candidatus\s+//o;
	    $genus =~ s/\s+$//o;
	    
	    $overview->{ genus } = $genus ;
	}
	
	if ($field->{Rank} eq q(species)) {
	    # set species
	    $species =  $field->{ScientificName};
	    
	    $species =~ s/^\s+//o;
	    $species =~ s/^Candidatus\s+//o;
	    $species =~ s/$genus\s+//;
	    $species =~ s/\s+$//o;
	    #$species =~ s/ii$/i/;
	    #$species =~ s/ae$/a/;
	    
	    $overview->{ species } = $species ;
	}
    }
    
    # set strain
    my $strain = "";
    if ($strain = $scientific_name) {
	$strain =~ s/^$genus\s+//;
	$strain =~ s/^$species\s+//;
	$overview->{ strain } = $strain;
    }
    
    return $overview;
}


# The above routine parses strings out of xml that can, and does, include
# escaped characters.  We need to convert to plain text. -- GJO
#
my %named_char = ( quot => '"', amp => '&', lt => '<', gt => '>', apos => "'" );

sub decode_html_chars
{
    join '', map { /&#(\d+);/      && ( $1 < 256 )         ? chr( $1 )
                 : /&([a-zA-Z]+);/ && $named_char{ lc $1 } ? $named_char{ lc $1 }
                 : $_
                 } split /(&[a-zA-Z]+|#\d+;)/, shift;
}




=head3 standard_genetic_code

    my $code = FIG::standard_genetic_code();

Return a hash containing the standard translation of nucleotide triples to proteins.
Methods such as L</translate> can take a translation scheme as a parameter. This method
returns the default translation scheme. The scheme is implemented as a reference to a
hash that contains nucleotide triplets as keys and has protein letters as values.

=cut

sub genetic_code {
    my ($ncbi_genetic_code_num) = @_;
    my $code = &standard_genetic_code();

    if (($ncbi_genetic_code_num ==  1) ||
	($ncbi_genetic_code_num == 11)
	) {
	#...Do nothing
    }
    elsif ($ncbi_genetic_code_num ==  4) {
	$code->{TGA} = 'W';
    }
    else {
	die "Sorry, only genetic codes 1, 4, and 11 are currently supported";
    }

    return $code;
}

#: Return Type $;
sub standard_genetic_code {

    my $code = {};

    $code->{"AAA"} = "K";
    $code->{"AAC"} = "N";
    $code->{"AAG"} = "K";
    $code->{"AAT"} = "N";
    $code->{"ACA"} = "T";
    $code->{"ACC"} = "T";
    $code->{"ACG"} = "T";
    $code->{"ACT"} = "T";
    $code->{"AGA"} = "R";
    $code->{"AGC"} = "S";
    $code->{"AGG"} = "R";
    $code->{"AGT"} = "S";
    $code->{"ATA"} = "I";
    $code->{"ATC"} = "I";
    $code->{"ATG"} = "M";
    $code->{"ATT"} = "I";
    $code->{"CAA"} = "Q";
    $code->{"CAC"} = "H";
    $code->{"CAG"} = "Q";
    $code->{"CAT"} = "H";
    $code->{"CCA"} = "P";
    $code->{"CCC"} = "P";
    $code->{"CCG"} = "P";
    $code->{"CCT"} = "P";
    $code->{"CGA"} = "R";
    $code->{"CGC"} = "R";
    $code->{"CGG"} = "R";
    $code->{"CGT"} = "R";
    $code->{"CTA"} = "L";
    $code->{"CTC"} = "L";
    $code->{"CTG"} = "L";
    $code->{"CTT"} = "L";
    $code->{"GAA"} = "E";
    $code->{"GAC"} = "D";
    $code->{"GAG"} = "E";
    $code->{"GAT"} = "D";
    $code->{"GCA"} = "A";
    $code->{"GCC"} = "A";
    $code->{"GCG"} = "A";
    $code->{"GCT"} = "A";
    $code->{"GGA"} = "G";
    $code->{"GGC"} = "G";
    $code->{"GGG"} = "G";
    $code->{"GGT"} = "G";
    $code->{"GTA"} = "V";
    $code->{"GTC"} = "V";
    $code->{"GTG"} = "V";
    $code->{"GTT"} = "V";
    $code->{"TAA"} = "*";
    $code->{"TAC"} = "Y";
    $code->{"TAG"} = "*";
    $code->{"TAT"} = "Y";
    $code->{"TCA"} = "S";
    $code->{"TCC"} = "S";
    $code->{"TCG"} = "S";
    $code->{"TCT"} = "S";
    $code->{"TGA"} = "*";
    $code->{"TGC"} = "C";
    $code->{"TGG"} = "W";
    $code->{"TGT"} = "C";
    $code->{"TTA"} = "L";
    $code->{"TTC"} = "F";
    $code->{"TTG"} = "L";
    $code->{"TTT"} = "F";

    return $code;
}

sub trans_tab {
    my($code) = @_;

    my $tt = &FIG::standard_genetic_code;
    if ($code == 4)
    {
	$tt->{'TGA'} = "W";
    }
    return $tt;
}

sub fr_to_go {
    my($self,$role) = @_;

    my $roleQ = quotemeta $role;
    my $rdbH = $self->db_handle;
    my $relational_db_response = $rdbH->SQL("SELECT go_id  FROM fr2go WHERE role = '$roleQ'");
    return map { $_->[0] } @{$relational_db_response};
}

=head3 translate

    my $aa_seq = &FIG::translate($dna_seq, $code, $fix_start);

Translate a DNA sequence to a protein sequence using the specified genetic code.
If I<$fix_start> is TRUE, will translate an initial C<TTG> or C<GTG> code to
C<M>. (In the standard genetic code, these two combinations normally translate
to C<V> and C<L>, respectively.)

=over 4

=item dna_seq

DNA sequence to translate. Note that the DNA sequence can only contain
known nucleotides.

=item code

Reference to a hash specifying the translation code. The hash is keyed by
nucleotide triples, and the value for each key is the corresponding protein
letter. If this parameter is omitted, the L</standard_genetic_code> will be
used.

=item fix_start

TRUE if the first triple is to get special treatment, else FALSE. If TRUE,
any value in the first position will be translated as C<M>.

=item RETURN

Returns a string resulting from translating each nucleotide triple into a
protein letter.

=back

=cut
#: Return Type $;
sub translate {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);

    my( $dna,$code,$start ) = @_;
    my( $i,$j,$ln );
    my( $x,$y );
    my( $prot );

    if (! defined($code)) {
        $code = &FIG::standard_genetic_code;
    }
    $ln = length($dna);
    $prot = "X" x ($ln/3);
    $dna =~ tr/a-z/A-Z/;

    for ($i=0,$j=0; ($i < ($ln-2)); $i += 3,$j++) {
        $x = substr($dna,$i,3);
        if ($y = $code->{$x}) {
            substr($prot,$j,1) = $y;
        }
    }

    if ($start) {
        substr($prot,0,1) = 'M';
    }
    return $prot;
}

=head3 reverse_comp

    my $dnaR = FIG::reverse_comp($dna);

or

    my $dnaR = $fig->reverse_comp($dna);

Return the reverse complement os the specified DNA sequence.

NOTE: for extremely long DNA strings, use L</rev_comp>, which allows you to
pass the strings around in the form of pointers.

=over 4

=item dna

DNA sequence whose reverse complement is desired.

=item RETURN

Returns the reverse complement of the incoming DNA sequence.

=back

=cut
#: Return Type $;
sub reverse_comp {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($seq) = @_;

    return ${&rev_comp(\$seq)};
}

=head3 rev_comp

    my $dnaRP = FIG::rev_comp(\$dna);

or

    my $dnaRP = $fig->rev_comp(\$dna);

Return the reverse complement of the specified DNA sequence. The DNA sequence
is passed in as a string reference rather than a raw string for performance
reasons. If this is unnecessary, use L</reverse_comp>, which processes strings
instead of references to strings.

=over 4

=item dna

Reference to the DNA sequence whose reverse complement is desired.

=item RETURN

Returns a reference to the reverse complement of the incoming DNA sequence.

=back

=cut
#: Return Type $;
sub rev_comp {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my( $seqP ) = @_;
    my( $rev  );

    $rev =  reverse( $$seqP );
    $rev =~ tr/A-Z/a-z/;
    $rev =~ tr/acgtumrwsykbdhv/tgcaakywsrmvhdb/;
    return \$rev;
}

# This routine was written by Gary to definitively handle the "scratch" subdirectory issue.
# It takes as parameters key-value pairs.  The relevant ones are
#
#     tmpdir => NameOfTmpDirectoryToBeUsed  [can be ommitted]
#     tmp    => TheNameOfTheTmpDirectoryToContainTheSubdirectory [can be ommitted]
#
# if tmpdir exists, save_tmp is set to "true".  You need to test this at the end
# of your script and blow away the directory unless save_tmp is true.
# if tmpdir does not exist, it will be created if possible.
#
# tmp is where to put tmpdir, if it is not specified.  if tmp is omitted, it
# will all be ok.
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  ( $tmp_dir, $save_tmp ) = temporary_directory( \%options )
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub temporary_directory
{
    my $options = shift;

    my $tmp_dir  = $options->{ tmpdir };
    my $save_tmp = $options->{ savetmp } || '';

    if ( $tmp_dir )
    {
        if ( -d $tmp_dir ) { $options->{ savetmp } = $save_tmp = 1 }
    }
    else
    {
        my $tmp = $options->{ tmp } && -d  $options->{ tmp } ?  $options->{ tmp }
                : $FIG_Config::temp && -d  $FIG_Config::temp ?  $FIG_Config::temp
                :                      -d '/tmp'             ? '/tmp'
                :                                              '.';
	$tmp_dir = sprintf( "$tmp/fig_tmp_dir.%05d.%09d", $$, int(1000000000*rand) );
    }

    if ( $tmp_dir && ! -d $tmp_dir )
    {
        mkdir $tmp_dir;
        if ( ! -d $tmp_dir )
        {
            print STDERR "FIG::temporary_directory could not create '$tmp_dir: $!'\n";
            $options->{ tmpdir } = $tmp_dir = undef;
        }
    }

    return ( $tmp_dir, $save_tmp );
}

sub verify_external_tool {
    my(@progs) = @_;

    my $prog;
    foreach $prog (@progs)
    {
        my @tmp = `which $prog`;
        if ($tmp[0] =~ /^no $prog/)
        {
            print STDERR $tmp[0];
            exit(1);
        }
    }
}

=head3 verify_dir

    FIG::verify_dir($dir);

or

    $fig->verify_dir($dir);

Insure that the specified directory exists.  If it must be created, the permissions will
be set to C<0777>.

=cut

sub verify_dir {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($dir) = @_;

    if (!defined($dir))
    {
        Confess("FIG::verify_dir: missing \$dir argument\n");
    }
    if ($dir eq "")
    {
        confess("FIG::verify_dir: refusing to create a directory named ''\n");
    }

    if (-d $dir) {
        return
    }
    if ($dir =~ /^(.*)\/[^\/]+$/ and $1 ne '') {
        &verify_dir($1);
    }
    if (!mkdir($dir,0777) && $! != Errno::EEXIST)
    {
	confess "Could not make directory $dir: $!";
    }
}

=head3 run

    FIG::run($cmd);

or

    $fig->run($cmd);

Run a command. If the command fails, the error will be traced.

=cut

sub run {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($cmd) = @_;

    if ($ENV{FIG_VERBOSE}) {
        my @tmp = `date`;
        chomp @tmp;
        print STDERR "$tmp[0]: running $cmd\n";
    }
    Trace("Running command: $cmd") if T(3);
    (system($cmd) == 0) || Confess("FAILED: $cmd");
}

=head3 run_gathering_output

    FIG::run_gathering_output($cmd, @args);

or

    $fig->run_gathering_output($cmd, @args);

Run a command, gathering the output. This is similar to the backtick
operator, but it does not invoke the shell. Note that the argument list
must be explicitly passed one command line argument per argument to
run_gathering_output.

If the command fails, the error will be traced.

=cut

sub run_gathering_output {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($cmd, @args) = @_;

    #
    # Run the command in a safe fork-with-pipe/exec.
    #

    my $pid = open(PROC_READ, "-|");

    if ($pid == 0)
    {
        exec { $cmd } $cmd, @args;
        die "could not execute $cmd @args: $!\n";
    }

    if (wantarray)
    {
        my @out;
        while (<PROC_READ>)
        {
            push(@out, $_);
        }
        if (!close(PROC_READ))
        {
            Confess("FAILED: $cmd @args with error return $?");
        }
        return @out;
    }
    else
    {
        my $out = '';

        while (<PROC_READ>)
        {
            $out .= $_;
        }
        if (!close(PROC_READ))
        {
            Confess("FAILED: $cmd @args with error return $?");
        }
        return $out;
    }
}

=head3 interpret_error_code

    ($exitcode, $signal, $msg) = &FIG::interpret_error_code($rc);

Determine if the given result code was due to a process exiting abnormally
or by receiving a signal.

=cut

sub interpret_error_code
{
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);

    my($rc) = @_;

    if (WIFEXITED($rc))
    {
	return (WEXITSTATUS($rc),  undef, "Exited with status " . WEXITSTATUS($rc));
    }
    elsif (WIFSIGNALED($rc))
    {
	return (undef, WTERMSIG($rc), "Terminated with signal " .  WTERMSIG($rc));
    }
    elsif (WIFSTOPPED($rc))
    {
	return (undef, WSTOPSIG($rc), "Stopped with signal " .  WSTOPSIG($rc));
    }
    else
    {
	return ($rc, undef, "Unknown return code $rc");
    }
}

=head3 find_fig_executable

C<< $path = FIG::find_fig_executable("index_sims_file") >>

Looks for the given executable first in $FIG_Config::ext_bin, then in
$FIG_Config::bin. Supports code running either in the original SEED
world which had C programs build as part of FigKernelScripts and the new
world which puts them into the common runtime.

=cut

sub find_fig_executable
{
    my($exe) = @_;
    my $path;
    if (-x ($path = "$FIG_Config::ext_bin/$exe"))
    {
	return $path;
    }
    elsif (-x ($path = "$FIG_Config::bin/$exe"))
    {
	return $path;
    }
    else
    {
	cluck "FIG executable '$exe' not found in standard locations";
	return $exe;
    }
}

=head3 augment_path

    FIG::augment_path($dirName);

Add a directory to the system path.

This method adds a new directory to the front of the system path. It looks in the
configuration file to determine whether this is Windows or Unix, and uses the
appropriate separator.

=over 4

=item dirName

Name of the directory to add to the path.

=back

=cut

sub augment_path {
    my ($dirName) = @_;
    if ($FIG_Config::win_mode) {
        $ENV{PATH} = "$dirName;$ENV{PATH}";
    } else {
        $ENV{PATH} = "$dirName:$ENV{PATH}";
    }
}

=head3 read_fasta_record

    my ($seq_id, $seq_pointer, $comment) = FIG::read_fasta_record(\*FILEHANDLE);

or

    my ($seq_id, $seq_pointer, $comment) = $fig->read_fasta_record(\*FILEHANDLE);

Read and parse the next logical record of a FASTA file. A FASTA logical record
consists of multiple lines of text. The first line begins with a C<< > >> symbol
and contains the sequence ID followed by an optional comment. (NOTE: comments
are currently deprecated, because not all tools handle them properly.) The
remaining lines contain the sequence data.

This method uses a trick to smooth its operation: the line terminator character
is temporarily changed to C<< \n> >> so that a single read operation brings in
the entire logical record.

=over 4

=item FILEHANDLE

Open handle of the FASTA file. If not specified, C<STDIN> is assumed.

=item RETURN

If we are at the end of the file, returns C<undef>. Otherwise, returns a
three-element list. The first element is the sequence ID, the second is
a pointer to the sequence data (that is, a string reference as opposed to
as string), and the third is the comment.

=back

=cut
#: Return Type @;
sub read_fasta_record {

    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my ($file_handle) = @_;
    my ($old_end_of_record, $fasta_record, @lines, $head, $sequence, $seq_id, $comment, @parsed_fasta_record);

    if (not defined($file_handle))  { $file_handle = \*STDIN; }

    $old_end_of_record = $/;
    $/ = "\n>";

    if (defined($fasta_record = <$file_handle>)) {
        chomp $fasta_record;
        @lines  =  split( /\n/, $fasta_record );
        $head   =  shift @lines;
        $head   =~ s/^>?//;
        $head   =~ m/^(\S+)/;
        $seq_id = $1;
        if ($head  =~ m/^\S+\s+(.*)$/)  { $comment = $1; } else { $comment = ""; }
        $sequence  =  join( "", @lines );
        @parsed_fasta_record = ( $seq_id, \$sequence, $comment );
    } else {
        @parsed_fasta_record = ();
    }

    $/ = $old_end_of_record;

    return @parsed_fasta_record;
}

=head3 display_id_and_seq

    FIG::display_id_and_seq($id_and_comment, $seqP, $fh);



Display a fasta ID and sequence to the specified open file. This method is designed
to work well with L</read_fasta_sequence> and L</rev_comp>, because it takes as
input a string pointer rather than a string. If the file handle is omitted it
defaults to STDOUT.

The output is formatted into a FASTA record. The first line of the output is
preceded by a C<< > >> symbol, and the sequence is split into 60-character
chunks displayed one per line. Thus, this method can be used to produce
FASTA files from data gathered by the rest of the system.

=over 4

=item id_and_comment

The sequence ID and (optionally) the comment from the sequence's FASTA record.
The ID

=item seqP

Reference to a string containing the sequence. The sequence is automatically
formatted into 60-character chunks displayed one per line.

=item fh

Open file handle to which the ID and sequence should be output. If omitted,
C<\*STDOUT> is assumed.

=back

=cut

sub display_id_and_seq {

    if (UNIVERSAL::isa($_[0],__PACKAGE__)) {
        shift @_;
        #Trace("Invalid call to display_id_and_seq.");
    }

    my( $id, $seqP, $fh ) = @_;

    if (! defined($fh) )  { $fh = \*STDOUT; }

    print $fh ">$id\n";
    &display_seq($seqP, $fh);
}

=head3 display_seq

    FIG::display_seq(\$seqP, $fh);

or

    $fig->display_seq(\$seqP, $fh);

Display a fasta sequence to the specified open file. This method is designed
to work well with L</read_fasta_sequence> and L</rev_comp>, because it takes as
input a string pointer rather than a string. If the file handle is omitted it
defaults to STDOUT.

The sequence is split into 60-character chunks displayed one per line for
readability.

=over 4

=item seqP

Reference to a string containing the sequence.

=item fh

Open file handle to which the sequence should be output. If omitted,
C<STDOUT> is assumed.

=back

=cut

sub display_seq {

    shift if UNIVERSAL::isa($_[0],__PACKAGE__);

    my ( $seqP, $fh ) = @_;
    my ( $i, $n, $ln );

    if (! defined($fh) )  { $fh = \*STDOUT; }

    $n = length($$seqP);
#   confess "zero-length sequence ???" if ( (! defined($n)) || ($n == 0) );
    for ($i=0; ($i < $n); $i += 60) {
        if (($i + 60) <= $n) {
            $ln = substr($$seqP,$i,60);
        } else {
            $ln = substr($$seqP,$i,($n-$i));
        }
        print $fh "$ln\n";
    }
}


=head3 flatten_dumper

    FIG::flatten_dumper( $perl_ref_or_object_1, ... );

    $fig->flatten_dumper( $perl_ref_or_object_1, ... );

Takes a list of perl references or objects, and "flattens" their Data::Dumper() output
so that it can be printed on a single line.

=cut

sub flatten_dumper {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my @x = @_;
    my $x;

    foreach $x (@x)
    {
	$x = Dumper($x);

	$x =~ s/\$VAR\d+\s+\=\s+//o;
	$x =~ s/\n//gso;
	$x =~ s/\s+/ /go;
	$x =~ s/\'//go;
#       $x =~ s/^[^\(\[\{]+//o;
#       $x =~ s/[^\)\]\}]+$//o;
    }

    return @x;
}


##########  I commented the pods on the following routines out, since they should not
##########  be part of the SOAP/WSTL interface
#=pod
#
#=head3 file2N
#
#usage: $n = $fig->file2N($file)
#
#In some of the databases I need to store filenames, which can waste a lot of
#space.  Hence, I maintain a database for converting filenames to/from integers.
#
#=cut
#
sub file2N :Scalar {
    my($self,$file) = @_;
    my($relational_db_response);

    my $rdbH = $self->db_handle;

    #
    # Strip the figdisk path from the file. N2file replaces it if the path
    # in the database is relative.
    #
    $file =~ s,^$FIG_Config::fig_disk/,,;

    if (($relational_db_response = $rdbH->SQL("SELECT fileno FROM file_table WHERE ( file = \'$file\')")) &&
        (@$relational_db_response == 1)) {
        return $relational_db_response->[0]->[0];
    } elsif (($relational_db_response = $rdbH->SQL("SELECT MAX(fileno) FROM file_table "))  && (@$relational_db_response == 1) && ($relational_db_response->[0]->[0])) {
        my $fileno = $relational_db_response->[0]->[0] + 1;
        if ($rdbH->SQL("INSERT INTO file_table ( file, fileno ) VALUES ( \'$file\', $fileno )")) {
            return $fileno;
        }
    } elsif ($rdbH->SQL("INSERT INTO file_table ( file, fileno ) VALUES ( \'$file\', 1 )")) {
        return 1;
    }
    return undef;
}

#=pod
#
#=head3 N2file
#
#usage: $filename = $fig->N2file($n)
#
#In some of the databases I need to store filenames, which can waste a lot of
#space.  Hence, I maintain a database for converting filenames to/from integers.
#
#=cut
#
sub N2file :Scalar
{
    my($self,$fileno) = @_;

    #
    # Cache outputs. This results in a huge savings of time when files are
    # accessed multiple times (as in when a bunch of sims are requested).
    #

    my $fcache = $self->cached("_n2file");

    my $fname;
    if (defined($fname = $fcache->{$fileno}))
    {
        return $fname;
    }

    my $rdbH = $self->db_handle;

    my $relational_db_response = $rdbH->SQL("SELECT file FROM file_table WHERE ( fileno = $fileno )");

    if ($relational_db_response and @$relational_db_response == 1)
    {
        $fname = $relational_db_response->[0]->[0];

        #
        # If $fname is relative, prepend the base of the fig_disk.
        # (Updated to use PERL's system-independent filename utilities.
        #

        $fname = File::Spec->rel2abs($fname, $FIG_Config::fig_disk);

        $fcache->{$fileno} = $fname;
        return $fname;
    }
    return undef;
}


#=pod
#
#=head3 openF
#
#usage: $fig->openF($filename)
#
#Parts of the system rely on accessing numerous different files.  The most obvious case is
#the situation with similarities.  It is important that the system be able to run in cases in
#which an arbitrary number of files cannot be open simultaneously.  This routine (with closeF) is
#a hack to handle this.  I should probably just pitch them and insist that the OS handle several
#hundred open filehandles.
#
#=cut
#
sub openF {
    my($self,$file) = @_;
    my($fxs,$x,@fxs,$fh);

    $fxs = $self->cached('_openF');
    if ($x = $fxs->{$file}) {
        $x->[1] = time();
        return $x->[0];
    }

    @fxs = keys(%$fxs);
    if (defined($fh = new FileHandle "<$file")) {
        if (@fxs >= 50) {
            @fxs = sort { $fxs->{$a}->[1] <=> $fxs->{$b}->[1] } @fxs;
            $x = $fxs->{$fxs[0]};
            undef $x->[0];
            delete $fxs->{$fxs[0]};
        }
        $fxs->{$file} = [$fh,time()];
        return $fh;
    }
    return undef;
}

#=pod
#
#=head3 closeF
#
#usage: $fig->closeF($filename)
#
#Parts of the system rely on accessing numerous different files.  The most obvious case is
#the situation with similarities.  It is important that the system be able to run in cases in
#which an arbitrary number of files cannot be open simultaneously.  This routine (with openF) is
#a hack to handle this.  I should probably just pitch them and insist that the OS handle several
#hundred open filehandles.
#
#=cut
#
sub closeF {
    my($self,$file) = @_;
    my($fxs,$x);

    if (($fxs = $self->{_openF}) && ($x = $fxs->{$file})) {
        undef $x->[0];
        delete $fxs->{$file};
    }
}

=head3 sapling

    my $sapDB = $fig->sapling();

Return a copy of the L<Sapling> database object. If one has already been
created, it will be re-used. Otherwise, one will be created and cached in
the FIG object.

=cut

sub sapling {
    # Get the parameters.
    my ($self) = @_;
    # Look for the cached object.
    my $retVal = $self->{sapling};
    # Did we find it?
    if (! defined $retVal) {
        # Get access to ERDB.
        require ERDB;
        # Connect to the sapling.
        $retVal = ERDB::GetDatabase('Sapling');
        # Cache it for future use.
        $self->{sapling} = $retVal;
    }
    # Return the result.
    return $retVal;
}

=head3 ec_name

    my $enzymatic_function = $fig->ec_name($ec);

Returns the enzymatic name corresponding to the specified enzyme code.

=over 4

=item ec

Code number for the enzyme whose name is desired. The code number is actually
a string of digits and periods (e.g. C<1.2.50.6>).

=item RETURN

Returns the name of the enzyme specified by the indicated code, or a null string
if the code is not found in the database.

=back

=cut

sub ec_name {
    my($self,$ec) = @_;

    ($ec =~ /^\d+\.\d+\.\d+\.\d+$/) || return "";
    my $rdbH = $self->db_handle;
    my $relational_db_response = $rdbH->SQL("SELECT name FROM ec_names WHERE ( ec = \'$ec\' )");

    return (@$relational_db_response == 1) ? $relational_db_response->[0]->[0] : "";
    return "";
}

=head3 all_roles

    my @roles = $fig->all_roles;

Return a list of the known roles. Currently, this is a list of the enzyme codes and names.

The return value is a list of list references. Each element of the big list contains an
enzyme code (EC) followed by the enzymatic name.

=cut

sub all_roles {
    my($self) = @_;

    my $rdbH = $self->db_handle;
    my $relational_db_response = $rdbH->SQL("SELECT ec,name FROM ec_names");

    return @$relational_db_response;
}

=head3 expand_ec

    my $expanded_ec = $fig->expand_ec($ec);

Expands "1.1.1.1" to "1.1.1.1 - alcohol dehydrogenase" or something like that.

=cut

sub expand_ec {
    my($self,$ec) = @_;
    my($name);

    return ($name = $self->ec_name($ec)) ? "$ec - $name" : $ec;
}

=head3 clean_tmp

    FIG::clean_tmp();

Delete temporary files more than two days old.

We store temporary files in $FIG_Config::temp.  There are specific classes of files
that are created and should be saved for at least a few days.  This routine can be
invoked to clean out those that are over two days old.

=cut

sub clean_tmp {

    my($file);
    if (opendir(TMP,"$FIG_Config::temp")) {
    #       change the pattern to pick up other files that need to be cleaned up
        my @temp = grep { $_ =~ /^(Geno|tmp)/ } readdir(TMP);
        foreach $file (@temp) {
            if (-M "$FIG_Config::temp/$file" > 2) {
                unlink("$FIG_Config::temp/$file");
            }
        }
    }
}

################ Routines to process genomes and genome IDs  ##########################


=head3 genomes

    my @genome_ids = $fig->genomes($complete, $restrictions, $domain);

Return a list of genome IDs. If called with no parameters, all genome IDs
in the database will be returned.

Genomes are assigned ids of the form X.Y where X is the taxonomic id maintained by
NCBI for the species (not the specific strain), and Y is a sequence digit assigned to
this particular genome (as one of a set with the same genus/species).  Genomes also
have versions, but that is a separate issue.

=over 4

=item complete

TRUE if only complete genomes should be returned, else FALSE.

=item restrictions

TRUE if only restriction genomes should be returned, else FALSE.

=item domain

Name of the domain from which the genomes should be returned. Possible values are
C<Bacteria>, C<Virus>, C<Eukaryota>, C<unknown>, C<Archaea>, and
C<Environmental Sample>. If no domain is specified, all domains will be
eligible.

=item RETURN

Returns a list of all the genome IDs with the specified characteristics.

=back

=cut
#: Return Type @;
sub genomes  :Remote :List {
    my( $self, $complete, $restrictions, $domain ) = @_;

    return $self->{gdata}->genomes($complete, $restrictions, $domain);

#   my $rdbH = $self->db_handle;
#
#   my @where = ();
#   if ($complete) {
#       push(@where, "( complete = \'1\' )")
#   }
#
#   if ($restrictions) {
#       push(@where, "( restrictions = \'1\' )")
#   }
#
#   if ($domain) {
#       push( @where, "( maindomain = '$domain' )" )
#   }
#
#   my $relational_db_response;
#   if (@where > 0) {
#       my $where = join(" AND ",@where);
#       $relational_db_response = $rdbH->SQL("SELECT genome  FROM genome where $where");
#   } else {
#       $relational_db_response = $rdbH->SQL("SELECT genome  FROM genome");
#   }
#   my @genomes = sort { $a <=> $b } map { $_->[0] } @$relational_db_response;
#   return @genomes;
}

sub genome_list {
  my( $self ) = @_;

  return $self->{gdata}->genome_list();

#   my $rdbH = $self->db_handle;
#   my $relational_db_response = $rdbH->SQL("SELECT genome, gname, maindomain FROM genome where complete=1");
#
#   return $relational_db_response;
}

=head3 genome_info

    my $info = $fig->genome_info();

Return an array reference of information from the genome table

=over 4

=item RETURN

This will return an array reference of genome table entries. All entries of the table will be
returned. The columns will be the following:

genome, gname, szdna, maindomain, pegs, rnas, complete, taxonomy

=back

=cut

sub genome_info {
    my ($self) = @_;
    return $self->{gdata}->genome_info();
#    my $rdbH = $self->db_handle;
#    return $rdbH->SQL("SELECT genome, gname, szdna, maindomain, pegs, rnas, complete, taxonomy FROM genome");
}

=head3 is_complete

    my $flag = $fig->is_complete($genome);

Return TRUE if the genome with the specified ID is complete, else FALSE.

=over 4

=item genome

ID of the relevant genome.

=item RETURN

Returns TRUE if there is a complete genome in the database with the specified ID,
else FALSE.

=back

=cut

sub is_complete {
    my($self,$genome) = @_;

    return $self->{gdata}->complete($genome);

    #  Previous version
    #
    # my $rdbH = $self->db_handle;
    # my $relational_db_response = $rdbH->SQL("SELECT genome  FROM genome where (genome = '$genome') AND (complete = '1')");
    # return (@$relational_db_response == 1)
}

=head3 is_genome

    my $flag = $fig->is_genome($genome);

Return TRUE if the specified genome exists, else FALSE.

=over 4

=item genome

ID of the genome to test.

=item RETURN

Returns TRUE if a genome with the specified ID exists in the data store, else FALSE.

=back

=cut

sub is_genome {
    my($self,$genome) = @_;
    my($y);

    return defined($self->{gdata}->complete($genome));

    #  Previous version
    #
    # my $is_genome = $self->cached("_is_genome");
    #
    # if (defined($y = $is_genome->{$genome}))
    # {
    #     return $y;
    # }
    #
    # my $rdbH = $self->db_handle;
    # my $relational_db_response = $rdbH->SQL("SELECT genome  FROM genome where (genome = '$genome')");
    # $y = (@$relational_db_response == 1);
    # $is_genome->{$genome} = $y;
    # return $y;
}

=head3 assert_genomes

    $fig->assert_genomes(gid, gid, ...);

Assert that the given list of genomes does exist, and allow is_genome() to succeed for them.

This is used in FIG-based computations in the context of the RAST genome-import code, so that
genomes that currently exist only in RAST are treated as present for the purposes of FIG.pm-based
code.

=cut

sub assert_genomes
{
    my($self, @genomes) = @_;

    my $assert = $self->cached("_is_genome");
    map { $assert->{$_} = 1 } @genomes;
}

=head3 genome_counts

    my ($arch, $bact, $euk, $vir, $env, $unk) = $fig->genome_counts($complete);

Count the number of genomes in each domain. If I<$complete> is TRUE, only complete
genomes will be included in the counts.

=over 4

=item complete

TRUE if only complete genomes are to be counted, FALSE if all genomes are to be
counted

=item RETURN

A six-element list containing the number of genomes in each of six categories--
Archaea, Bacteria, Eukaryota, Viral, Environmental, and Unknown, respectively.

=back

=cut

sub genome_counts {
    my($self,$complete) = @_;
    my($x,$relational_db_response);

#     my $rdbH = $self->db_handle;

#     if ($complete) {
#         $relational_db_response = $rdbH->SQL("SELECT genome, maindomain FROM genome where complete = '1'");
#     } else {
#         $relational_db_response = $rdbH->SQL("SELECT genome,maindomain FROM genome");
#     }

    $relational_db_response = $self->{gdata}->get_results_for_genome_counts($complete);


    my ($arch, $bact, $euk, $vir, $env, $unk) = (0, 0, 0, 0, 0, 0);
    if (@$relational_db_response > 0) {
        foreach $x (@$relational_db_response) {
            if    ($x->[1] =~ /^archaea/i)  { ++$arch }
            elsif ($x->[1] =~ /^bacter/i)   { ++$bact }
            elsif ($x->[1] =~ /^eukar/i)    { ++$euk }
            elsif ($x->[1] =~ /^vir/i)      { ++$vir }
            elsif ($x->[1] =~ /^env/i)      { ++$env }
            else  { ++$unk }
        }
    }

    return ($arch, $bact, $euk, $vir, $env, $unk);
}


=head3 genome_domain

    my $domain = $fig->genome_domain($genome_id);

Find the domain of a genome.

=over 4

=item genome_id

ID of the genome whose domain is desired.

=item RETURN

Returns the name of the genome's domain (archaea, bacteria, etc.), or C<undef> if
the genome is not in the database.

=back

=cut

sub genome_domain {
    my($self,$genome) = @_;
    my $relational_db_response;
    my $rdbH = $self->db_handle;

    if ($genome) {

	return $self->{gdata}->maindomain($genome);

#         if (($relational_db_response = $rdbH->SQL("SELECT genome,maindomain FROM genome WHERE ( genome = \'$genome\' )"))
#             && (@$relational_db_response == 1)) {
#             # die Dumper($relational_db_response);
#             return $relational_db_response->[0]->[1];
#         }
    }
    return undef;
}


=head3 genome_pegs

    my $num_pegs = $fig->genome_pegs($genome_id);

Return the number of protein-encoding genes (PEGs) for a specified
genome.

=over 4

=item genome_id

ID of the genome whose PEG count is desired.

=item RETURN

Returns the number of PEGs for the specified genome, or C<undef> if the genome
is not indexed in the database.

=back

=cut

sub genome_pegs {
    my($self,$genome) = @_;

    return $self->{gdata}->pegs($genome);

    my $relational_db_response;
    my $rdbH = $self->db_handle;

    if ($genome) {
        if (($relational_db_response = $rdbH->SQL("SELECT pegs FROM genome WHERE ( genome = \'$genome\' )"))
            && (@$relational_db_response == 1)) {
            return $relational_db_response->[0]->[0];
        }
    }
    return undef;
}


=head3 genome_rnas

    my $num_rnas = $fig->genome_rnas($genome_id);

Return the number of RNA-encoding genes for a genome.
"$genome_id" is indexed in the "genome" database, and 'undef' otherwise.

=over 4

=item genome_id

ID of the genome whose RNA count is desired.

=item RETURN

Returns the number of RNAs for the specified genome, or C<undef> if the genome
is not indexed in the database.

=back

=cut

sub genome_rnas {
    my($self,$genome) = @_;

    return $self->{gdata}->rnas($genome);

    my $relational_db_response;
    my $rdbH = $self->db_handle;

    if ($genome) {
        if (($relational_db_response = $rdbH->SQL("SELECT rnas FROM genome WHERE ( genome = \'$genome\' )"))
            && (@$relational_db_response == 1)) {
            return $relational_db_response->[0]->[0];
        }
    }
    return undef;
}


=head3 genome_szdna

    my $szdna = $fig->genome_szdna($genome_id);

Return the number of DNA base-pairs in a genome's contigs.

=over 4

=item genome_id

ID of the genome whose base-pair count is desired.

=item RETURN

Returns the number of base pairs in the specified genome's contigs, or C<undef>
if the genome is not indexed in the database.

=back

=cut

sub genome_szdna {
    my($self,$genome) = @_;
    my $relational_db_response;
    my $rdbH = $self->db_handle;

    if ($genome) {
        if (($relational_db_response =
            $rdbH->SQL("SELECT szdna FROM genome WHERE ( genome = \'$genome\' )"))
            && (@$relational_db_response == 1)) {

            return $relational_db_response->[0]->[0];

        }
    }
    return undef;
}

=head3 genome_version

    my $version = $fig->genome_version($genome_id);

Return the version number of the specified genome.

Versions are incremented for major updates.  They are put in as major
updates of the form 1.0, 2.0, ...

Users may do local "editing" of the DNA for a genome, but when they do,
they increment the digits to the right of the decimal.  Two genomes remain
comparable only if the versions match identically.  Hence, minor updating should be
committed only by the person/group responsible for updating that genome.

We can, of course, identify which genes are identical between any two genomes (by matching
the DNA or amino acid sequences).  However, the basic intent of the system is to
support editing by the main group issuing periodic major updates.

=over 4

=item genome_id

ID of the genome whose version is desired.

=item RETURN

Returns the version number of the specified genome, or C<undef> if the genome is not in
the data store or no version number has been assigned.

=back

=cut

sub genome_version :Scalar {
    my($self,$genome) = @_;

    my(@tmp);
    if ((-s "$FIG_Config::organisms/$genome/VERSION") &&
        (@tmp = `cat $FIG_Config::organisms/$genome/VERSION`) &&
        ($tmp[0] =~ /^(\S+)$/)) {
        return $1;
    }
    return undef;
}

=head3 genome_md5sum

    my $md5sum = $fig->genome_md5sum($genome_id);

Returns the MD5 checksum of the specified genome.

The checksum of a genome is defined as the checksum of its signature file. The signature
file consists of tab-separated lines, one for each contig, ordered by the contig id.
Each line contains the contig ID, the length of the contig in nucleotides, and the
MD5 checksum of the nucleotide data, with uppercase letters forced to lower case.

The checksum is indexed in the database. If you know a genome's checksum, you can use
the L</genome_with_md5sum> method to find its ID in the database.

=over 4

=item genome

ID of the genome whose checksum is desired.

=item RETURN

Returns the specified genome's checksum, or C<undef> if the genome is not in the
database.

=back

=cut

sub genome_md5sum :Scalar {
    my($self,$genome) = @_;
    my $relational_db_response;
    my $rdbH = $self->db_handle;

    if ($genome) {
        if (($relational_db_response =
             $rdbH->SQL("SELECT md5sum FROM genome_md5sum WHERE ( genome = \'$genome\' )"))
            && (@$relational_db_response == 1)) {
            return $relational_db_response->[0]->[0];
        }
    }
    return undef;
}

=head3 genome_with_md5sum

    my $genome = $fig->genome_with_md5sum($cksum);

Find a genome with the specified checksum.

The MD5 checksum is computed from the content of the genome (see L</genome_md5sum>). This method
can be used to determine if a genome already exists for a specified content.

=over 4

=item cksum

Checksum to use for searching the genome table.

=item RETURN

The ID of a genome with the specified checksum, or C<undef> if no such genome exists.

=back

=cut

sub genome_with_md5sum :Scalar {
    my($self,$cksum) = @_;
    my $relational_db_response;
    my $rdbH = $self->db_handle;

    if (($relational_db_response =
         $rdbH->SQL("SELECT genome FROM genome_md5sum WHERE ( md5sum = \'$cksum\' )"))
         && (@$relational_db_response == 1)) {
        return $relational_db_response->[0]->[0];
    }

    return undef;
}

=head3 contig_md5sum

    my $cksum = $fig->contig_md5sum($genome, $contig);

Return the MD5 checksum for a contig. The MD5 checksum is computed from the content
of the contig. This method retrieves the checksum stored in the database. The checksum
can be compared to the checksum of an external contig as a cheap way of seeing if they
match.

=over 4

=item genome

ID of the genome containing the contig.

=item contig

ID of the relevant contig.

=item RETURN

Returns the checksum of the specified contig, or C<undef> if the contig is not in the
database.

=back

=cut

sub contig_md5sum :Scalar {
    my($self, $genome, $contig) = @_;
    my $relational_db_response;
    my $rdbH = $self->db_handle;

    if ($genome) {
        if (($relational_db_response =
             $rdbH->SQL(qq(SELECT md5 FROM contig_md5sums WHERE (genome = ? AND contig = ?)), undef, $genome, $contig))
             && (@$relational_db_response == 1)) {
            return $relational_db_response->[0]->[0];
        }
    }
    return undef;
}


# returns all contigs for a given md5sum
sub md5sum_to_contig_genome :Scalar {
    my($self, $md5) = @_;
    my $relational_db_response;
    my $rdbH = $self->db_handle;
    my $error = 0;

    if ($md5) {
        if (($relational_db_response =
             $rdbH->SQL(qq(SELECT contig, genome FROM contig_md5sums WHERE md5 = ?), undef, $md5))
             && (@$relational_db_response == 1)) {
            return @{$relational_db_response->[0]};
	    # return $relational_db_response->[0]->[0];
        }
 	elsif (@$relational_db_response > 1){
 	    print STDERR "ERROR, checksum is not unique.\n";
	    $error = "Checksum is not unique.\n";
	    foreach my $row (@$relational_db_response){
		print STDERR join("\t",@$row ),"\n";
		$error .= join("\t",@$row )."\n";
	    }
	    return (undef, undef , $error);
 	}
    }

    return (undef, undef , $error);
}


=head3 md5_of_peg

    my $cksum = $fig->md5_of_peg( $peg );

Return the MD5 checksum for a peg. The MD5 checksum is computed from the
uppercase sequence of the protein.  This method retrieves the checksum stored
in the database.

=over 4

=item peg

FIG ID of the peg.

=item RETURN

Returns the checksum of the specified contig as a hex string, or C<undef> if
the peg is not in the database.

=back

=cut

#
#  Return md5sum of a peg as hex string:
#
#    my $md5sum = $fig->md5_of_peg( $peg )
#
sub md5_of_peg {
    my( $self, $peg ) = @_;
    return undef if ! $peg;

    #  Try to find it in the DBMS

    my $rdbH = $self->db_handle;
    my $dbms_response = $rdbH->SQL( "SELECT md5 FROM protein_sequence_MD5 WHERE id = '$peg'" );

    return $dbms_response->[0]->[0] if $dbms_response && @$dbms_response;

    #  Try to make it from the translation

    my $sequence = $self->get_translation( $peg );
    return undef if ( ! $sequence );

    #  Got a sequence, find the md5, save it in the DBMS, and return it

    my $md5 = Digest::MD5::md5_hex( uc $sequence );
    $rdbH->SQL( "INSERT INTO protein_sequence_MD5 ( id, md5 ) VALUES ( '$peg', '$md5' )" );

    return $md5;
}

sub md5_of_peg_bulk {
    my( $self, $pegs ) = @_;

    my @pegs_left = @$pegs;

    my %res;

    my $batch = 1000;

    #  Try to find it in the DBMS

    my $rdbH = $self->db_handle;
    my @to_compute;

    while (@pegs_left)
    {
	my @pegs = splice(@pegs_left, 0, $batch);
	my %pegs;
	$pegs{$_} = 1 foreach @pegs;

	my $in = join(", ", map { "?" } @pegs);
	my $dbms_response = $rdbH->SQL(qq(SELECT id, md5 FROM protein_sequence_MD5
					  WHERE id IN ($in)), undef, @pegs);

	for my $ent (@$dbms_response)
	{
	    my($peg, $md5) = @$ent;
	    $res{$peg} = $md5;
	    delete $pegs{$peg};
	}

	#
	# If there are any left, compute and insert.
	#
	push(@to_compute, keys %pegs);
    }
#    print  STDERR "Need to compute @to_compute\n";

    if (@to_compute)
    {
	my $sth = $rdbH->{_dbh}->prepare(qq(INSERT INTO protein_sequence_MD5 (id, md5)
					    VALUES (?, ?)));

	for my $peg (@to_compute)
	{
	    my $sequence = $self->get_translation( $peg );
	    if ($sequence)
	    {
		#  Got a sequence, find the md5, save it in the DBMS, and return it

		my $md5 = Digest::MD5::md5_hex( uc $sequence );
		$sth->execute($peg, $md5);
		$res{$peg} = $md5;
	    }
	}
    }

    return \%res;
}

=head3 get_representative_genome

	my $rep_id = get_representative_genome($id)

return the representative genome of the set that $id is in

=over 4

=item genome_id

ID of the genome used for set lookup

=item RETURN

Return the representative genome of the set that $id is in, 0 if not found

=back

=cut


sub get_representative_genome {
    my($self, $id) = @_;
    my $repH;

    if (! ($repH = $self->{_repG})) {
	my @tab = map { [split(/\t/,$_)] } `cat $FIG_Config::data/Global/genome.sets`;
	my $x = shift @tab;
	while ($x)
	{
	    my $set  = $x->[0];
	    my $repG = $x->[1];
	    while ($x && ($x->[0] eq $set))
	    {
		$repH->{$x->[1]} = $repG;
		$x = shift @tab;
	    }
	}
	$self->{_repG} = $repH;
    }
    return $repH->{$id};
}

=head3 pegs_with_md5

    my @pegs = $fig->pegs_with_md5( $md5 );

Return all pegs with sequence matching the check sum.  Thus,

    my @pegs = $fig->pegs_with_md5( $fig->md5_of_peg( $peg ) );

produces all pegs with sequence identical the query peg.

=over 4

=item md5

The md5 checksum as a hex string (32 characters).

=item RETURN

Returns the list of pegs matching the given md5 checksum.

=back

=cut

#
#  Return all peg ids with a given md5sum:
#
#    my @fids = $fig->pegs_with_md5( $md5_in_hex )
#
sub pegs_with_md5
{
    my( $self, $md5 ) = @_;

    return grep { /^fig\|/ } $self->prots_with_md5($md5);
}

=head3 prots_with_md5

    my @fids = $fig->prots_with_md5( $md5 );

Return all proteins with sequence matching the check sum, including
non fig ids.

=over 4

=item md5

The md5 checksum as a hex string (32 characters).

=item RETURN

Returns the list of protein ids matching the given md5 checksum.

=back

=cut

#
#  Return all protein ids with a given md5sum:
#
#    my @fids = $fig->prots_with_md5( $md5_in_hex )
#
sub prots_with_md5
{
    my( $self, $md5 ) = @_;
    return () if ! $md5;

    $md5 = lc $md5;

    #
    # Pull sequences from the protein_sequence_MD5 table.
    #

    my $rdbH = $self->db_handle;
    my $relational_db_response =
             $rdbH->SQL( "SELECT id FROM protein_sequence_MD5 WHERE md5 = '$md5'" );

    my %prots = map { $_->[0], 1 } @$relational_db_response;

    #
    # Add in the ones from the peg.synonyms.
    #

    map { $prots{$_->[0]} = 1 } $self->mapped_prot_ids("gnl|md5|$md5");

    return keys %prots;
}


=head3 genus_species

    my $gs = $fig->genus_species($genome_id);

Return the genus, species, and possibly also the strain of a specified genome.

This method converts a genome ID into a more recognizble species name. The species name
is stored directly in the genome table of the database. Essentially, if the strain is
present in the database, it will be returned by this method, and if it's not present,
it won't.

=over 4

=item genome_id

ID of the genome whose name is desired.

=item RETURN

Returns the scientific species name associated with the specified ID, or C<undef> if the
ID is not in the database.

=back

=cut
#: Return Type $;
sub genus_species :Scalar {
    my ($self,$genome) = @_;
    my $ans;

    # if the genome is marked deleted and you are trying to undelete it, this will not be set
    if ($self->{gdata}->gname($genome)) {return $self->{gdata}->gname($genome)}

    my $genus_species = $self->cached('_genus_species');
    unless (scalar(keys(%$genus_species))) {
        my $rdbH = $self->db_handle;
        my $relational_db_response = $rdbH->SQL("SELECT genome,gname  FROM genome");
        my $pair;
        foreach $pair (@$relational_db_response) {
            $genus_species->{$pair->[0]} = $pair->[1];
        }
    }

    $ans = $genus_species->{$genome};
    if ((! $ans) && open(GEN,"<$FIG_Config::organisms/$genome/GENOME"))
    {
	$ans = <GEN>;
	close(GEN);
	chomp $ans;
	$genus_species->{$genome} = $ans;
    }

    if ($ans)
    {
	$ans =~ s/^\s+//o;
	$ans =~ s/^Candidatus\s+//o;
	$ans =~ s/\s+/ /go;
	$ans =~ s/\s+$//o;
    }

    return $ans;
}

=head3 set_genus_species

    my $gs = $fig->set_genus_species($genome_id, $genus_species_strain);

Sets the contents of the GENOME file of the specified genome ID

Does not (currently) update the relational DB.
UPDATE:
edited by RAE: this now sets the gname column in the 'genome' table.

=over 4

=item genome_id

ID of the genome whose name is desired.

=item genus_species_strain

The new biological name that will correspond to the  genome_id.

=item RETURN

Returns C<1> if the write was successful, and C<undef> if write fails.

=back

=cut
#: Return Type $;
sub set_genus_species :Scalar {
    my ($self, $genome, $genus_species_strain) = @_;
    chomp $genus_species_strain;

    my $genome_file = "$FIG_Config::organisms/$genome/GENOME";

    if (!-f $genome_file) {
	warn "$genome_file doe not exist";
	return undef;
    }
    else {
	if (system("cp -p $genome_file $genome_file~")) {
	    warn "Could not back up $genome_file";
	    return undef;
	}
	else {
	    if (not open(GENOME, ">$genome_file")) {
		warn "Could not write-open $genome_file";
		return undef;
	    }
	    else {
		print GENOME "$genus_species_strain\n";
		close(GENOME) || warn "Could not close genome file $genome_file";

		# RAE: Added an update to the RDB
		my $rdbH = $self->db_handle;

		my $relational_db_response = $rdbH->SQL("UPDATE genome SET gname= ? WHERE genome = ?", undef, $genus_species_strain, $genome);
		return 1;
	    }
	}
    }
}


=head3 org_of

    my $org = $fig->org_of($prot_id);

Return the genus/species name of the organism containing a protein. Note that in this context
I<protein> is not a certain string of amino acids but a protein encoding region on a specific
contig.

For a FIG protein ID (e.g. C<fig|134537.1.peg.123>), the organism and strain
information is always available. In the case of external proteins, we can usually
determine an organism, but not anything more precise than genus/species (and
often not that). When the organism name is not present, a null string is returned.

=over 4

=item prot_id

Protein or feature ID.

=item RETURN

Returns the displayable scientific name (genus, species, and strain) of the organism containing
the identified PEG. If the name is not available, returns a null string. If the PEG is not found,
returns C<undef>.

=back

=cut

sub org_of {
    my($self,$prot_id) = @_;
    my $relational_db_response;
    my $rdbH = $self->db_handle;

    if ($prot_id =~ /^fig\|/) {
        return  $self->is_deleted_fid( $prot_id) ? undef
            : $self->genus_species( $self->genome_of( $prot_id ) ) || "";
    }

    if (($relational_db_response =
            $rdbH->SQL("SELECT org FROM external_orgs WHERE ( prot = \'$prot_id\' )")) &&
            (@$relational_db_response >= 1)) {
        $relational_db_response->[0]->[0] =~ s/^\d+://;
        return $relational_db_response->[0]->[0];
    }
    return "";
}

=head3 orgid_of_orgname

    my $genomeID = $fig->orgid_of_orgname($genomeName);

Return the ID of the genome corresponding to the specified organism name, or a
null string if the genome is not found.

=over 4

=item genomeName

Name of the organism, consisting of the organism's genus, species, and
unique characterization, separated by spaces.

=item RETURN

Returns the genome ID number for the named organism, or an empty string if
the genome is not found.

=back

=cut

sub orgid_of_orgname {
  my($self,$genome_name) = @_;
  my $relational_db_response;
  my $rdbH = $self->db_handle;

  my $genome_nameQ = quotemeta $genome_name;

  if (($relational_db_response =
       $rdbH->SQL("SELECT genome FROM genome WHERE gname='$genome_nameQ'")) &&
      (@$relational_db_response >= 1)) {
    return $relational_db_response->[0]->[0];
  }
  return "";
}

=head3 orgname_of_orgid

    my $genomeName = $fig->orgname_of_orgid($genomeID);

Return the name of the genome corresponding to the specified organism ID.

=over 4

=item genomeID

ID of the relevant genome.

=item RETURN

Returns the name of the organism, consisting of the organism's genus, species, and
unique characterization, separated by spaces, or a null string if the genome is not
found.

=back

=cut

sub orgname_of_orgid {
  my($self,$genome) = @_;
  my $relational_db_response;
  my $rdbH = $self->db_handle;

  if (($relational_db_response =
       $rdbH->SQL("SELECT gname FROM genome WHERE genome='$genome'")) &&
      (@$relational_db_response >= 1)) {
    return $relational_db_response->[0]->[0];
  }
  return "";
}

=head3 genus_species_domain

    my ($gs, $domain) = $fig->genus_species_domain($genome_id);

Returns a genome's genus and species (and strain if that has been properly
recorded) in a printable form, along with its domain. This method is similar
to L</genus_species>, except it also returns the domain name (archaea,
bacteria, etc.).

=over 4

=item genome_id

ID of the genome whose species and domain information is desired.

=item RETURN

Returns a two-element list. The first element is the species name and the
second is the domain name.

=back

=cut

sub genus_species_domain {
    my ($self, $genome) = @_;

    my $genus_species_domain = $self->cached('_genus_species_domain');
    if ( ! $genus_species_domain->{ $genome } ) {
        my $rdbH = $self->db_handle;
        my $relational_db_response = $rdbH->SQL("SELECT genome,gname,maindomain FROM genome");
        my $triple;
        foreach $triple ( @$relational_db_response ) {
            $genus_species_domain->{ $triple->[0] } = [ $triple->[1], $triple->[2] ];
        }
    }
    my $gsdref = $genus_species_domain->{ $genome };
    return $gsdref ? @$gsdref : ( "", "" );
}

=head3 domain_color

    my $web_color = FIG::domain_color($domain);

Return the web color string associated with a specified domain. The colors are
extremely subtle (86% luminance), so they absolutely require a black background.
Archaea are slightly cyan, bacteria are slightly magenta, eukaryota are slightly
yellow, viruses are slightly silver, environmental samples are slightly gray,
and unknown or invalid domains are pure white.

=over 4

=item domain

Name of the domain whose color is desired.

=item RETURN

Returns a web color string for the specified domain (e.g. C<#FFDDFF> for
bacteria).

=back

=cut

my %domain_color = ( AR => "#DDFFFF", BA => "#FFDDFF", EU => "#FFFFDD",
                     VI => "#DDDDDD", EN => "#BBBBBB" );

sub domain_color {
    my ($domain) = @_;
    defined $domain || return "#FFFFFF";
    return $domain_color{ uc substr($domain, 0, 2) } || "#FFFFFF";
}

=head3 org_and_color_of

    my ($org, $color) = $fig->org_and_domain_of($prot_id);

Return the best guess organism and domain html color string of an organism.
In the case of external proteins, we can usually determine an organism, but not
anything more precise than genus/species (and often not that).

=over 4

=item prot_id

Relevant protein or feature ID.

=item RETURN

Returns a two-element list. The first element is the displayable organism name, and the second
is an HTML color string based on the domain (see L</domain_color>).

=back

=cut

sub org_and_color_of {
    my($self,$prot_id) = @_;
    my $relational_db_response;
    my $rdbH = $self->db_handle;

    if ($prot_id =~ /^fig\|/) {
        my( $gs, $domain ) = $self->genus_species_domain($self->genome_of($prot_id));
        return ( $gs, domain_color( $domain ) );
    }

    if (($relational_db_response =
         $rdbH->SQL("SELECT org FROM external_orgs WHERE ( prot = \'$prot_id\' )")) &&
         (@$relational_db_response >= 1)) {
        return ($relational_db_response->[0]->[0], "#FFFFFF");
    }
    return ("", "#FFFFFF");
}

=head3 partial_genus_matching

Return a list of genome IDs that match a partial genus.

For example partial_genus_matching("Listeria") will return all genome IDs that begin with Listeria, and this can also be restricted to complete genomes with another argument like this partial_genus_matching("Listeria", 1)

=cut

sub partial_genus_matching {
 my ($self, $gen, $complete)=@_;
 return grep {$self->genus_species($_) =~ /$gen/i} $self->genomes($complete);
}


=head3 abbrev

    my $abbreviated_name = FIG::abbrev($genome_name);

or

    my $abbreviated_name = $fig->abbrev($genome_name);

Abbreviate a genome name to 10 characters or less.

For alignments and such, it is very useful to be able to produce an abbreviation of genus/species.
That's what this does.  Note that multiple genus/species might reduce to the same abbreviation, so
be careful (disambiguate them, if you must).

The abbreviation is formed from the first three letters of the species name followed by the
first three letters of the genus name followed by the first three letters of the species name and
then the next four nonblank characters.

=over 4

=item genome_name

The name to abbreviate.

=item RETURN

An abbreviated version of the specified name.

=back

=cut

sub abbrev :Scalar {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($genome_name) = @_;

    $genome_name =~ s/^(\S{3})\S+/$1./;
    $genome_name =~ s/^(\S+)\s+(\S{3})\S+/$1$2./;
    $genome_name =~ s/ //g;
    if (length($genome_name) > 10) {
        $genome_name = substr($genome_name,0,10);
    }
    return $genome_name;
}

=head3 wikipedia_link

    my $wikipedia_link = $fig->wikipedia_link($genome_name);

Check if Wikipedia has a page about this genome. If so, return it's url.

=over 4

=item genome_name

The genome to find.

=item RETURN

The url of the wikipedia page.

=back

=cut

sub wikipedia_link {
    my ($self, $organism_name) = @_;

    if ($self->db_handle->table_exists('genome_wikipedia_link'))
    {
	my @organism_tokens = split(/\s/, $organism_name);
	my $gs = join(" ", @organism_tokens[0..1]);

	my $res = $self->db_handle->SQL(qq(SELECT url FROM genome_wikipedia_link WHERE gname = ?), undef,
					$gs);
	if (@$res)
	{
	    return $res->[0]->[0];
	}
	else
	{
	    return undef;
	}
    }
    else
    {
	return FIGRules::wikipedia_link($organism_name);
    }
}

=head3 organism_directory

    my $organism_directory = $fig->organism_directory($genome_id);

Get the directory that contains the organism data. This is just like the
FIGV version.

=over 4

=item genome_id

The id of the organism, e.g. 83333.1.

=item RETURN

A string containing the path to the organism directory.

=back

=cut

sub organism_directory {
    my ($self, $org_id) = @_;
    my $org_dir = "$FIG_Config::organisms/$org_id";
    return ((-d $org_dir) ? $org_dir : undef);
}

=head3 ncbi_contig_description

C<<my $name = ncbi_contig_description($contig_id)>>

Looks up the NCBI description line for this contig identifier. Values are cached
in the directory $FIG_Config::var/ncbi_contigs.

=cut

sub ncbi_contig_description
{
    my($self, $id) = @_;

    my $cache_dir = "$FIG_Config::fig/var/ncbi_contigs";
    &FIG::verify_dir($cache_dir);

    my $cache_file = "$cache_dir/$id";
    if (open(CF, $cache_file))
    {
        $_ = <CF>;
        close(CF);
        chomp;
        if ($_ ne '')
        {
            return $_;
        }
    }

    my $last_lookup = $self->{_ncbi_last_lookup};
    if ($last_lookup =~ /\d+/)
    {
        my $wait = $last_lookup + 3 - time;
        if ($wait > 0)
        {
            warn "waiting $wait for lookup\n";
            sleep($wait);
        }
    }
    $self->{_ncbi_last_lookup} = time;

    my $ua = new LWP::UserAgent;
    my $utils = "http://www.ncbi.nlm.nih.gov/entrez/eutils";

    my %params = (db => 'genome',
                  usehistory => 'y',
                  term => $id);
    my $res = $ua->get("$utils/esearch.fcgi?" . join("&", map { "$_=$params{$_}" } keys %params));
    if (not $res->is_success)
    {
        warn "esearch failed: " . $res->content;
        return;
    }
    %params = (db => 'genome',
                  usehistory => 'y',
                  term => $id);

    my $esearch_result = $res->content;
    $esearch_result =~
        m|<Count>(\d+)</Count>.*<QueryKey>(\d+)</QueryKey>.*<WebEnv>(\S+)</WebEnv>|s;

    my $Count    = $1;
    my $QueryKey = $2;
    my $WebEnv   = $3;

    %params = (rettype => 'summary',
               retmode => 'text',
               db => 'genome',
               query_key => $QueryKey,
               WebEnv => $WebEnv);

    $res = $ua->get("$utils/efetch.fcgi?" . join("&", map { "$_=$params{$_}" } keys %params));
    if (not $res->is_success)
    {
        warn "esearch failed: " . $res->content;
        return;
    }

    my $txt = $res->content;
    my($start, $ident);
    while ($txt =~ /([^\n]*)\n/sg)
    {
        my $l = $1;
        if ($l =~ /^\d+:\s+/)
        {
            $start = 1;
        }
        elsif ($start)
        {
            $ident = $l;
            last;
        }
    }
    print "Got ident $ident\n";
    if (open(CF, ">$cache_file"))
    {
        print CF "$ident\n";
        close(CF);
    }
    return $ident;
}


################ Routines to process Features and Feature IDs  ##########################

=head3 ftype

    my $type = FIG::ftype($fid);

or

    my $type = $fig->ftype($fid);

Returns the type of a feature, given the feature ID.  This just amounts
to lifting it out of the feature ID, since features have IDs of the form

        fig|x.y.f.n

where
        x.y is the genome ID
        f   is the type of feature
        n   is an integer that is unique within the genome/type

=over 4

=item fid

FIG ID of the feature whose type is desired.

=item RETURN

Returns the feature type (e.g. C<peg>, C<rna>, C<pi>, or C<pp>), or C<undef> if the
feature ID is not a FIG ID.

=back

=cut

sub ftype {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my $feature_id = $_[0] || '';

    return ($feature_id =~ /^fig\|\d+\.\d+\.([^\.]+)/) ? $1 : undef;
}

=head3 genome_of

    my $genome_id = $fig->genome_of($fid);

or

    my $genome_id = FIG::genome_of($fid);

Return the genome ID from a feature ID.

=over 4

=item fid

ID of the feature whose genome ID is desired.

=item RETURN

If the feature ID is a FIG ID, returns the genome ID embedded inside it; otherwise, it
returns C<undef>.

=back

=cut


sub genome_of :Scalar {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    # my $feature_id = (@_ == 1) ? $_[0] : $_[1];
    my $feature_id = $_[0] || '';

    return ($feature_id =~ /^fig\|(\d+\.\d+)/) ? $1 : undef;
}

=head3 genome_and_peg_of

    my ($genome_id, $peg_number = FIG::genome_and_peg_of($fid);

    my ($genome_id, $peg_number = $fig->genome_and_peg_of($fid);

Return the genome ID and peg number from a feature ID.

=over 4

=item prot_id

ID of the feature whose genome and PEG number as desired.

=item RETURN

Returns the genome ID and peg number associated with a feature if the feature
is represented by a FIG ID, else C<undef>.

=back

=cut

sub genome_and_peg_of {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my $prot_id = (@_ == 1) ? $_[0] : $_[1];

    if ($prot_id =~ /^fig\|(\d+\.\d+)\.peg\.(\d+)/) {
        return ($1, $2);
    }
    return undef;
}



=head3 by_fig_id

    my @sorted_by_fig_id = sort { FIG::by_fig_id($a,$b) } @fig_ids;

Compare two feature IDs.

This function is designed to assist in sorting features by ID. The sort is by
genome ID followed by feature type and then feature number.

=over 4

=item a

First feature ID.

=item b

Second feature ID.

=item RETURN

Returns a negative number if the first parameter is smaller, zero if both parameters
are equal, and a positive number if the first parameter is greater.

=back

=cut

sub by_fig_id {
    my($a,$b) = @_;
    my($g1,$g2,$t1,$t2,$n1,$n2);
    if (($a =~ /^fig\|(\d+\.\d+).([^\.]+)\.(\d+)$/) && (($g1,$t1,$n1) = ($1,$2,$3)) &&
         ($b =~ /^fig\|(\d+\.\d+).([^\.]+)\.(\d+)$/) && (($g2,$t2,$n2) = ($1,$2,$3))) {
        ($g1 <=> $g2) or ($t1 cmp $t2) or ($n1 <=> $n2);
    } else {
        $a cmp $b;
    }
}



=head3 by_fig_id

C<< my @sorted_by_location = sort { FIG::by_locus($a,$b) } @locations; >>

Compare two locations.

This function is designed to assist in sorting features by location.
The sort is by contig ID, followed by left boundary, then by right bounday,
then by strand.

=over 4

=item a

First location.

=item b

Second location.

=item RETURN

Returns a negative number if the first location is to the left, of the second,
zero if both locations are identical, and a positive number if the first location
is to the right of the second.

=back

=cut

sub by_locus
{
    my ($a, $b) = @_;

    my ($A_contig, $A_beg, $A_end) = &boundaries_of($a);
    my ($B_contig, $B_beg, $B_end) = &boundaries_of($b);

    return (  ($A_contig cmp $B_contig)
	   || (&min($A_beg, $A_end) <=> &min($B_beg, $B_end))
	   || (&max($A_beg, $A_end) <=> &max($B_beg, $B_end))
	   || ( ($A_beg <=> $A_end) <=> ($B_beg <=> $B_end) )
	   );
}



=head3 by_genome_id

    my @sorted_by_genome_id = sort { FIG::by_genome_id($a,$b) } @genome_ids;

Compare two genome IDs.

This function is designed to assist in sorting genomes by ID.

=over 4

=item a

First genome ID.

=item b

Second genome ID.

=item RETURN

Returns a negative number if the first parameter is smaller, zero if both parameters
are equal, and a positive number if the first parameter is greater.

=back

=cut

sub by_genome_id {
    my($a,$b) = @_;
    my($g1,$g2,$s1, $s2);
    if (($a =~ /^(\d+)\.(\d+)$/) && (($g1, $s1) = ($1, $2)) &&
        ($b =~ /^(\d+)\.(\d+)$/) && (($g2, $s2) = ($1, $2))) {
        ($g1 <=> $g2) or ($s1 <=> $s2);
    } else {
        $a cmp $b;
    }
}


=head3 next_feature

    my $feature = $fig->next_feature( \%options );

Locate the next feature (optionally filtered by type) in a contig. The start
position for the search can be defined by supplying genome, contig and
position, or by supplying a feature id.  Feature locations are defined by
their midpoint.  If a fid is supplied with contig and position, the latter
are used to resolve ambiguities in the desired segement of a feature with
a complex location.

=over 4

=item options

Options:

after =>  $fid
after => \@fids

Id(s) of features that should preceed the returned feature.  This is a local
operation, and is only meant to resolve features that are otherwise tied in
location.

contig => $contig

Name of contig of features.

exclude =>  $id
exclude => \@ids

Id(s) of features to exclude.  Note that features listed with the 'after'
option are also excluded (and that is most commonly the desired behavior).

fid => $fid

Alternative to supplying a location.  It is possible to supply a fid and
contig and position, which allows disambiguating the desired segment of
a feature with a complex location.

genome => $genome

Name of genome of features.

position => $position

Feature midpoint must be >= $position.  Note that this can be any multiple
of 1/2.  If the supplied value is negative, the position is taken from the
right end of the contig.

type =>  $type
type => \@types

Type(s) of desired feature (default is any type).

=item RETURN

Feature id or undef.

=back

=cut

sub next_feature
{
    my ( $self, $options ) = @_;
    return undef unless ref( $options ) eq 'HASH';

    my $fid      = $options->{ fid };
    my $genome   = $options->{ genome };
    my $contig   = $options->{ contig };
    my $position = $options->{ position };

    if ( ! $genome )
    {
      return undef unless $fid;
      ( $genome ) = $fid =~ /^fig\|(\d+\.\d+)\./;
    }

    if ( ! $contig || ! $position )
    {
      return undef unless $fid;
      my ( $region ) = grep { ( ! $contig ) || ( $_->[0] eq $contig ) }
                       FIG::boundaries_of_2( scalar $self->feature_location( $fid ) );
      my ( $beg, $end );
      ( $contig, $beg, $end ) = @{ $region || [] };
      return undef unless ( $contig && $beg && $end );
      $position = 0.5 * ( $beg + $end );
    }

    my $length = $self->contig_ln( $genome, $contig );
    return undef unless $length;

    #  Negative position counts from left end of contig

    if ( $position < 0 ) { $position += $length + 1 }
    my $pos2 = 2 * $position;

    my $after = $options->{ after };
    my @after = ref( $after ) eq 'ARRAY' ?  @$after
              : $after                   ? ( $after )
              : ();
    my %after = map { $_ => 1 } @after;
    $after{ $fid } = 1 if $fid;

    my $exclude = $options->{ exclude };
    my @exclude = ref( $exclude ) eq 'ARRAY' ?  @$exclude
                : $exclude                   ? ( $exclude )
                : ();
    my %exclude = map { $after{ $_ } ? () : ( $_ => 1 ) } @exclude;

    my $type  = $options->{ type };
    my @types = ref( $type ) eq 'ARRAY' ?  @$type
              : $type                   ? ( $type )
              : ();
    my %types = map { $_ => 1 } @types;

    my $rdbH = $self->db_handle;
    my $minV = int( $position + 0.5 );  #  Round up if not integer
    my $maxV = $minV + 9999;
    my $relational_db_response;

    while ( $minV <= $length )
    {
        my $query = "SELECT id "
                  . "FROM features "
                  . "WHERE ( minloc <= $maxV ) "
                    . "AND ( maxloc >= $minV ) "
                    . "AND ( genome = \'$genome\' ) "
                    . "AND ( contig = \'$contig\' );";
        ( $relational_db_response = $rdbH->SQL( $query ) ) or return undef;

        if ( @$relational_db_response >= 1 )
        {
            my @tmp = map  { my $loc  = $self->feature_location( $_ );
                             my @loc2 = filtered_location( $loc, $contig, $position, 0.5*$length );
                             my ( $type, $numb ) = $_ =~ /\.([^.]+)\.(\d+)$/;
                             $loc2[0] && ( $loc2[1]+$loc2[2] >= $pos2 ) ? [ $_, @loc2, $type, $numb+0 ] : ()
                           }
                      grep { my ( $type ) = $_ =~ /^fig\|\d+\.\d+\.([^.]+)\.\d+$/;
                             $type && !$exclude{$_} && ( !@types || $types{$type} )
                           }
                      map  { $_->[0] }           # extract the id
                      @$relational_db_response;
	    # @tmp contains [$id, $contig, $beg, $end, $type, $numb]
            if ( @tmp )
            {
                #  This sort produces an unambiguous ordering of all features.
                #  In a more general ordering, the contig would also be sorted,
                #  but presently we know that they are all on one contig.
                @tmp = sort {    ($a->[2]+$a->[3]) <=>    ($b->[2]+$b->[3])  # midpoint
                           || min($a->[2],$a->[3]) <=> min($b->[2],$b->[3])  # left end
                           ||  lc $a->[4]          cmp  lc $b->[4]           # type
                           ||     $a->[5]          <=>     $b->[5]           # number
                            }
                       @tmp;

                #  Process the 'after' ids:
                my $data;
                foreach ( reverse @tmp )
                {
                    last if ( $after{ $_->[0] } );
                    $data = $_;
                }

                #  An id was found.  Make sure that its midpoint was within
                #  the range of coordinates surveyed.  If not, there might be
                #  a smaller feature with a closer midpoint, but which was not
                #  found in the DB query.

                if ( $data )
                {
                    my $mid  = 0.5 * ( $data->[2] + $data->[3] );

                    #======== This is the only return with a defined id ========

                    return $data->[0] if $mid <= $maxV;

                    $maxV = int( $mid + 0.5 );  #  Might have missed a shorter feature
                }
                else  # Nothing yet; look further.
                {
                    $minV  = $maxV + 1;
                    $maxV += 10000;
                }
            }
            else  # Nothing yet; look further.
            {
                $minV  = $maxV + 1;
                $maxV += 10000;
            }
        }
        else  # Nothing yet; look further.
        {
            $minV  = $maxV + 1;
            $maxV += 10000;
        }
    }

    return undef;
}


=head3 previous_feature

    my $feature = $fig->previous_feature( \%options );

Locate the previous feature (optionally filtered by type) in a contig. The
start position for the search can be defined by supplying genome, contig and
position, or by supplying a feature id.  Feature locations are defined by
their midpoint.  If a fid is supplied with contig and position, the latter
are used to resolve ambiguities in the desired segement of a feature with
a complex location.

=over 4

=item options

Options:

before =>  $fid
before => \@fids

Id(s) of features that should follow the returned feature.  This is a local
operation, and is only meant to resolve features that are otherwise tied in
location.

contig => $contig

Name of contig of features.

exclude =>  $id
exclude => \@ids

Id(s) of features to exclude.  Note that features listed with the 'before'
option are also excluded (and that is most commonly the desired behavior).

fid => $fid

Alternative to supplying a location.  It is possible to supply a fid and
contig and position, which allows disambiguating the desired segment of
a feature with a complex location.

genome => $genome

Name of genome of features.

position => $position

Feature midpoint must be >= $position.  Note that this can be any multiple
of 1/2.  If the supplied value is negative, the position is taken from the
right end of the contig.

type =>  $type
type => \@types

Type(s) of desired feature (default is any type).

=item RETURN

Feature id or undef.

=back

=cut

sub previous_feature
{
    my ( $self, $options ) = @_;
    return undef unless ref( $options ) eq 'HASH';

    my $fid      = $options->{ fid };
    my $genome   = $options->{ genome };
    my $contig   = $options->{ contig };
    my $position = $options->{ position };

    if ( ! $genome )
    {
      return undef unless $fid;
      ( $genome ) = $fid =~ /^fig\|(\d+\.\d+)\./;
    }

    if ( ! $contig || ! $position )
    {
      return undef unless $fid;
      my ( $region ) = grep { ( ! $contig ) || ( $_->[0] eq $contig ) }
                       FIG::boundaries_of_2( scalar $self->feature_location( $fid ) );
      my ( $beg, $end );
      ( $contig, $beg, $end ) = @{ $region || [] };
      return undef unless ( $contig && $beg && $end );
      $position = 0.5 * ( $beg + $end );
    }

    my $length = $self->contig_ln( $genome, $contig );
    return undef unless $length;

    #  Negative position counts from left end of contig

    if ( $position < 0 ) { $position += $length + 1 }
    my $pos2 = 2 * $position;

    my $before = $options->{ before };
    my @before = ref( $before ) eq 'ARRAY' ?  @$before
               : $before                   ? ( $before )
               : ();
    my %before = map { $_ => 1 } @before;
    $before{ $fid } = 1 if $fid;

    my $exclude = $options->{ exclude };
    my @exclude = ref( $exclude ) eq 'ARRAY' ?  @$exclude
                : $exclude                   ? ( $exclude )
                : ();
    my %exclude = map { $before{ $_ } ? () : ( $_ => 1 ) } @exclude;

    my $type  = $options->{ type };
    my @types = ref( $type ) eq 'ARRAY' ?  @$type
              : $type                   ? ( $type )
              : ();
    my %types = map { $_ => 1 } @types;

    my $rdbH = $self->db_handle;
    my $maxV = int( $position );   #  Round down if not integer
    my $minV = $maxV - 9999;
    my $relational_db_response;

    while ( $maxV >= 1 )
    {
        my $query = "SELECT id "
                  . "FROM features "
                  . "WHERE ( minloc <= $maxV ) "
                    . "AND ( maxloc >= $minV ) "
                    . "AND ( genome = \'$genome\' ) "
                    . "AND ( contig = \'$contig\' );";
        ( $relational_db_response = $rdbH->SQL( $query ) ) or return undef;

        if ( @$relational_db_response >= 1 )
        {
            my @tmp = map  { my $loc  = $self->feature_location( $_ );
                             my @loc2 = filtered_location( $loc, $contig, $position, 0.5*$length );
                             my ( $type, $numb ) = $_ =~ /\.([^.]+)\.(\d+)$/;
                             ( $loc2[0] && ( $loc2[1]+$loc2[2] <= $pos2 ) ) ? [ $_, @loc2, $type, $numb+0 ] : ()
                           }
                      grep { my ( $type ) = $_ =~ /^fig\|\d+\.\d+\.([^.]+)\.\d+$/;
                             $type && !$exclude{$_} && ( !@types || $types{$type} )
                           }
                      map  { $_->[0] }           # extract the id
                      @$relational_db_response;

            if ( @tmp )
            {
                #  This sort produces an unambiguous ordering of all features.
                #  In a more general ordering, the contig would also be sorted,
                #  but presently we know that they are all on one contig.
                @tmp = sort {    ($a->[2]+$a->[3]) <=>    ($b->[2]+$b->[3])  # midpoint
                           || min($a->[2],$a->[3]) <=> min($b->[2],$b->[3])  # left end
                           ||  lc $a->[4]          cmp  lc $b->[4]           # type
                           ||     $a->[5]          <=>     $b->[5]           # number
                            }
                       @tmp;
                # if ( 1 ) { print STDERR Dumper( \@tmp ) }

                #  Process the 'before' ids:
                my $data;
                foreach ( @tmp )
                {
                    last if ( $before{ $_->[0] } );
                    $data = $_;
                }

                #  An id was found.  Make sure that its midpoint was within
                #  the range of coordinates surveyed.  If not, there might be
                #  a smaller feature with a closer midpoint, but which was not
                #  found in the DB query.

                if ( $data )
                {
                    my $mid  = 0.5 * ( $data->[2] + $data->[3] );

                    #======== This is the only return with a defined id ========

                    return $data->[0] if $mid >= $minV;

                    $minV = int( $mid );  #  Might have missed a shorter feature
                }
                else  # Nothing yet; look further.
                {
                    $maxV  = $minV - 1;
                    $minV -= 10000;
                }
            }
            else  # Nothing yet; look further.
            {
                $maxV  = $minV - 1;
                $minV -= 10000;
            }
        }
        else  # Nothing yet; look further.
        {
            $maxV  = $minV - 1;
            $minV -= 10000;
        }
    }

    return undef;
}


#  Process a (possibly complex) location for the most applicable segment(s).
#  Mostly this is for dealing locations that wrap through the origin of a
#  contig.
#
#     ( $contig, $beg, $end ) = filtered_location( $loc, $contig, $focus, $range )
#
sub filtered_location
{
    my ( $loc, $contig, $focus, $range ) = @_;
    my @regions = grep { $_->[0] eq $contig } boundaries_of_2( $loc );
    if ( @regions < 2 ) { return @{ $regions[0] || [] } }

    #  Okay, life was not simple.  Let's see if we can throw out the most
    #  distant parts.

    my $min_mid = $focus - $range;
    my $max_mid = $focus + $range;
    @regions = grep { my $mid = 0.5 * ( $_->[1] + $_->[2] );
                      $mid >= $min_mid && $mid <= $max_mid
                    }
               @regions;
    if ( @regions < 2 ) { return @{ $regions[0] || [] } }

    #  This should be very rare, and returning an empty list would be
    #  reasonable.  However, let's fall back to an interval that covers
    #  everything remaining.

    my ( $beg, $end ) = ( shift @regions )[1,2];
    foreach ( @regions )
    {
        my ( $b1, $e1 ) = @$_[1,2];
        if ( $beg < $end ) { $beg = $b1 if $b1 < $beg; $end = $e1 if $e1 > $end }
        else               { $beg = $b1 if $b1 > $beg; $end = $e1 if $e1 < $end }
    }

    ( $contig, $beg, $end )
}


=head3 genes_in_region

    my ($features_in_region, $beg1, $end1) = $fig->genes_in_region($genome, $contig, $beg, $end, size_limit);

Locate features that overlap a specified region of a contig. This includes features that begin or end
outside that region, just so long as some part of the feature can be found in the region of interest.

It is often important to be able to find the genes that occur in a specific region on
a chromosome.  This routine is designed to provide this information.  It returns all genes
that overlap positions from I<$beg> through I<$end> in the specified contig.

The I<$size_limit> parameter limits the search process. It is presumed that no features are longer than the
specified size limit. A shorter size limit means you'll miss some features; a longer size limit significantly
slows the search process. For prokaryotes, a value of C<10,000> (the default) seems to work best. With some
optimization of the code, I see minimal degradation with C<100,000> (the new default), and this will pick-up
essentially all known prophage, pathogenicity islands and polyketide synthetases.

=over 4

=item genome

ID of the genome containing the relevant contig.

=item contig

ID of the relevant contig.

=item beg

Position of the first base pair in the region of interest.

=item end

Position of the last base pair in the region of interest.

=item size_limit

Maximum allowable size for a feature. If omitted, C<100,000> is assumed. Note that this is not actually
the maximum size, but no gene smaller than this can be missed.

=item RETURN

Returns a three-element list. The first element is a reference to a list of the feature IDs found. The second
element is the position of the leftmost base pair in any feature found. This may be well before the region of
interest begins or it could be somewhere inside. The third element is the position of the rightmost base pair
in any feature found. Again, this can be somewhere inside the region or it could be well to the right of it.

=back

=cut
#: Return Type @;
sub genes_in_region
{
    # require Time::HiRes;
    # my $t1 = Time::HiRes::time();
    my($self, $genome, $contig, $beg, $end, $pad, $allowed_types) = @_;
    defined( $self ) && $genome && defined( $contig ) && defined( $beg ) && defined( $end )
        or return ([], undef, undef);

    my $rdbH = $self->db_handle
        or return ([], undef, undef);

    $pad = 100000 if ! $pad;
    my $minV = $beg - $pad;
    my $tq = "";
    my @qs = ();
    if (ref($allowed_types) eq 'ARRAY')
    {
	$tq = "AND (type IN (" . join(", ", map { "?" } @$allowed_types) . ")) ";
    }
    else
    {
	$allowed_types = [];
    }
    my $relational_db_response = $rdbH->SQL( "SELECT id,location FROM features "
                                           . " WHERE ( genome =  ? ) "
                                           . "   AND ( contig =  ? ) "
					   . $tq
                                           . "   AND ( minloc <= ? ) AND ( minloc > ? ) "
                                           . "   AND ( maxloc >= ? );",
                                             undef, $genome, $contig, @$allowed_types, $end, $minV, $beg
                                           ) || [];
    @$relational_db_response
        or return ( [], undef, undef );

    my @kept;
    my ( $l, $u ) = ( 10000000000, 0 );
    foreach ( @$relational_db_response )
    {
        my ( $fid, $loc ) = @$_;
        my @segs = grep { $_->[3] <= $end && $_->[4] >= $beg }
                   location_segments( $loc, $contig );
        @segs or next;

        my $keep = [ $fid, ( @{shift @segs} )[3,4], 0 ];
        foreach ( @segs )
        {
            $keep->[1] = $_->[3] if $_->[3] < $_->[1];
            $keep->[2] = $_->[4] if $_->[4] > $_->[2];
        }

        $keep->[3] = $keep->[1] + $keep->[2];
        $l = $keep->[1] if $keep->[1] < $l;
        $u = $keep->[2] if $keep->[2] > $u;

        push @kept, $keep;
    }

    my @fids = map  { $_->[0] }                #  Keep just the fid
               sort { $a->[3] <=> $b->[3] }    #  Sort by midpoints
               @kept;

    my $del = $self->is_deleted_fid_bulk(@fids);
    @fids = grep { ! $del->{$_} } @fids;

    # my $t2 = Time::HiRes::time();
    # printf STDERR "Elapse = %.4f sec\n", $t2-$t1;

    @fids ?  ( \@fids, $l, $u ) : ( [], undef, undef );
}


#=============================================================================
#  To make this work as desired, we really need to consider the consistently
#  oriented and ordered segments of feature locations. We can then ask whether
#  any part of a segment properly overlaps the region in question. If a contig
#  is supplied, the segments are filtered to those on the given contig.
#
#     @segments = $fig->location_segments( $location )
#     @segments = $fig->location_segments( $location, $contig )
#
#     Each segment is [ $contig, $beg, $end, $min, $max ]
#
#=============================================================================
sub location_segments
{
    my ( $loc, $contig ) = @_;
    $loc or return ();

    #  Location parts are [ $contig, $beg, $end, $dir, 2*$midpoint ]
    #  No parts of length 1 are allowed (the direction is ambiguous).
    my @parts = map { /^(.+)_(\d+)_(\d+)$/ && $2 != $3 ? [ $1, $2, $3, $3<=>$2, $2+$3 ] : () }
                split /,/, $loc;
    @parts or return ();

    my $seg = shift @parts;
    my @segs = ( $seg );
    foreach my $part ( @parts )
    {
        if ( ( $seg->[0] eq $part->[0] )                   # Same contig
          && ( $seg->[3] == $part->[3] )                   # Same direction
          && ( $seg->[3] == ( $part->[4] <=> $seg->[4] ) ) # Movement in correct dir
           )
        {
            $seg->[2] = $part->[2];           # Move endpoint
            $seg->[4] = $seg->[1]+$seg->[2];  # Move midpoint
        }
        else
        {
            push @segs, ($seg = $part); # Make this a new segment
        }
    }

    @segs = grep { $_->[0] eq $contig } @segs  if defined $contig;

    map { [ @$_[0..2], min_max( @$_[1,2] ) ] } @segs;
}


sub min_max { $_[0] <= $_[1] ? ( $_[0], $_[1] ) : ( $_[1], $_[0] ) }



#  These will be part of the fix to genes_in_region.  -- GJO
#  Well, manybe not.

=head3 regions_spanned

    my ( [ $contig, $beg, $end ], ... ) = $fig->regions_spanned( $loc );

or

    my ( [ $contig, $beg, $end ], ... ) = FIG::regions_spanned( $loc );

The location of a feature in a scalar context is

    contig_b1_e1, contig_b2_e2, ...   [one contig_b_e for each segment]

This routine takes as input a scalar location in the above form
and reduces it to one or more regions spanned by the gene. This
involves combining regions in the location string that are on the
same contig and going in the same direction. Unlike L</boundaries_of>,
which returns one region in which the entire gene can be found,
B<regions_spanned> handles wrapping through the orgin, features
split over contigs and exons that are not ordered nicely along
the chromosome (ugly but true).

=over 4

=item loc

The location string for a feature.

=item RETURN

Returns a list of list references. Each inner list contains a contig ID, a starting
position, and an ending position. The starting position may be numerically greater
than the ending position (which indicates a backward-traveling gene). It is
guaranteed that the entire feature is covered by the regions in the list.

=back

=cut

sub regions_spanned {
    shift if UNIVERSAL::isa( $_[0], __PACKAGE__ );
    my( $location ) = ( @_ == 1 ) ? $_[0] : $_[1];
    defined( $location ) || return undef;

    my @regions = ();

    my ( $cur_contig, $cur_beg, $cur_end, $cur_dir );
    my ( $contig, $beg, $end, $dir );
    my @segs = split( /\s*,\s*/, $location );  # should not have space, but ...
    @segs || return undef;

    #  Process the first segment

    my $seg = shift @segs;
    ( ( $cur_contig, $cur_beg, $cur_end ) = ( $seg =~ /^(\S+)_(\d+)_(\d+)$/ ) )
       || return undef;
    $cur_dir = ( $cur_end >= $cur_beg ) ? 1 : -1;

    foreach $seg ( @segs ) {
        ( ( $contig, $beg, $end ) = ( $seg =~ /^(\S+)_(\d+)_(\d+)$/ ) ) || next;
        $dir = ( $end >= $beg ) ? 1 : -1;

        #  Is this a continuation?  Update end

        if ( ( $contig eq $cur_contig )
        && ( $dir == $cur_dir )
        && ( ( ( $dir > 0 ) && ( $end > $cur_end ) )
            || ( ( $dir < 0 ) && ( $end < $cur_end ) ) )
        )
        {
            $cur_end = $end;
        }

        #  Not a continuation.  Report previous and update current.

        else
        {
            push @regions, [ $cur_contig, $cur_beg, $cur_end ];
            ( $cur_contig, $cur_beg, $cur_end, $cur_dir )
            = ( $contig, $beg, $end, $dir );
        }
    }

    #  There should alwasy be a valid, unreported region.

    push @regions, [ $cur_contig, $cur_beg, $cur_end ];

    return wantarray ? @regions : \@regions;
}

=head3 filter_regions

    my  @regions = FIG::filter_regions( $contig, $min, $max,  @regions );

or

    my \@regions = FIG::filter_regions( $contig, $min, $max,  @regions );

or

    my @regions = FIG::filter_regions( $contig, $min, $max, \@regions );

or

    my \@regions = FIG::filter_regions( $contig, $min, $max, \@regions );

Filter a list of regions to those that overlap a specified section of a
particular contig. Region definitions correspond to those produced
by L</regions_spanned>.  That is, C<[>I<contig>C<,>I<beg>C<,>I<end>C<]>.
In the function call, either I<$contig> or I<$min> and I<$max> can be
undefined (permitting anything). So, for example,

    my @regions = FIG::filter_regions(undef, 1, 5000, $regionList);

would return all regions in C<$regionList> that overlap the first
5000 base pairs in any contig. Conversely,

    my @regions = FIG::filter_regions('NC_003904', undef, undef, $regionList);

would return all regions on the contig C<NC_003904>.

=over 4

=item contig

ID of the contig whose regions are to be passed by the filter, or C<undef>
if the contig doesn't matter.

=item min

Leftmost position of the region used for filtering. Only regions which contain
at least one base pair at or beyond this position will be passed. A value
of C<undef> is equivalent to zero.

=item max

Rightmost position of the region used for filtering. Only regions which contain
at least one base pair on or before this position will be passed. A value
of C<undef> is equivalent to the length of the contig.

=item regionList

A list of regions, or a reference to a list of regions. Each region is a
reference to a three-element list, the first element of which is a contig
ID, the second element of which is the start position, and the third
element of which is the ending position. (The ending position can be
before the starting position if the region is backward-traveling.)

=item RETURN

In a scalar context, returns a reference to a list of the filtered regions.
In a list context, returns the list itself.

=back

=cut

sub filter_regions {
    my ( $contig, $min, $max, @regions ) = @_;

    @regions || return ();
    ( ref( $regions[0] ) eq "ARRAY" ) || return undef;

    #  Is it a region list, or a reference to a region list?

    if ( ref( $regions[0]->[0] ) eq "ARRAY" ) { @regions = @{ $regions[0] } }

    if ( ! defined( $contig ) )
    {
        ( defined( $min ) && defined( $max ) ) || return undef;
    }
    else       # with a defined contig name, allow undefined range
    {
        defined( $min ) || ( $min =          1 );
        defined( $max ) || ( $max = 1000000000 );
    }
    ( $min <= $max ) || return ();

    my ( $c, $b, $e );
    my @filtered = grep { ( @$_ >= 3 )               #  Allow extra fields?
                       && ( ( $c, $b, $e ) = @$_ )
                       && ( ( ! defined( $contig ) ) || ( $c eq $contig ) )
                       && ( ( $e >= $b ) || ( ( $b, $e ) = ( $e, $b ) ) )
                       && ( ( $b <= $max ) && ( $e >= $min ) )
                        } @regions;

    return wantarray ? @filtered : \@filtered;
}

=head3 close_genes

    my @features = $fig->close_genes($fid, $dist);

Return all features within a certain distance of a specified other feature.

This method is a quick way to get genes that are near another gene. It calls
L</boundaries_of> to get the boundaries of the incoming gene, then passes
the region computed to L</genes_in_region>.

So, for example, if the specified I<$dist> is 500, the method would select
a region that extends 500 base pairs to either side of the boundaries for
the gene I<$fid>, and pass it to C<genes_in_region> for analysis. The
features returned would be those that overlap the selected region. Note
that the flaws inherent in C<genes_in_region> are also inherent in this
method: if a feature is more than 10000 base pairs long, it may not
be caught even though it has an overlap in the specified region.

=over 4

=item fid

ID of the relevant feature.

=item dist

Desired maximum distance.

=item RETURN

Returns a list of feature IDs for genes that overlap or are close to the boundaries
for the specified incoming feature.

=back

=cut

sub close_genes {
    my($self,$fid,$dist) = @_;

#   warn "In close_genes, self=$self, fid=$fid";

    my $loc = $self->feature_location($fid);
    if ($loc)
    {
        my($contig,$beg,$end) = $self->boundaries_of($loc);
        if ($contig && $beg && $end)
        {
            my $min = &min($beg,$end) - $dist;
            my $max = &max($beg,$end) + $dist;
            my $feat;
            ($feat,undef,undef) = $self->genes_in_region(&FIG::genome_of($fid),$contig,$min,$max);
            return @$feat;
        }
    }
    return ();
}

## returns (Subsys,[pegs]) for the Subsystem containing the most pegs

sub close_in_subsystem {
    my($self,$peg) = @_;

    my %close;
    my %sub = map { $_ => 1 } grep { $self->usable_subsystem($_) } $self->peg_to_subsystems($peg);
    foreach my $peg1 ($self->close_genes($peg,5000))
    {
	next if ($peg eq $peg1);
	foreach my $sub1 ($self->peg_to_subsystems($peg1))
	{
	    if ($sub{$sub1})
	    {
		$close{$sub1}->{$peg1} = 1;
	    }
	}
    }

    my $sofar;
    my $best;
    foreach my $sub1 (keys(%close))
    {
	my $x = $close{$sub1};
	my @others = keys(%$x);
	if (@others > $best)
	{
	    $sofar = [@others];
	    $best = $sub1;
	}
    }
    return $best ? ($best,$sofar) : undef;
}

=head3 adjacent_genes

    my ($left_fid, $right_fid) = $fig->adjacent_genes($fid, $dist);

Return the IDs of the genes immediately to the left and right of a specified
feature.

This method gets a list of the features within the specified distance of
the incoming feature (using L</close_genes>), and then chooses the two
closest to the feature found. If the incoming feature is on the + strand,
these are features to the left and the right. If the incoming feature is
on the - strand, the features will be returned in reverse order.

=over 4

=item fid

ID of the feature whose neighbors are desired.

=item dist

Maximum permissible distance to the neighbors.

=item RETURN

Returns a two-element list containing the IDs of the features on either side
of the incoming feature.

=back

=cut

sub adjacent_genes
{
    my ($self, $fid, $dist) = @_;
    my (@close, $strand, $i);

#   warn "In adjacent_genes, self=$self, fid=$fid";


    $strand = $self->strand_of($fid);

    $dist   = $dist || 2000;
    @close  = $self->close_genes($fid, $dist);
    for ($i=0; $i < @close; ++$i) { last if ($close[$i] eq $fid); }

    # RAE note that if $self->strand_of($close[$i-1]) ne $strand then left/right neighbors
    # were never set! oops!

    # I think the concept of Left and right is confused here. In my mind, left and right
    # are independent of strand ?? E.g. take a look at PEG fig|196600.1.peg.1806
    # this is something like
    #
    #     ---> <--1805--- --1806-->   <--1807--  <----
    #
    # 1805 is always the left neighbor, no?

    my ($left_neighbor, $right_neighbor) = ($close[$i-1], $close[$i+1]);

#    if (0) # this was if ($i > 0) I just skip this whole section!
#    {
#       if ($self->strand_of($close[$i-1]) eq $strand) { $left_neighbor  = $close[$i-1]; }
#    }

    if ($i < $#close)
    {
        if ($self->strand_of($close[$i+1]) eq $strand) { $right_neighbor = $close[$i+1]; }
    }

    # ...return genes in transcription order...
    if ($strand eq '-')
    {
        ($left_neighbor, $right_neighbor) = ($right_neighbor, $left_neighbor);
    }

    return ($left_neighbor, $right_neighbor) ;
}

=head3 compute_genome_similarity

Compute a rough estimate of "similarity" between genomes using the following algorithm:

	1.  You need at least five "genes" from each genome (let's work with incomplete as well as complete).  You get these by

		a. Taking up to 5 of the "universal genes"
		b. supplemented by genes (starting from 1) that are greater than 300 aa

	2.  For each gene from the set consider the set of similarities for it.
		For each match that covers over 200 aa of the gene,

			if the % identify > 70, count a "too-similar{Genome2}"
			else count a "not-too-similar{Genome2}"

	     For each Genome2, if the "too-similar{Genome2}" count > "not-too-similar{Genome2}" count,
				the Genome1-Genome2 matches are too similar.
	  	   else, they are not

Used for filtering candidate PCHs in remove_clustered_pchs2.pl.

=over 4

=item univ_hash

Hash where the keys are the annotations for the universal proteins to be used
in the similarity computation.

=item match_len

Minimum length of similarity match required to be considered for genome similarity.

=item num_genes

Number of genes to consider for the com.putation.

=item RETURN

List of lists of the form [genome2, is-similar, count of too-similar hits, count of not-too-similar hist]

=back

=cut

sub compute_genome_similarity
{
    my($fig, $genome, $univ_hash, $match_len, $num_genes) = @_;

    #
    # set defaults
    #

    $match_len = 200 unless defined($match_len);
    $num_genes = 5 unless defined($num_genes);

    if (!defined($univ_hash))
    {
	my @univ = (
		    "Phenylalanyl-tRNA synthetase beta chain (EC 6.1.1.20)",
		    "Prolyl-tRNA synthetase (EC 6.1.1.15)",
		    "Phenylalanyl-tRNA synthetase alpha chain (EC 6.1.1.20)",
		    "Histidyl-tRNA synthetase (EC 6.1.1.21)",
		    "Arginyl-tRNA synthetase (EC 6.1.1.19)",
		    "Tryptophanyl-tRNA synthetase (EC 6.1.1.2)",
		    "Preprotein translocase secY subunit (TC 3.A.5.1.1)",
		    "Tyrosyl-tRNA synthetase (EC 6.1.1.1)",
		    "Methionyl-tRNA synthetase (EC 6.1.1.10)",
		    "Threonyl-tRNA synthetase (EC 6.1.1.3)",
		    "Valyl-tRNA synthetase (EC 6.1.1.9)"
		   );
	$univ_hash = {};
	map { $univ_hash->{$_} = 1 } @univ;
    }

    my $dbh = $fig->db_handle();

    my @genes;

    #
    # First try to find genes in the universal list.
    #

    my $univ_genes = $fig->find_features_by_annotation($genome, $univ_hash);

    @genes = values %$univ_genes;
#    print "found " . int(@genes) . " universal genes\n";
    if (@genes < $num_genes)
    {
	#
	# Need to pull some genes from the beginning of the genome.
	#
	push(@genes, $fig->find_features_from_start_of_genome('peg', $genome, $num_genes - @genes, 300));
    }
    $#genes = $num_genes - 1 if @genes > $num_genes;

    my @sims = $fig->sims(\@genes, undef, undef, 'fig');
    my(%counts);
    for my $sim (@sims)
    {
	next unless $sim->ln1 > $match_len;
	my $g2 = $fig->genome_of($sim->id2);
	if ($sim->iden > 70)
	{
	    $counts{$g2}->{1}++;
	}
	else
	{
	    $counts{$g2}->{0}++;
	}
    }
    my @out;
    for my $g2 (sort keys %counts)
    {
	my $too_count = $counts{$g2}->{1};
	my $not_too_count = $counts{$g2}->{0};
	my $too_similar = ($too_count > $not_too_count) ? 1 : 0;
	push(@out, [$g2, $too_similar, $too_count, $not_too_count]);
#	print "$g2: $too_similar $too_count $not_too_count\n";
    }
    return @out;
}

sub find_features_from_start_of_genome
{
    my($fig, $ftype, $genome, $num, $min_len) = @_;

    my @genes;

    my @all_pegs = $fig->all_features($genome, $ftype);
    while (@genes < $num and @all_pegs)
    {
	my $peg = shift @all_pegs;
	my $loc = $fig->feature_location($peg);
	my $len = abs($fig->beg_of($loc) - $fig->end_of($loc));
	if ($len > $min_len)
	{
	    push(@genes, $peg);
	}
    }

    return @genes;
}


sub find_features_by_annotation
{
    my($fig, $genome, $anno_hash) = @_;

    my $af = "$FIG_Config::organisms/$genome/assigned_functions";

    my $res = {};

    if (!open(F, "<$af"))
    {
	warn "cannot open $af: $!\n";
	return $res;
    }

    while (<F>)
    {
	chomp;
	my($id, $func) = split(/\t/);

	if ($anno_hash->{$func})
	{
	    $res->{$func} = $id;
	}
    }
    return $res;
}

=head3 feature_location

    my $loc = $fig->feature_location($fid);

or

    my @loc = $fig->feature_location($fid);;

Return the location of a feature. The location consists
of a list of (contigID, begin, end) triples encoded
as strings with an underscore delimiter. So, for example,
C<NC_002755_100_199> indicates a location starting at position
100 and extending through 199 on the contig C<NC_002755>. If
the location goes backward, the start location will be higher
than the end location (e.g. C<NC_002755_199_100>).

In a scalar context, this method returns the locations as a
comma-delimited string

    NC_002755_100_199,NC_002755_210_498

In a list context, the locations are returned as a list

    (NC_002755_100_199, NC_002755_210_498)

=over 4

=item fid

ID of the feature whose location is desired.

=item RETURN

Returns the locations of a feature, either as a comma-delimited
string or a list.

=back

=cut

sub feature_location :Scalar :List {
    my($self,$feature_id) = @_;
    my($relational_db_response,$locations,$location);

#   warn "In feature_location, self=$self, fid=$feature_id";

    if ($self->is_deleted_fid($feature_id)) { return undef }

    $locations = $self->cached('_location');
    if (! ($location = $locations->{$feature_id})) {

	if ($self->{memcache})
	{
	    $locations->{$feature_id} = $location = $self->{memcache}->get("l:$feature_id");
	    #	    print STDERR "got from memcache $location\n";
	}

	if (!$location)
	{
	    my $rdbH = $self->db_handle;
	    if (($relational_db_response = $rdbH->SQL("SELECT location FROM features WHERE ( id = \'$feature_id\' )")) &&
		(@$relational_db_response == 1)) {
                $locations->{$feature_id} = $location = $relational_db_response->[0]->[0];
		# print STDERR "got from db $location\n";
		$self->{memcache}->set("l:$feature_id", $location) if $self->{memcache};
	    }
	}
    }

    if ($location) {
        return wantarray() ? split(/,/,$location) : $location;
    }
    return undef;
}

sub feature_location_bulk {
    my($self,$feature_ids) = @_;

    my @out;
    my @ids;
    my $locations = $self->cached('_location');
    for my $id (@$feature_ids)
    {
	if ($self->is_deleted_fid($id))
	{
	    push(@out, [$id, undef]);
	}
	elsif (my $location = $locations->{$id})
	{
	    push(@out, [$id, $location]);
	}
	else
	{
	    push(@ids, $id);
	}
    }
    #
    # Search in memcache first.
    #
    if ($self->{memcache})
    {
	my $mcout  = $self->{memcache}->get_multi(map { "l:$_" } @ids);
#	print "memcache get " . Dumper($mcout);
	push @out, map { my $k = $_; s/^l://;  [$_, $mcout->{$k}] } keys %$mcout;

	@ids = grep { !$mcout->{"l:$_"} } @ids;
    }

    #
    # Query for any remaining.
    #

    my $dbh = $self->db_handle->{_dbh};

    my @update;
    my %need = map { $_ => 1 } @ids;
    while (@ids)
    {
	my @batch = splice(@ids, 0, 1000);

	my $prots = join(", ", map { $dbh->quote($_) } @batch);

	my $r = $dbh->selectall_arrayref(qq(SELECT id, location
					    FROM features
					    WHERE id in ($prots)));
	for my $row (@$r)
	{
	    my($id, $loc) = @$row;
	    $locations->{$id} = $loc;
	    push(@out, [$id, $loc]);
	    push(@update, ["l:$id", $loc]);
	    delete $need{$id};
        }
    }
    push @out, map { [$_ => undef] } keys %need;

    if ($self->{memcache} && @update)
    {
#	print "memcache put " . Dumper(\@update);
	$self->{memcache}->set_multi(@update);
    }

    return @out;
}

=head3 contig_of

    my $contigID = $fig->contig_of($location);

Return the ID of the contig containing a location.

This method only works with SEED-style locations (I<contigID>C<_>I<beg>C<_>I<end>).
For more comprehensive location parsing, use the B<Location> object.

=over 4

=item location

A SEED-style location (I<contigID>C<_>I<beg>C<_>I<end>), or a comma-delimited list
of SEED-style locations. In the latter case, only the first location in the list will
be processed.

=item RETURN

Returns the contig ID from the first location in the incoming string.

=back

=cut

sub contig_of
{
    my ($self, $locus) = @_;

    $locus =~ m/^([^,]+)_\d+_\d+/;

    return $1;
}

=head3 beg_of

    my $beg = $fig->beg_of($location);

Return the beginning point of a location.

This method only works with SEED-style locations (I<contigID>C<_>I<beg>C<_>I<end>).
For more comprehensive location parsing, use the B<Location> object.

=over 4

=item location

A SEED-style location (I<contigID>C<_>I<beg>C<_>I<end>), or a comma-delimited list
of SEED-style locations. In the latter case, only the first location in the list will
be processed.

=item RETURN

Returns the beginning point from the first location in the incoming string.

=back

=cut

sub beg_of
{
    my ($self, $locus) = @_;

    $locus =~ m/^[^,]+_(\d+)_\d+/;

    return $1;
}

=head3 end_of

    my $end = $fig->end_of($location);

Return the ending point of a location.

This method only works with SEED-style locations (I<contigID>C<_>I<beg>C<_>I<end>).
For more comprehensive location parsing, use the B<Location> object.

=over 4

=item location

A SEED-style location (I<contigID>C<_>I<beg>C<_>I<end>), or a comma-delimited list
of SEED-style locations. In the latter case, only the first location in the list will
be processed.

=item RETURN

Returns the contig ID from the first location in the incoming string.

=back

=cut

sub end_of
{
    my ($self, $locus) = @_;

    $locus =~ m/\S+_\d+_(\d+)$/;

    return $1;
}

=head3 upstream_of

    my $dna = $fig->upstream_of($peg, $upstream, $coding);

Return the DNA immediately upstream of a feature. This method contains code lifted from
the C<upstream.pl> script.

=over 4

=item peg

ID of the feature whose upstream DNA is desired.

=item upstream

Number of base pairs considered upstream.

=item coding

Number of base pairs inside the feature to be included in the upstream region.

=item RETURN

Returns the DNA sequence upstream of the feature's begin point and extending into the coding
region. Letters inside a feature are in upper case and inter-genic letters are in lower case.
A hyphen separates the true upstream letters from the coding region.

=back

=cut
#: Return Type $;
sub upstream_of {
    # Get the parameters.
    my ($self, $peg, $upstream, $coding) = @_;
    # Declare the work variables.
    my ($gene_before, $inter_genic, $c_seq) = ("", "", "");
    # Compute the upstream.
    my ($contig,$beg,$end) = $self->boundaries_of(scalar $self->feature_location($peg));
    my $genome = $self->genome_of($peg);
    my $retVal = FIGRules::Upstream($self, $genome, "${contig}_${beg}_${end}", $upstream, $coding);
    # Return the result.
    return $retVal;
}

=head3 strand_of

    my $strand = $fig->contig_of($location);

Return the strand (C<+> or C<->) of a location.

This method only works with SEED-style locations (I<contigID>C<_>I<beg>C<_>I<end>).
For more comprehensive location parsing, use the B<Location> object.

=over 4

=item location

A comma-delimited list of SEED-style location (I<contigID>C<_>I<beg>C<_>I<end>).

=item RETURN

Returns C<+> if the list describes a forward-oriented location, and C<-> if the list
described a backward-oriented location.

=back

=cut

sub strand_of
{
    my ($self, $fid) = @_;
    my ($beg, $end);

#   warn "In strand_of, self=$self, fid=$fid";

    (undef, $beg, $end) = $self->boundaries_of($self->feature_location($fid));

    if ($beg < $end) { return '+'; } else { return '-'; }
}

=head3 find_contig_with_checksum

    my $contigID = $fig->find_contig_with_checksum($genome, $checksum);

Find a contig in the given genome with the given checksum.

This method is useful for determining if a particular contig has already been
recorded for the given genome. The checksum is computed from the contig contents,
so a matching checksum indicates that the contigs may have the same content.

=over 4

=item genome

ID of the genome whose contigs are to be examined.

=item checksum

Checksum value for the desired contig.

=item RETURN

Returns the ID of a contig in the given genome that has the caller-specified checksum,
or C<undef> if no such contig exists.

=back

=cut

sub find_contig_with_checksum
{
    my($self, $genome, $checksum) = @_;

    #
    # This implementation scans all the contig files for the organism; when
    # we have contig checksums in the database this will simplify
    # significantly.
    #
    # For some efficiency, we cache the checksums we compute along the way since
    # it's probably likely we'll poke at other contigs for this organism.
    #

    my $gdir = "$FIG_Config::organisms/$genome";

    my $cached_cksum = $self->cached('_contig_checksum');

    if (opendir(my $dh, $gdir))
    {
        for my $file (map { "$gdir/$_" } grep { $_ =~ /^contigs\d*$/ } readdir($dh))
        {
            local $/ = "\n>";

            if (open(my $fh, "<$file"))
            {
                while (<$fh>)
                {
                    chomp;

                    #
                    # Pull the contig identifier from the first line.
                    # We need the >? to handle the first line in the file;
                    # the others are removed by the chomp above because
                    # $/ is set to "\n>";
                    #

                    if (s/^>?\s*(\S+)([^\n]*)\n//)
                    {
                        my $ident = $1;
                        my $contig_txt = $_;

                        $contig_txt =~ s/\s//sg;
                        $contig_txt = uc($contig_txt);

                        #
                        # See if we have a cached value.
                        #

                        my $this_checksum;
                        $this_checksum = $cached_cksum->{$genome, $ident};
                        if (!$this_checksum)
                        {

                            my($rd, $wr, $pid);

                            if (!($pid = open2($rd, $wr, "cksum")))
                            {
                                Confess("Cannot run open2 cksum: $!");
                            }

                            $wr->write($contig_txt, length($contig_txt));

                            close($wr);

                            $_ = <$rd>;
                            close($rd);
                            waitpid $pid, 0;
                            chomp;

                            my @vals = split(/\s+/, $_);
                            $this_checksum = $vals[0];
                            $cached_cksum->{$genome, $ident} = $this_checksum;
                        }
                        if ($this_checksum == $checksum)
                        {
                            return $ident;
                        }
                    }
                }
            }
        }
    }
}

=head3 contig_checksum

    my $checksum = $fig->contig_checksum($genome, $contig);

or

    my @checksum = $fig->contig_checksum($genome, $contig);

Return the checksum of the specified contig. The checksum is computed from the
contig's content in a parallel process. The process returns a space-delimited list
of numbers. The numbers can be split into a real list if the method is invoked in
a list context. For b

=cut

sub contig_checksum
{
    my($self, $genome, $contig) = @_;

    my $contig_txt = $self->read_contig($genome, $contig);

    $contig_txt =~ s/\s//sg;
    $contig_txt = uc($contig_txt);

    my($rd, $wr, $pid);

    if (!($pid = open2($rd, $wr, "cksum")))
    {
        Confess("Cannot run open2 cksum: $!");
    }

    $wr->write($contig_txt, length($contig_txt));

    close($wr);

    $_ = <$rd>;
    close($rd);
    waitpid $pid, 0;

    chomp;
    my @vals = split(/\s+/, $_);
    if (wantarray)
    {
        return @vals;
    }
    else
    {
        return $vals[0];
    }
}

=head3 read_contig

Read a single contig from the contigs file.

=cut
sub read_contig
{
    my($self, $genome, $contig) = @_;

    #
    # Read the contig. The database has it in a set of chunks, but we
    # just use the seek to the starting point, and read up to the next "\n>".
    #

    my $ret = $self->db_handle->SQL(qq!SELECT fileno, seek FROM contig_seeks
                                       WHERE genome = '$genome' and contig = '$contig' and
                                       startn = 0!);
    if (!$ret or @$ret != 1)
    {
        return undef;
    }

    my($fileno, $seek) = @{$ret->[0]};
    my $file = $self->N2file($fileno);

    my($fh, $contig_txt);

    if (!open($fh, "<$file"))
    {
        warn "contig_checksum: could not open $file: $!\n";
        return undef;
    }

    seek($fh, $seek, 0);

    {
        local $/ = "\n>";

        $contig_txt = <$fh>;
        chomp($contig_txt);
    }

    return $contig_txt;
}

=head3 boundaries_of

usage: ($contig,$beg,$end,$strand) = $fig->boundaries_of($loc)

The location of a feature in a scalar context is

    contig_b1_e1,contig_b2_e2,...   [one contig_b_e for each exon]

This routine takes as input such a location and reduces it to a single
description of the entire region containing the gene.

=cut

sub boundaries_of {
    if ( UNIVERSAL::isa($_[0],__PACKAGE__) ) {
        shift;
    } else {
        # my ($package, $filename, $line) = caller;
        # warn "Deprecated boundaries_of called from $filename line $line.";
    }

    #  This test is almost certainly unnecessary, given the shift above
    #  If it is not so, then there is a bad violation of the call format
    # my $location = ( (@_ == 1) ? $_[0] : $_[1] );
    my $location = $_[0] || '';  # this makes location defined, no matter what.
    $location =~ s/\s+//g;       # a more general solution to whitespace
    my @exons = split( /,/, $location );
    my( $contig1, $beg ) = $exons[ 0] =~ /^(\S+)_(\d+)_\d+$/;
    my( $contig2, $end ) = $exons[-1] =~ /^(\S+)_\d+_(\d+)$/;

    #  If $beg and $end are true, then $contig1 and $contig2 are defined
    if ( ! ( $beg && $end && ($contig1 eq $contig2) ) )
    {
        Cluck("Could not parse loc=$location.") if T(0);
        return ();
    }

    my $strand = ($beg <= $end) ? qq(+) : qq(-);
    return ($contig1, $beg, $end, $strand);
}


#  Return the number of overlapping nucleotides between 2 location values

sub location_overlap
{
    if ( UNIVERSAL::isa( $_[0], __PACKAGE__ ) ) { shift }

    my( $loc1, $loc2 ) = @_;
    ( defined( $loc1 ) && defined( $loc2 ) ) or return undef;

    my ( $c1, $b1, $e1 ) = &boundaries_of( $loc1 );
    my ( $c2, $b2, $e2 ) = &boundaries_of( $loc2 );

    my $overlap = 0;
    if ( $c1 eq $c2 )
    {
        my ( $min1, $max1 ) = &minmax( $b1, $e1 );
        my ( $min2, $max2 ) = &minmax( $b2, $e2 );
        if ( ( $min1 <= $max2 ) && ( $max1 >= $min2 ) )
        {
            $overlap = &min( $max1, $max2 ) - &max( $min1, $min2 ) + 1;
        }
    }

    return $overlap;
}

#  Return the minimum followed by the maximum of a pair of numbers

sub minmax { $_[0] < $_[1] ? @_[0,1] : @_[1,0] }


=head3 boundaries_of_2

    \@regions = $fig->boundaries_of_2( $location );

     @regions = $fig->boundaries_of_2( $location );

Locations can be a list of intervals (contig_beg_end), but the intervals
need not be on a single contig, contiguous, or in a consistent orientation
(e.g., a feature that wraps from the end to the beginning of a genome, or a
trans-spliced protein).  This function defines a region of a gene a sequence
parts of the location that are on the same contig, in the same orientation,
and with end points that progress along the contig in the same direction as
the individual parts.  This function is a generalization of boundaries_of().
The latter function returns undef if the first and last contigs are not
the same, and returns a location spanning nearly the entire contig it the
location spans the origin.

=over 4

=item location

contig1_beg1_end1,contig2_beg2_end2,...

=item @regions

   ( [contig, beg, end], [contig, beg, end], ... )

where consecutive location intervals with a common contig, direction, and
consistent direction of progression along the contig are merged.  The vast
majority of genes will be reduced to the a single region, which is that
returned by boundaries_of().  That is, most of the time:

    boundaries_of( $loc )

is the same as

    @{ boundaries_of_2( $loc )->[0] || [] }

=back

=cut

sub boundaries_of_2
{
    if ( UNIVERSAL::isa( $_[0], __PACKAGE__ ) ) { shift }

    my ( $loc ) = @_;
    my @parts = grep { $_->[0] && $_->[1] && $_->[2] }
                map  { [ /^(.+)_(\d+)_(\d+)$/ ] }
                split /,/, $loc;
    @parts or return wantarray ? () : [];

    my ( $contig, $beg, $end ) = @{ shift @parts };
    my $dir = ( $end <=> $beg );
    my @regions;
    foreach ( @parts )
    {
        my ( $c1, $b1, $e1 ) = @$_;
        my $d1 = ( $e1 <=> $b1 );
        if ( ( $c1 eq $contig ) && ( $d1 == $dir ) && ( ( $e1 <=> $end ) == $dir ) )
        {
            $end = $e1;
        }
        else
        {
            push @regions, [ $contig, $beg, $end ];
            ( $contig, $beg, $end, $dir ) = ( $c1, $b1, $e1, $d1 );
        }
    }
    push @regions, [ $contig, $beg, $end ];
    wantarray() ? @regions : \@regions;
}


=head3 all_features_detailed

    my $featureList = $fig->all_features_detailed($genome);

Returns a list of all features in the designated genome, with their location, alias,
and type information included. This is used in the GenDB import and Sprout load to
speed up the process.

Deleted features are not returned!

=over 4

=item genome

ID of the genome whose features are desired.

=item RETURN

Returns a reference to a list of tuples. Each tuple consists of four elements: (1) the feature
ID, (2) the feature location (as a comma-delimited list of location specifiers), (3) the feature
aliases (as a comma-delimited list of named aliases), and (4) the feature type.

=back

=cut
#: Return Type $@@;
sub all_features_detailed {
    my($self,$genome) = @_;

    my $rdbH = $self->db_handle;
    my $relational_db_response = $rdbH->SQL("SELECT id, location, aliases, type FROM features WHERE  (genome = \'$genome\')");
        my @features;
    foreach my $tuple (@$relational_db_response)
    {
        push @features, $tuple unless ($self->is_deleted_fid($tuple->[0]));
    }
    @features = sort { &FIG::by_fig_id($a->[0],$b->[0]) } @features;
    return \@features;
}

=head3 all_features_detailed_fast

    my $featureList = $fig->all_features_detailed($genome, $min, $max, $contig);

Returns a list of all features in the designated genome, with various useful information
included.

Deleted features are not returned!

=over 4

=item genome

ID of the genome whose features are desired.

=item min (optional)

If specified, the minimum contig location of interest. Features not entirely to the right
of this location are ignored.

=item max (optional)

If specified, the maximum contig location of interest. Features not entirely to the left
of this location are ignore.

=item contig (optional)

If specified, the contig of interest. Features not on this contig are ignored.

=item RETURN

Returns a reference to a list of tuples. Each tuple consists of four elements: (1) the feature
ID, (2) the feature location (as a comma-delimited list of location specifiers), (3) the feature
aliases (as a comma-delimited list of named aliases), (4) the feature type, (5) the leftmost
index of the feature's first location, (6) the rightmost index of the feature's last location,
(7) the current functional assignment, (8) the user who made the assignment, and (9) the
quality of the assignment (which is usually an empty string).

=back

=cut

# does the same as the above, except using the advantage of a join statement
# and including minloc and maxloc as well as the function, annotator and quality
sub all_features_detailed_fast {
  my($self,$genome, $min, $max, $contig) = @_;

  my $minmax = "";
  if (defined($min) && defined($max)) {
    $minmax = "AND ((minloc <= $min AND maxloc >= $min) OR (minloc <= $max AND maxloc >= $max) OR (minloc >= $min AND maxloc <= $max)) ";
  }

  my $contig_line = "";
  if (defined($contig)) {
    $contig_line = "AND features.contig = '" . $contig . "' ";
  }

  my $rdbH = $self->db_handle;

  my $relational_db_response = $rdbH->SQL("SELECT DISTINCT id, location, aliases, type, minloc, maxloc, assigned_function, made_by, quality FROM (SELECT id, location, aliases, type, minloc, maxloc FROM features LEFT OUTER JOIN deleted_fids ON features.id = deleted_fids.fid WHERE features.genome = '" . $genome . "' " . $contig_line . $minmax . "AND fid IS NULL) AS t1 LEFT OUTER JOIN assigned_functions on t1.id = assigned_functions.prot");

  return $relational_db_response || ();

# SELECT id, location, aliases, type, minloc, maxloc, assigned_function, made_by, quality FROM (SELECT id, location, aliases, type, minloc, maxloc FROM features LEFT OUTER JOIN deleted_fids ON features.id = deleted_fids.fid WHERE features.genome = '83333.1' AND ((minloc < 1 AND maxloc > 1) OR (minloc < 4639221 AND maxloc > 4639221) OR (minloc > 1 AND maxloc < 4639221)) AND fid IS NULL) AS t1 LEFT OUTER JOIN assigned_functions on t1.id = assigned_functions.prot;
}


sub contig_lengths {
  my ($self, $genome) = @_;

  my $contig_lengths;

  my $rdbH = $self->db_handle;
  my $relational_db_response = $rdbH->SQL("SELECT contig, len FROM contig_lengths WHERE genome='$genome'");

  foreach my $contig (@$relational_db_response) {
    $contig_lengths->{$contig->[0]} = $contig->[1];
  }

  return $contig_lengths;
}

sub contig_lengths_bulk {
  my ($self, $genome_ids) = @_;

  my $contig_lengths = {};

  if (ref($genome_ids) eq 'ARRAY' && @$genome_ids > 0)
  {
      my $in = join(", ", map { "?" } @$genome_ids);

      my $rdbH = $self->db_handle;
      my $relational_db_response = $rdbH->SQL(qq(SELECT genome, contig, len
						 FROM contig_lengths
						 WHERE genome IN ($in)), undef, @$genome_ids);
      for my $ent (@$relational_db_response)
      {
	  my($g, $c, $l) = @$ent;
	  $contig_lengths->{$g}->{$c} = $l;
      }
  }

  return $contig_lengths;
}

=head3 all_features

    my @fidList = $fig->all_features($genome,$type);

Returns a list of all feature IDs of a specified type in the designated genome.  You would
usually use just

    $fig->pegs_of($genome) or
    $fig->rnas_of($genome)

which simply invoke this routine.

=over 4

=item genome

ID of the genome whose features are desired.

=item type (optional)

Type of feature desired (peg, rna, etc.). If omitted, all features are returned.

=item RETURN

Returns a list of the IDs for the desired features.

=back

=cut

sub all_features {
    my($self,$genome,$type) = @_;

    my $rdbH = $self->db_handle;
    my $where = "(genome = \'$genome\'";
    if ($type) {
        $where .= " AND (type = \'$type\')";
    }
    $where .= ")";
    my $relational_db_response = $rdbH->SQL("SELECT id FROM features WHERE $where");

    if (@$relational_db_response > 0)
    {
        return sort { &FIG::by_fig_id($a,$b) }
	       grep { ! $self->is_deleted_fid($_) } map { $_->[0] } @$relational_db_response;
    }
    return ();
}

=head3 all_valid_protein_features_in_seed

    my $feature_list = $fig->all_valid_protein_features_in_SEED($complete)

Returns a list of all protein features in the current SEED that are not deleted
and do not belong to a deleted genome.

If $complete is true, limit to complete genomes.

Records returned contain

    domain
    complete
    feature_id
    id index in genome
    feature type
    genome id
    location
    contig
    min location on contig
    max location on contig
    aliases
    md5 checksum
    assigned funtion
=cut

sub all_valid_protein_features_in_SEED
{
    my($self, $complete) = @_;

    my @where;

    push(@where, "type = 'peg'");

    if ($complete)
    {
	push(@where, "genome.complete = 1");
    }

    my $where = join(" AND ", map { "($_)" } @where);

    my $rdbH = $self->db_handle();

    my $res = $rdbH->SQL(qq(SELECT genome.maindomain, genome.complete, features.*, md5, assigned_function
			    FROM features JOIN genome ON features.genome = genome.genome
			         LEFT OUTER JOIN deleted_fids ON features.id = deleted_fids.fid
			         JOIN protein_sequence_MD5 ON protein_sequence_MD5.id = features.id
			         LEFT OUTER JOIN assigned_functions ON assigned_functions.prot = features.id
			    WHERE deleted_fids.fid IS NULL AND
			    	$where));
    return $res;
}

sub essentiality_data {
  my($self,$genome,$experiment, $value) = @_;

  my $rdbH = $self->db_handle;

  my $defined_val = "";
  if (defined($value)) {
    if (ref($value) eq "ARRAY") {
      my $vals;
      foreach my $val (@$value) {
	push(@$vals, "val='" . $val . "'");
      }
      $defined_val = " AND (" . join(" OR ", @$vals) . ")";
    } else {
      $defined_val = " AND val='" . $value . "'";
    }
  }

  my $statement = "SELECT prot, aliases, assigned_function, val, minloc FROM (SELECT CONCAT('fig|', genome, '.', ftype, '.', id) AS pid, val FROM attribute WHERE genome='" . $genome . "' AND tag='" . $experiment . "'" . $defined_val . ") AS t1 LEFT OUTER JOIN assigned_functions on t1.pid = assigned_functions.prot LEFT OUTER JOIN features ON t1.pid = features.id ORDER BY minloc";

  my $relational_db_response = $rdbH->SQL($statement);

  my $return;
  foreach my $row (@$relational_db_response) {
      my $retval = $rdbH->SQL("SELECT DISTINCT subsystem from subsystem_index WHERE protein='" . $row->[0] . "'");
      my $subsystems;
      foreach my $subsystem (@$retval) {
	push(@$subsystems, $subsystem->[0]);
      }
      push(@$row, $subsystems || []);
      push(@$return, $row);
  }

  return $return || ();
}

=head3 genome_id_to_genome_object

  my $obj = $fig->genome_id_to_genome_object($genome_id)

Return a KBase-style genome object for the given genome ID.
    
=cut

sub genome_id_to_genome_object
{
    my($self, $genome_id, $options) = @_;

    $options = {} unless ref($options) eq 'HASH';

    if (ref($self) ne 'FIGV' && ! $self->is_genome($genome_id))
    {
	return undef;
    }

    my $gobj = GenomeTypeObject->new();
print STDERR "Load from ", $self->organism_directory($genome_id), "\n";
    my @tax_id;
    if (open(my $tfh, "<", $self->organism_directory($genome_id) . "/TAXONOMY_ID"))
    {
	my $tax_id = <$tfh>;
	chomp $tax_id;
	@tax_id = (taxonomy_id => $tax_id,
		ncbi_taxonomy_id => $tax_id);
    }

    my $gc = 11;
    if (open(my $gcfh, "<", $self->organism_directory($genome_id) . "/GENETIC_CODE"))
    {
	$gc = <$gcfh>;
	chomp $gc;
    }

    $gobj->set_metadata({
	id => $genome_id,
	scientific_name => $self->genus_species($genome_id),
	domain => $self->genome_domain($genome_id),
	source => 'SEED',
	source_id => $genome_id,
	genetic_code => $gc,
	taxonomy => $self->taxonomy_of($genome_id),
        @tax_id,
    });

    my $feats = $self->all_features_detailed_fast($genome_id);

    for my $f (@$feats)
    {
	my($fid, $loc_str, $aliases, $type, $min, $max, $func, $who, $qual) = @$f;
	my $trans = $self->get_translation($fid);
	my @loc = map { my $b = BasicLocation->new($_); [$b->Contig, $b->Begin, $b->Dir, $b->Length] } split(/,/, $loc_str);
	$gobj->add_feature({
	    -id => $fid,
	    -type => $type,
	    -location => [@loc],
	    -function => $func,
	    -annotator => $who,
	    ($trans ? (-protein_translation => $trans) : ()),
	    -aliases => [split(/,/, $aliases)],
	   });
    }

    #
    # We will pull the contigs from the contigs file, not the database, since we
    # need them in their entirety.
    #

    if (!$options->{skip_contigs})
    {
	my $fh;
	if (open($fh, "<", $self->organism_directory($genome_id) . "/contigs"))
	{
	    my $ctgs = [];
	    while (my($id, $def, $seq) = read_next_fasta_seq($fh))
	    {
		push(@$ctgs, { id => $id, dna => $seq });
	    }
	    $gobj->add_contigs($ctgs);
	    close($fh);
	}
	else
	{
	    warn "genome_id_to_genome_object: cannot get contigs for $genome_id: $!";
	}
    }

    #
    # Read closest.genomes if it's there.
    #

    if (open(my $cfh, "<", $self->organism_directory($genome_id) . "/closest.genomes"))
    {
	my $c = [];
	while (<$cfh>)
	{
	    chomp;
	    my($id, $score, $name) = split(/\t/);
	    push(@$c, { genome_id => $id, closeness_measure => $score, genome_name => $name, analysis_method => 'SEED_dir' });
	}
	$gobj->{close_genomes} = $c;
	close($cfh);
    }

    return $gobj->prepare_for_return();
}

=head3 pegs_of

usage: $fig->pegs_of($genome)

Returns a list of all PEGs in the specified genome.  Note that order is not
specified.

=cut

sub pegs_of {
    my($self,$genome) = @_;

    return $self->all_features($genome,"peg");
}


=head3 rnas_of

usage: $fig->rnas_of($genome)

Returns a list of all RNAs for the given genome.

=cut

sub rnas_of {
    my($self,$genome) = @_;

    return $self->all_features($genome,"rna");
}

=head3 feature_aliases

usage: @aliases = $fig->feature_aliases($fid)  OR
       $aliases = $fig->feature_aliases($fid)

Returns a list of aliases (gene IDs, arbitrary numbers assigned by authors, etc.) for the feature.
These must come from the tbl files, so add them there if you want to see them here.

In a scalar context, the aliases come back with commas separating them.

=cut

sub feature_aliases {
    my($self,$feature_id) = @_;
    my($rdbH,$relational_db_response,@aliases,$aliases,%aliases,$x);

    if ($self->is_deleted_fid($feature_id)) { return undef }

    $rdbH = $self->db_handle;
    @aliases = ();
    if (($relational_db_response = $rdbH->SQL("SELECT aliases FROM features WHERE ( id = \'$feature_id\' )")) &&
        (@$relational_db_response == 1))
    {
        $aliases = $relational_db_response->[0]->[0];
        %aliases = map { $_ => 1 } split(/,/,$aliases);
    }

    if (($relational_db_response = $rdbH->SQL("SELECT alias FROM ext_alias WHERE ( id = \'$feature_id\' )")) &&
        (@$relational_db_response > 0))
    {
        foreach $x (@$relational_db_response)
        {
            $aliases{$x->[0]} = 1;
        }
    }

    @aliases = sort grep { not /^xxx\d+$/ } keys(%aliases);

    return wantarray() ? @aliases : join(",",@aliases);
}

sub feature_aliases_in_tbl {
    my($self,$feature_id) = @_;
    my($rdbH,$relational_db_response,@aliases,$aliases,%aliases,$x);

    if ($self->is_deleted_fid($feature_id)) { return undef }

    $rdbH = $self->db_handle;
    @aliases = ();
    if (($relational_db_response = $rdbH->SQL("SELECT aliases FROM features WHERE ( id = \'$feature_id\' )")) &&
        (@$relational_db_response == 1))
    {
        $aliases = $relational_db_response->[0]->[0];
        %aliases = map { $_ => 1 } split(/,/,$aliases);
    }

    @aliases = sort grep { not /^xxx\d+$/ } keys(%aliases);

    return wantarray() ? @aliases : join(",",@aliases);
}

sub feature_aliases_bulk {
    my($self,$id_list,$no_del_check) = @_;
    my($rdbH,$relational_db_response,@aliases,$aliases,%aliases,$x,$id);

    my(@ids);

    if ($no_del_check)
    {
        @ids = @$id_list;
    }
    else
    {
        @ids = grep { not $self->is_deleted_fid($_) } @$id_list;
    }

    my $out = {};
    #
    # Search in memcache first.
    #
    if ($self->{memcache})
    {
	my $mcout  = $self->{memcache}->get_multi(map { "a:$_" } @ids);
	map { s/^a://;  $out->{$_} = $mcout->{"a:$_"} } keys %$mcout;

	@ids = grep { !$out->{$_} } @ids;
    }

    my $cond = join(" or ", map { "id = '$_'" } @ids);

    if($cond){
	$rdbH = $self->db_handle;

	my $res = $rdbH->SQL("SELECT id, aliases FROM features WHERE ( $cond )");

	%aliases = ();
	for my $ent (@$res)
	{
	    ($id, $aliases) = @$ent;
	    map { $aliases{$id}->{$_} = 1 }  split(/,/,$aliases);
	}

	$res = $rdbH->SQL(qq(SELECT id, alias
			     FROM ext_alias
			     WHERE ( $cond )));

	for my $ent (@$res)
	{
	    my $alias;
	    ($id, $alias) = @$ent;
	    $aliases{$id}->{$alias} = 1;
	}

	my @update;
	for my $id (keys(%aliases))
	{
	    $out->{$id} = [sort grep { not /^xxx\d+$/ } keys(%{$aliases{$id}})];

	    if ($self->{memcache})
	    {
		push(@update, ["a:$id", $out->{$id}]);
	    }
	}
	if ($self->{memcache} && @update)
	{
	    my $res = $self->{memcache}->set_multi(@update);
	    # print STDERR  "Updated with memcached: " . Dumper(\@update, $res);
	}
    }

    return $out;
}

sub fids_to_patric
{
    my($self,$id_list) = @_;
    my($rdbH,$relational_db_response,@aliases,$aliases,%aliases,$x,$id);

    my @ids = @$id_list;

    my $out = {};
    #
    # Search in memcache first.
    #
    if ($self->{memcache})
    {
	my $mcout  = $self->{memcache}->get_multi(map { "patric:$_" } @ids);
	map { s/^patric://;  $out->{$_} = $mcout->{"patric:$_"} } keys %$mcout;

	@ids = grep { !$out->{$_} } @ids;
    }

    if (@ids && $self->db_handle->table_exists('patric_map'))
    {
	my $qs = join(", ", map { "?" } @ids);
	my $cond = "fid IN ($qs)";

	my @update;
	$rdbH = $self->db_handle;

	my $res = $rdbH->SQL("SELECT fid, patric FROM patric_map WHERE ( $cond )", undef, @ids);

	for my $ent (@$res)
	{
	    my($id, $patric) = @$ent;
	    $out->{$id} = $patric;
	    push(@update, ["patric:$id", $patric]);
	}

	if ($self->{memcache} && @update)
	{
	    my $res = $self->{memcache}->set_multi(@update);
	    # print STDERR  "Updated with memcached: " . Dumper(\@update, $res);
	}
    }

    return $out;
}

=head3 uniprot_aliases

    my @aliases = $fig->uniprot_aliases($fid)
    OR
    my $aliases = $fig->uniprot_aliases($fid)

Return the uniprot aliases (SwissProt, TREMBL and UniProt) for a PEG.

The aliases returned may be from a different organism than the organism of
the input feature $fid.

A call to get_corresponding_ids is done first and will return the same-sequence
same-genome ids. If none are found, mapped_prot_ids is called which will give the
same-sequence ids.

If you need to know which form of alias is being returned, call these methods directly.

Only one id is returned for every accession found.
Example 1: If both uni|Q8FLC2 and uni|Q8FLC2_ECOL6 are found in the database, only uni|Q8FLC2
will be returned.
Example 2: If sp|P75616 and uni|P75616 are found in the database, only sp|P75616 will be
returned. The order of preference here is sp before tr before uni.

=over 4

=item fid

Feature ID of the PEG whose aliases are desired.

=item RETURN

Depending on the context of the call, either a list of aliases (sp, tr and uni) is
returned, or a comma-separated string. If no aliases are found, the empty list or
string will be returned.

=back

=cut

sub uniprot_aliases {
    my($self, $feature_id) = @_;

    my @unis = ();

    # get_corresponding_ids preserves the organism mapping from fid to uniprot id, so use these results when available
    @unis = map {'uni|' . $_->[0]} grep {$_->[1] eq 'UniProt'} $self->get_corresponding_ids($feature_id, 1);

    # if get_corresponding_ids does not return anything, use the synonyms.
    # this is sequence based, so the organism may differ
    if ( ! @unis ) {
	# get list of [id, length], where id is sp, tr, uni and length is the amino acid sequence length
	@unis = map {$_->[0]} grep {$_->[0] =~ /^(sp|tr|uni)/} $self->mapped_prot_ids($feature_id);
    }

    # create a hash with the ids for every accession, e.g. the id sp|P75616 has accession P75616
    my %uni_acc;
    foreach my $uni ( @unis ) {
	if ( $uni =~ /\|([^_]+)/ ) {
	    push @{ $uni_acc{$1} }, $uni;
	}
    }

    my @unique_acc_unis = ();

    # get a single id for each accession
    # the ids may contain 'sp|', 'tr|' or 'uni|', sorting brings sp to the front followed by tr followed by uni.
    # add only the first id from the sort to the list to be returned.
    foreach my $acc ( keys %uni_acc ) {
	push @unique_acc_unis, (sort @{ $uni_acc{$acc} })[0];
    }

    return wantarray() ? @unique_acc_unis : join(",", @unique_acc_unis);
}

=head3 uniprot_aliases_bulk

    my $hash = $fig->uniprot_aliases_bulk(\@fids, $no_del_check);

Return a hash mapping the specified feature IDs to lists of their uniprot
aliases.

=over 4

=item fids

A list of FIG feature IDs.

=item no_del_check

If TRUE, deleted feature IDs B<will not> be removed from the feature ID list
before processing. The default is FALSE, which means deleted feature IDs
B<will> be removed before processing.

=item RETURN

Returns a hash mapping each feature ID to a list of its uniprot aliases.

=back

=cut

sub uniprot_aliases_bulk {
    my($self,$id_list,$no_del_check) = @_;
    my(@ids);

    if ($no_del_check)
    {
        @ids = @$id_list;
    }
    else
    {
        @ids = grep { not $self->is_deleted_fid($_) } @$id_list;
    }

    my $cond = join(" or ", map { "id = '$_'" } @ids);

    my $rdbH = $self->db_handle;

    my $res = $rdbH->SQL(qq(SELECT id, alias
                            FROM ext_alias
                            WHERE ( $cond ) AND alias like 'uni|%'));

    my %aliases;
    for my $ent (@$res)
    {
        my($id, $alias) = @$ent;
        $aliases{$id}->{$alias} = 1;
    }

    my $out = {};

    for my $id (keys(%aliases))
    {
        $out->{$id} = [sort grep { not /^xxx\d+$/ } keys(%{$aliases{$id}})];
    }

    return $out;
}


############################################
#
#  map SEED internal references for external databases into globally "acceptable" format
#  see: http://www.ncbi.nlm.nih.gov/collab/db_xref.html for format definition
#
sub rewrite_db_xrefs {

    my ($self, $alias) = @_;

    if ( $alias =~ /^gi/ ) { 	# /db_xref="GI:1234567890"
	$alias =~ s/^gi\|/GI:/;
	return $alias;
    }
    elsif ($alias =~ /^uni/) { # /db_xref=" UniProtKB/TrEMBL:Q00177"
	$alias =~ s%uni\|%UniProtKB/TrEMBL:%;
	return $alias
	}
    elsif ($alias =~ /^kegg/i){
	$alias =~ s/kegg\|/KEGG:/i;
	$alias =~ s/^(.*):/$1+/;
	return  $alias
	}
    elsif ($alias =~ /^sp\|/) { # /db_xref="UniProtKB/Swiss-Prot:P12345"
	$alias =~ s%sp\|%UniProtKB/Swiss-Prot:%;
	return  $alias
	}
    else { # unsupported external alias, return empty string
	return '';
    }

}

=head3 rewrite_db_xrefs_brc

Convert an alias to a db_xref. This uses the BRC format db_xref, which is a conglomeration of NCBI, GO, and BioMoby.

This method will return a correctly formatted db_ref if the argument is one of our currently recognized formats, otherwise it returns undef.

This example code should provide the functions you want

foreach my $alias ($fig->feature_aliases($peg))
{
	if (my $dbxref=$fig->rewrite_db_xrefs_brc($alias)) {print "The dbxref is $dbxref\n"}
	else {print "The alias is $alias\n"}
}


For a list of approved dbxrefs, see http://www.brc-central.org/cgi-bin/brc-central/dbxref_list.cgi

=cut

sub rewrite_db_xrefs_brc {

	my ($self, $alias) = @_;

	if ($alias =~ /^COG\:/ || $alias =~ /GeneID\:/ ||  $alias =~ /^CDD\:/ || $alias =~ /^Locus_Tag\:/)
	{ # these are valid db_xrefs and don't need changing
		return $alias;
	}
	elsif ($alias =~ /^NP\_/ || $alias =~ /^YP\_/ || $alias =~ /^ZP\_/)
	{
		$alias =~ s/^/RefSeq_Prot:/;
		return $alias;
	}
	elsif ($alias =~ s/^eric\|/ERIC\:/) 	{return $alias}
	elsif ($alias =~ s/^gi\|/NCBI_gi:/) 	{return $alias}
	elsif ($alias =~ s/^uni\|/UniProtKB:/) 	{return $alias}
	elsif ($alias =~ s/^kegg\|(.*?)\:/KEGG\:$1\+/i) 	{return $alias}
	elsif ($alias =~ s/^sp\|/Swiss-Prot:/) 	{return $alias}
	elsif ($alias =~ s/^tr\|/TrEMBL:/) 	{return $alias}
	elsif ($alias =~ s/^tigr\|/TIGR_CMR:/) 	{return $alias}
	#
	# 2007-08-15: the validator is saying Locus_Tag isn't acceptable.
	# elsif ($alias =~ s/^LocusTag/Locus_Tag/) {return $alias}
	elsif ($alias =~ s/^img\|/IMG:/) 	{return $alias}
	else
	{
		return undef;
	}
}


=head3 by_all_aliases

usage: @pegs = $fig->by_all_aliases($alias);

This combines by_alias and by_raw_alias to get a single list of pegs for a given alias. NOTE: feature_aliases does both calls,
therefore, this should provide round-robin access with feature_aliases.

=cut

sub by_all_aliases {
	my ($self, $alias) = @_;
	my %pegs = map {$_=>1} ($self->by_alias($alias), $self->by_raw_alias($alias));

	return wantarray() ? keys %pegs : [keys %pegs];
}

=head3 by_alias

usage: $peg = $fig->by_alias($alias)

Returns a FIG id if the alias can be converted.  Right now we convert aliases
of the form NP_* (RefSeq IDs), gi|* (GenBank IDs), sp|* (Swiss Prot), uni|* (UniProt),
kegg|* (KEGG) and maybe a few more

=cut

sub by_alias {
    my($self,$alias,$genome) = @_;
    my($rdbH,$relational_db_response);
    my ($peg, $flag) = FIGRules::NormalizeAlias($alias);
    if ($flag) {
        return $peg;
    } else {
        my $genomeQ = $genome ? quotemeta $genome : "";
        $rdbH = $self->db_handle;

        if (($relational_db_response = $rdbH->SQL("SELECT id FROM ext_alias WHERE ( alias = ? )", undef, $peg)) &&
            (@$relational_db_response > 0)) {

            if (@$relational_db_response == 1) {
                $peg = $relational_db_response->[0]->[0];
                return wantarray() ? ($peg) : $peg;
            } elsif (wantarray()) {
                return map { $_->[0] } @$relational_db_response;
            }
        }
        return wantarray() ? () : "";
    }
}

=head3 by_raw_alias

usage: $peg = $fig->by_raw_alias($alias)

Returns all FIG ids having the given alias.  Unlike by_alias, we do not attempt any
kind of normalization.  I'm not sure this function is needed, but by_alias searches
only in ext_alias table whereas here I'm searching in the features table.  ext_alias
does not have all external aliases which is keeping my code from working.  In particular,
it lacks EnsemblGene.  It would be nice to combine these two functions.  -Ed
=cut

sub by_raw_alias {
    my($self,$alias) = @_;
    my($rdbH,$relational_db_response);
    my ($peg);

    $rdbH = $self->db_handle;

	# RAE: initially this LIKE was %,$alias,% but that doesn't work for first/last aliases in the list
    if (($relational_db_response = $rdbH->SQL("SELECT id FROM features WHERE aliases LIKE \'%$alias%\'")) && (@$relational_db_response > 0)) {

        if (@$relational_db_response == 1) {
            $peg = $relational_db_response->[0]->[0];
            return wantarray() ? ($peg) : $peg;
        } elsif (wantarray()) {
            return map { $_->[0] } @$relational_db_response;
        }
    }
    return wantarray() ? () : "";
}

sub to_alias {
    my($self,$fid,$type) = @_;

    my @aliases = $self->feature_aliases($fid);
    if ($type)
    {
        @aliases = grep { $_ =~ /^$type\|/ } @aliases;
    }

    if (wantarray())
    {
        return @aliases;
    }
    elsif (@aliases > 0)
    {
        return $aliases[0];
    }
    else
    {
        return "";
    }
}

=head3 possibly_truncated

usage: $fig->possibly_truncated($feature_id) or $fig->possibly_truncated($genome, $loc)


Returns the empty string if the feature or location is not near either end of a contig.

Returns 'stop' if the feature or location is on the 'plus' strand and near the end of a contig,
or is on the 'minus' starnd and near the beginning of the contig.

Returns 'start' if the feature or location is on the 'plus' strand and near the beginning of a contig,
or is on the 'minus' starnd and near the end of the contig.

Possibly truncated STOPs have return priority over possibly truncated STARTs.

=cut

sub possibly_truncated {
    my($self, $arg1, $arg2) = @_;
    my($fid,  $loc,  $genome);

    if (($arg1 =~ m/^fig\|\d+\.\d+\.([^\.]+)\.\d+$/) && (not defined($arg2))) {
        $fid    = $arg1;
        $loc    = $self->feature_location($fid);
        $genome = $self->genome_of($fid);
    }
    elsif (($arg1 =~ m/^\d+\.\d+$/) && ($arg2 =~ m/\S+_\d+_\d+/)) {
        $genome = $arg1;
        $loc    = $arg2;
    }
    else {
        confess "Invalid Arguments ", join(", ", @_), " to FIG::possibly_truncated";
    }

    my ($contig, $beg, $end) = $self->boundaries_of($loc);
    if (! (defined($contig) && defined($beg) && defined($end))) { return 0 }

    if    ($self->possibly_truncated_stop($genome, $contig, $beg, $end)) {
	#...Truncated STOP is worse than truncated START...
	return qq(stop);
    }
    elsif ($self->possibly_truncated_start($genome, $contig, $beg, $end)) {
	return qq(start);
    }
    else {
	return qq();
    }
}

sub possibly_truncated_start {
    my ($self, $genome, $contig, $beg, $end) = @_;
    my $strand = ($end <=> $beg);

    if (  (($strand > 0) && ($beg < 300))
       || (($strand < 0) && ($beg > ($self->contig_ln($genome,$contig) - 300)))
       ) {
	return qq(start);
    }
    return qq();
}

sub possibly_truncated_stop {
    my ($self, $genome, $contig, $beg, $end) = @_;
    my $strand = ($end <=> $beg);

    if (  (($strand > 0) && ($end > ($self->contig_ln($genome,$contig) - 300)))
       || (($strand < 0) && ($end < 300))
       ) {
	return qq(stop);
    }
    return qq();
}



=head3 possible_frameshift

=over 4

=item USAGE:

C<< my $fs = $fig->possible_frameshift($peg); >>

=item RETURNS:

A pointer to a list of the form [ContigName,BegOfRegionContaining,EndOfContainingRegion,DNAofContaining,TemplatePEGid]

boolean C<FALSE> otherwise.

=back

=cut

use gjoparseblast;

sub possible_frameshift {
    my($self,$peg) = @_;

    my($tmp_dir) = $FIG_Config::temp;

    my $tmp_dna  = "$tmp_dir/tmp_dna.$$.fasta";
    my $tmp_prot = "$tmp_dir/tmp_prot.$$.fasta";

    #...Skip tests and return '0' if truncated...
    if (! $self->possibly_truncated($peg))
    {
	#...Get best precomputed BLAST hit if E-value < 1.0e-20:
	my @sims = $self->sims($peg,5,1.0e-20,"fig");
	while ((@sims > 0) && $self->possibly_truncated($sims[0]->id2)) { shift @sims }

	#...If a sim was returned:
	if (my $sim = shift @sims)
	{
	    #...Get best hit FID and boundaries:
	    my $peg2 = $sim->id2;
	    my $ln1  = $sim->ln1;
	    my $ln2  = $sim->ln2;
	    my $b2   = $sim->b2;
	    my $e2   = $sim->e2;

	    #...Convert from AA to BP, and pad out w/ 100 bp guard region:
	    my $adjL = 100 + (($b2-1) * 3);
	    my $adjR = 100 + (($ln2 - $e2) * 3);

	    if ($ENV{DEBUG}) { print STDERR "adjL = $adjL adjR = $adjR ln1 = $ln1 peg2 = $peg2 ln2 = $ln2\n" }
	    #...If hit is more than 20% longer than query:

	    if ($ln2 > (1.2 * $ln1))
	    {
		#...Get and parse query location:
		my $loc = $self->feature_location($peg);

		if ($loc =~ /^(\S+)_(\d+)_(\d+)/)
		{
		    my $genome_id = &FIG::genome_of($peg);
		    my $contig = $1;
		    my $beg    = $2;
		    my $end    = $3;

		    #...Extract DNA subsequence, including guard regions:
		    my($begA,$endA,$dna);
		    if ($beg < $end)
		    {
			$begA = &FIG::max(1, $beg - $adjL);
			$endA = &FIG::min($end+$adjR, $self->contig_ln($genome_id,$contig));
			$dna  = $self->dna_seq($genome_id,join("_",($contig,$begA,$endA)));
		    }
		    else
		    {
			$endA = &FIG::max(1, $beg - $adjL);
			$begA = &FIG::min($end+$adjR, $self->contig_ln($genome_id,$contig));
			$dna  = $self->dna_seq($genome_id,join("_",($contig,$begA,$endA)));
		    }

		    if (defined($dna) && (length($dna) > 90))
		    {
			#...Open tmp-file and write FASTA containing DNA subregion to be BLASTed:
			open( TMP, ">$tmp_dna") || die "could not open $tmp_dna";
			print TMP  ">dna\n$dna\n";
			close(TMP);

			#...Fetch its translation, and print to tmp FASTA file for BLASTing:
			my $prot  = $self->get_translation($peg2);

			if (defined($prot) && (length($prot) > 30))
			{
			    open( TMP, ">$tmp_prot") || die "could not open $tmp_prot";
			    print TMP  ">tmp_prot\n$prot\n";
			    close(TMP);

			    #...Build BLAST nucleotide database for extracted DNA region,
			    #   and TBLASTN $peg2 against the DNA:
			    &FIG::run("formatdb -i $tmp_dna -pF");
			    open(BLAST,"blastall -i $tmp_prot -d $tmp_dna -p tblastn -FF -e 1.0e-20 |")
				|| die "could not blast";

			    #...Parse the TBLASTN output; find and sort HSPs by left boundary:
			    my $db_seq_out = &gjoparseblast::next_blast_subject(\*BLAST,1);
			    if ($ENV{DEBUG}) { print STDERR &Dumper(['blast output',$db_seq_out]) }
			    my @hsps       = sort { $a->[0] <=> $b->[0] }
			                      map { [$_->[9], $_->[10], $_->[12], $_->[13]] }
			                      grep { $_->[1] < 1.0e-20 }
			                      @ { $db_seq_out->[6] };

			    #...Extract HSP boundary pairs:
			    my @prot = map { [$_->[0], $_->[1]] } @hsps;
			    my @dna  = map { [$_->[2], $_->[3]] } @hsps;
			    if ($ENV{DEBUG}) { print STDERR &Dumper(\@prot,\@dna) }

			    #...If the "cover" of the HSPs covers more than 90% of $peg2 w gaps < 3 AA,
			    #   and the "cover" of the HPSs cover more than 90% of the extracted DNA
			    #   w/ gaps < 9 bp (but not a multiple of 3), suspect a possible frameshift:
			    if (&covers(\@prot,length($prot),3,0) && &covers(\@dna,3*length($prot),9,1))
			    {
				unlink($tmp_dna,$tmp_prot);
				return [$contig,$begA,$endA,$dna,$peg2];
			    }
			}
		    }
		}
	    }
	}
    }
    unlink($tmp_dna,$tmp_prot);
    return 0;
}

sub covers {
    my($hsps,$ln,$diff,$must_shift) = @_;

    if ($ENV{DEBUG}) { print STDERR &Dumper(['hsps',$hsps,'ln',$ln,'diff',$diff,'must_shift',$must_shift]) }
    my $hsp1 = shift @$hsps;
    my $hsp2;
    my $merged = 0;
    while ($hsp1 && ($hsp2 = shift @$hsps) &&
	   ($must_shift ? &diff_frames($hsp1,$hsp2) : 1) &&
	   ($hsp1 = &merge($hsp1,$hsp2,$diff)))
    {
	$merged = 1;
	if ($ENV{DEBUG}) { print STDERR &Dumper(['merged',$hsp1]) }
    }
    return ($merged && $hsp1 && (($hsp1->[1] - $hsp1->[0]) > (0.9 * $ln)));
}

sub diff_frames {
    my($hsp1,$hsp2) = @_;
    return ((($hsp1->[1]+1) % 3) != ($hsp2->[0] % 3));
}



=head3 merge

Merge two HSPs unless their overlap or separation is too large.

RETURNS: Merged boundaries if merger succeeds, and C<undef> if merger fails.

=cut

sub merge {
    my($hsp1,$hsp2,$diff) = @_;

    my($b1,$e1) = @$hsp1;
    my($b2,$e2) = @$hsp2;
    return (($e2 > $e1) && (($b2-$e1) <= $diff)) ? [$b1,$e2] : undef;
}

sub near_end {
    my($self,$genome,$contig,$x) = @_;

    return (($x < 300) || ($x > ($self->contig_ln($genome,$contig) - 300)));
}

sub is_real_feature {
    my($self,$fid) = @_;
    my($relational_db_response);

    if ($self->is_deleted_fid($fid)) { return 0 }

    my $rdbH = $self->db_handle;
    return (($relational_db_response = $rdbH->SQL("SELECT id FROM features WHERE ( id = \'$fid\' )")) &&
            (@$relational_db_response == 1)) ? 1 : 0;
}


sub is_real_feature_bulk {
    my( $self, $fids ) = @_;
    return wantarray ? () : []  unless $self && $fids && ref( $fids ) eq 'ARRAY';

    #  Remove deleted fids

    my @fids = grep { ! $self->is_deleted_fid( $_ ) } @$fids;

    #  Try to find feature in DBMS

    my $batch = 1000;
    my $rdbH = $self->db_handle;
    my @real;
    while ( @fids )
    {
	my @ids = splice( @fids, 0, $batch );
	my $in = join( ", ", map { "?" } @ids );
	my $dbms_response = $rdbH->SQL(qq(SELECT id FROM features WHERE id IN ($in)), undef, @ids);
        push @real, map { $_->[0] } @$dbms_response if $dbms_response;
    }

    wantarray ? @real : \@real;
}


=head3 map_peg_to_ids

C<<my $gnum, $pnum = $fig->map_peg_to_ids($peg)>>

Map a peg ID to a pair of numbers describing that peg.

In order to conserve storage and increase performance for some operations (the
functional coupling computation, for instance), we provide a mechanism by which a full peg
(of the form fig|X.Y.peg.Z) is mapped to a pair of integers: a genome number and a PEG
index. We maintain a table genome_mapping that retains the mapping between genome ID
and local genome number. No effort is expended to ensure this mapping is at all coherent
between SEED instances; this is purely a local mechanism for performance enhancement.

=over 4

=item $peg

ID of the peg to be mapped.

=item RETURN

A pair of numbers ($gnum, $pnum)

=back

=cut

sub map_peg_to_ids
{
    my($self, $peg) = @_;

    my $mapperc = $self->cached("_genome_mapper");
    my $mapper = $mapperc->{mapper_obj};
    if (!defined($mapper))
    {
        $mapper = new GenomeIDMap($self);
        $mapperc->{mapper_obj} = $mapper;
    }

    return $mapper->map_peg_to_nums($peg);
}

sub map_ids_to_peg
{
    my($self, @ids) = @_;

    my $mapperc = $self->cached("_genome_mapper");
    my $mapper = $mapperc->{mapper_obj};
    if (!defined($mapper))
    {
        $mapper = new GenomeIDMap($self);
        $mapperc->{mapper_obj} = $mapper;
    }

    return $mapper->map_nums_to_peg(@ids);
}

sub map_genome_to_id
{
    my($self, $genome) = @_;

    my $mapperc = $self->cached("_genome_mapper");
    my $mapper = $mapperc->{mapper_obj};
    if (!defined($mapper))
    {
        $mapper = new GenomeIDMap($self);
        $mapperc->{mapper_obj} = $mapper;
    }

    return $mapper->map_genome_id_to_gnum($genome);
}

sub map_id_to_genome
{
    my($self, $id) = @_;

    my $mapperc = $self->cached("_genome_mapper");
    my $mapper = $mapperc->{mapper_obj};
    if (!defined($mapper))
    {
        $mapper = new GenomeIDMap($self);
        $mapperc->{mapper_obj} = $mapper;
    }

    return $mapper->map_gnum_to_gid($id);
}

################ Routines to process abstract functional coupling for PEGs  ##########################

=head3 abstract_coupled_to

    my @coupled_to = $fig->abstract_coupled_to($peg);

Return a list of functionally coupled PEGs.


=over 4

=item peg

ID of the protein encoding group whose functionally-coupled proteins are desired.

=item RETURN

Returns a list of 4-tuples, each consisting of the ID of a coupled
PEG, a score, a "type" which indicates the method that produced the
score, and "extra data" in the form of a pointer to a list.  If there
are no PEGs functionally coupled to the incoming PEG, it will return
an empty list. If the PEG data is not present, it will return an empty list.

=back

=cut

sub abstract_coupled_to {
    my($self, $peg) = @_;

    my $rdbH = $self->db_handle;
    if (! $rdbH->table_exists('afc')) { return () }

    my $relational_db_response = $rdbH->SQL("SELECT peg2,score,type,extra FROM afc
                                             WHERE  peg1 = \'$peg\' ");
    return sort { ($b->[1] <=> $a->[1]) or ($a->[0] cmp $b->[0]) or ($a->[2] cmp $b->[2]) }
           map { [$_->[0],$_->[1],$_->[2],[split(/\t/,$_->[3])]] }
           @$relational_db_response;
}

################ Routines to process functional coupling for PEGs  ##########################

=head3 coupled_to

    my @coupled_to = $fig->coupled_to($peg);

Return a list of functionally coupled PEGs.

The new form of coupling and evidence computation is based on precomputed data.
The old form took minutes to dynamically compute things when needed.  The old form
still works, if the directory B<Data/CouplingData> is not present.  If it is present,
it theis assumed to contain comprehensive coupling data in the form of precomputed scores
and PCHs.

If B<Data/CouplingData> is present, this routine returns a list of 2-tuples [Peg,Sc].  It
returns the empty list if the peg is not coupled. It returns C<undef> if B<Data/CouplingData>
is not there.

=over 4

=item peg

ID of the protein encoding group whose functionally-coupled proteins are desired.

=item RETURN

Returns a list of 2-tuples, each consisting of the ID of a coupled PEG and a score. If
there are no PEGs functionally coupled to the incoming PEG, it will return an empty
list. If the PEG data is not present, it will return C<undef>.

=back

=cut

use FC;
sub coupled_to {
    my($self, $peg) = @_;
    if ($FIG_Config::use_SAPserver_fc) {
	my $sapO = SAPserver->new();
	my $featureHash =   $sapO->conserved_in_neighborhood({
                                -ids => [$peg]
                            });
        my $l = $featureHash->{$peg};
        my @ret;
        foreach my $pair (@$l) {
		push @ret, [$pair->[1], $pair->[0]];
	}
	return @ret;
    }
    if ($FIG_Config::use_pch_server)
    {
	return $self->net_coupled_to($peg);
    }

     if ($FIG_Config::use_sapling_fc)
     {
	 my $db = $self->sapling;
         return &FC::coupled_to($db,$peg);
     }

    my $rdbH = $self->db_handle;
    if (! $rdbH->table_exists('fc_pegs')) { return undef }

    my $relational_db_response = $rdbH->SQL("SELECT peg2,score FROM fc_pegs
                                             WHERE  peg1 = \'$peg\' ");
    return @$relational_db_response;
}

sub net_coupled_to {
    my($self, $peg) = @_;
    return FIGRules::NetCouplingData('coupled_to', id1 => $peg);
}

sub coupled_to_batch {
    my($self, @peg) = @_;

    return () unless @peg;

    if ($FIG_Config::use_SAPserver_fc) {
	my $sapO = SAPserver->new();
	my $featureHash =   $sapO->conserved_in_neighborhood({
                                -ids => \@peg
                            });
	my @ret;
        foreach my $peg (keys %$featureHash) {
		my $l = $featureHash->{$peg};
		foreach my $pair (@$l) {
			push(@ret, [$peg, $pair->[1], $pair->[0]]);
		}
        }

	return @ret;

    }
    if ($FIG_Config::use_pch_server)
    {
	return $self->net_coupled_to_batch(@peg);
    }

    if ($FIG_Config::use_sapling_fc)
    {
	my $db = $self->sapling;
	return &FC::coupled_to_batch($db,@peg);
     }

    my $rdbH = $self->db_handle;
    if (! $rdbH->table_exists('fc_pegs')) { return undef }

    my $cond = join(", ", map { "'$_'" } @peg);
    my $relational_db_response = $rdbH->SQL("SELECT peg1, peg2,score FROM fc_pegs
                                             WHERE  peg1 in ($cond) ");
    return @$relational_db_response;
}

sub net_coupled_to_batch {
    my($self, @peg) = @_;
    return FIGRules::NetCouplingData('coupled_to_batch', id1 => \@peg);
}

sub net_in_pch_pin_with_and_evidence {
    my($self, $peg) = @_;
    return FIGRules::NetCouplingData('in_pch_pin_with_and_evidence', id1 => $peg);
}

sub net_coupling_evidence
{
    my($self, $peg1, $peg2) = @_;
    return FIGRules::NetCouplingData('coupling_evidence', id1 => $peg1, id2 => $peg2);
}

sub net_coupling_and_evidence {
    my($self, $peg) = @_;
    my @rawList = FIGRules::NetCouplingData('coupling_and_evidence', id1 => $peg);
    # The return is supposed to be a list of 3-tuples, where the third element is
    # another list of 3-tuples. Instead, it comes back as n-tuples. We fix things
    # below.
    my @retVal = ();
    for my $rawTuple (@rawList) {
        my ($score, $p2, @rest) = @{$rawTuple};
        my @ev;
        while (my @x = splice(@rest, 0, 2)) {
            push(@ev, \@x);
        }
        push(@retVal, [$score, $p2, \@ev]);
    }
    return @retVal;
}

sub net_coupling_and_evidence_batch {
    my($self, $pegL) = @_;
    my @rawList = FIGRules::NetCouplingData('coupling_and_evidence_batch', id1 => $pegL);
    # The return is supposed to be a list of 3-tuples, where the third element is
    # another list of 3-tuples. Instead, it comes back as n-tuples. We fix things
    # below.
    my @retVal = ();
    for my $rawTuple (@rawList) {
        my ($peg1, $score, $p2, @rest) = @{$rawTuple};
        my @ev;
        while (my @x = splice(@rest, 0, 2)) {
            push(@ev, \@x);
        }
        push(@retVal, [$peg1, $score, $p2, \@ev]);
    }
    return @retVal;
}

sub net_bbhs {
    my ($self, $peg, $cutoff) = @_;
    return FIGRules::BBHData($peg, $cutoff);
}


=head3 coupling_evidence

usage: @evidence = $fig->coupling_evidence($peg1,$peg2)

The new form of coupling and evidence computation is based on precomputed data.
The old form took minutes to dynamically compute things when needed.  The old form
still works, ikf the directory Data/CouplingData is not present.  If it is present,
it is assumed to contain comprehensive coupling data in the form of precomputed scores
and PCHs.

If Data/CouplingData is present, this routine returns a list of 3-tuples [Peg3,Peg4,Rep].
Here, Peg3 is similar to Peg1, Peg4 is similar to Peg2, and Rep == 1 iff this is a
"representative pair".  That is, we take all pairs and create a representative set
in which each pair is not "too close" to any other pair in the representative set.
Think of "too close" as being roughly 95% identical at the DNA level.  This keeps (usually)
a single pair from a set of different genomes from the same species.

It returns the empty list if the peg is not coupled.  It returns undef, if Data/CouplingData
is not there.

=cut

sub coupling_evidence {
    my($self,$peg1,$peg2) = @_;
    ### Must update to use sapling
    if ($FIG_Config::use_pch_server)
    {
	return $self->net_coupling_evidence($peg1, $peg2);
    }

     if ($FIG_Config::use_sapling_fc)
     {
	 my $db = $self->sapling;
	 return &FC::coupling_evidence($db,$peg1,$peg2);
     }

    my $rdbH = $self->db_handle;
    if (! $rdbH->table_exists('pchs')) { return undef }

    my $relational_db_response = $rdbH->SQL("SELECT peg3,peg4,rep FROM pchs
                                             WHERE  peg1 = \'$peg1\' AND
                                                    peg2 = \'$peg2\'");
    return @$relational_db_response;
}

=head3 coupling_and_evidence

usage: @coupling_data = $fig->coupling_and_evidence($fid,$bound,$sim_cutoff,$coupling_cutoff,$keep_record)

A computation of couplings and evidence starts with a given peg and produces a list of
3-tuples.  Each 3-tuple is of the form

    [Score,CoupledToFID,Evidence]

Evidence is a list of 2-tuples of FIDs that are close in other genomes (producing
a "pair of close homologs" of [$peg,CoupledToFID]).  The maximum score for a single
PCH is 1, but "Score" is the sum of the scores for the entire set of PCHs.

NOTE: once the new version of precomputed coupling is installed (i.e., when Data/CouplingData
is filled with the precomputed relations), the parameters on computing evidence are ignored.

If $keep_record is true, the system records the information, asserting coupling for each
of the pairs in the set of evidence, and asserting a pin from the given $fd through all
of the PCH entries used in forming the score.

=cut

sub coupling_and_evidence {
    my($self,$peg1,$bound,$sim_cutoff,$coupling_cutoff,$keep_record,$try_old) = @_;
    ### Must update to use sapling
    if ($FIG_Config::use_pch_server)
    {
	return $self->net_coupling_and_evidence($peg1);
    }

    my $rdbH = $self->db_handle;
    if ($rdbH->table_exists('fc_pegs') && $self->is_complete(&FIG::genome_of($peg1)))
    {
        my @ans = ();
        my $tuple;
        foreach $tuple ($self->coupled_to($peg1))
        {
            my($peg2,$sc) = @$tuple;
            my $evidence = [map { [$_->[0],$_->[1]] } $self->coupling_evidence($peg1,$peg2)];
            push(@ans,[$sc,$peg2,$evidence]);
        }
        if ((@ans > 0) || (! $try_old))
        {
            return @ans;
        }
    }
    return &coupling_and_evidence1($self,$peg1,$bound,$sim_cutoff,$coupling_cutoff,$keep_record);
}

sub coupling_and_evidence_batch {
    my($self,$peg1L,$bound,$sim_cutoff,$coupling_cutoff,$keep_record,$try_old) = @_;
    ### Must update to use sapling
    if ($FIG_Config::use_pch_server)
    {
	return $self->net_coupling_and_evidence_batch($peg1L);
    }

    my $rdbH = $self->db_handle;

    if ($rdbH->table_exists('fc_pegs'))
    {
        my @ans = ();
        my $tuple;
        foreach $tuple ($self->coupled_to_batch(@$peg1L))
        {
            my($peg1, $peg2,$sc) = @$tuple;
            my $evidence = [map { [$_->[0],$_->[1]] } $self->coupling_evidence($peg1,$peg2)];
            push(@ans,[$sc,$peg2,$evidence]);
        }
        if ((@ans > 0) || (! $try_old))
        {
            return @ans;
        }
    }
    return undef;
}

sub coupling_and_evidence1 {
    my($self,$feature_id,$bound,$sim_cutoff,$coupling_cutoff,$keep_record) = @_;
    my($neighbors,$neigh,$similar1,$similar2,@hits,$sc,$ev,$genome1);

    if ($self->is_deleted_fid($feature_id)) { return () }

    if ($feature_id =~ /^fig\|(\d+\.\d+)/)
    {
        $genome1 = $1;
    }
    else
    {
        return ();
    }
    my $locations = $self->feature_location($feature_id);
    my($contig,$beg,$end) = $self->boundaries_of($locations);
    if (! $contig) { return () }

    ($neighbors,undef,undef) = $self->genes_in_region($self->genome_of($feature_id),
                                                      $contig,
                                                      &min($beg,$end) - $bound,
                                                      &max($beg,$end) + $bound);
    if (@$neighbors == 0) { return () }
    $similar1 = $self->acceptably_close($feature_id,$sim_cutoff);
    @hits = ();

    foreach $neigh (grep { $_ =~ /peg/ } @$neighbors)
    {
        next if ($neigh eq $feature_id);
        $similar2 = $self->acceptably_close($neigh,$sim_cutoff);
        ($sc,$ev) = $self->coupling_ev($genome1,$similar1,$similar2,$bound);
        if ($sc >= $coupling_cutoff)
        {
            push(@hits,[$sc,$neigh,$ev]);
        }
    }
    if ($keep_record)
    {
        $self->add_chr_clusters_and_pins($feature_id,\@hits);
    }
    return sort { $b->[0] <=> $a->[0] } @hits;
}

sub fast_coupling {
    my($self,$peg,$bound,$coupling_cutoff) = @_;
    my($genome,$genome1,$genome2,$peg1,$peg2,$peg3,%maps,$loc,$loc1,$loc2,$loc3);
    my($pairs,$sc,%ev);

    if ($self->is_deleted_fid($peg)) { return undef }

    my @ans = ();

    $genome = &genome_of($peg);
    foreach $peg1 ($self->in_pch_pin_with($peg)) {
        $peg1 =~ s/,.*$//;
        if ($peg ne $peg1) {
            $genome1 = &genome_of($peg1);
            $maps{$peg}->{$genome1} = $peg1;
        }
    }

    $loc = [$self->boundaries_of(scalar $self->feature_location($peg))];
    foreach $peg1 ($self->in_cluster_with($peg)) {
        if ($peg ne $peg1) {
    #       print STDERR "peg1=$peg1\n";
            $loc1 = [$self->boundaries_of(scalar $self->feature_location($peg1))];
            if (&close_enough($loc,$loc1,$bound)) {
                foreach $peg2 ($self->in_pch_pin_with($peg1)) {
                    $genome2 = &genome_of($peg2);
                    if (($peg3 = $maps{$peg}->{$genome2}) && ($peg2 ne $peg3)) {
                        $loc2 = [$self->boundaries_of(scalar $self->feature_location($peg2))];
                        $loc3 = [$self->boundaries_of(scalar $self->feature_location($peg3))];
                        if (&close_enough($loc2,$loc3,$bound)) {
                            push @{$ev{$peg1}}, [$peg3,$peg2];
                        }
                    }
                }
            }
        }
    }
    foreach $peg1 (keys(%ev)) {
        $pairs = $ev{$peg1};
        my @pegMap = $peg, map { $_->[0] } @$pairs;
        $sc = $self->score(\@pegMap);
        if ($sc >= $coupling_cutoff) {
            push(@ans,[$sc,$peg1]);
        }
    }
    return sort { $b->[0] <=> $a->[0] } @ans;
}

sub score {
    my($self,$pegs) = @_;
    my(@ids);

    if ($self->{_no9s_scoring})
    {
        @ids = map { $self->maps_to_id($_) } grep { ! $self->is_environmental($_) } @$pegs;
    }
    else
    {
        @ids = map { $self->maps_to_id($_) } @$pegs;
    }
    return &score1($self,\@ids) - 1;
}

sub score1 {
    my($self,$pegs) = @_;
    my($sim);

    my($iden_cutoff) = 97;
    my($iden_cutoff_gap) = 100 - $iden_cutoff;

    my($first,@rest) = @$pegs;
    my $count = 1;
    my %hits = map { $_ => 1 } @rest;
    my @ordered = sort { $b->[0] <=> $a->[0] }
                  map { $sim = $_; [$sim->iden,$sim->id2] }
                  grep { $hits{$_->id2} }
                  $self->sims($first,1000,1,"raw");
    my %ordered = map { $_->[1] => 1 } @ordered;
    foreach $_ (@rest)
    {
        if (! $ordered{$_})
        {
            push(@ordered,[0,$_]);
        }
    }

    while ((@ordered > 0) && ($ordered[0]->[0] >= $iden_cutoff))
    {
        shift @ordered ;
    }
    while (@ordered > 0)
    {
        my $start = $ordered[0]->[0];
        $_ = shift @ordered;
        my @sub = ( $_->[1] );
        while ((@ordered > 0) && ($ordered[0]->[0] > ($start-$iden_cutoff_gap)))
        {
            $_ = shift @ordered;
            push(@sub, $_->[1]);
        }

        if (@sub == 1)
        {
            $count++;
        }
        else
        {
            $count += &score1($self,\@sub);
        }
    }
    return $count;
}

=head3 add_chr_clusters_and_pins

usage: $fig->add_chr_clusters_and_pins($peg,$hits)

The system supports retaining data relating to functional coupling.  If a user
computes evidence once and then saves it with this routine, data relating to
both "the pin" and the "clusters" (in all of the organisms supporting the
functional coupling) will be saved.

$hits must be a pointer to a list of 3-tuples of the sort returned by
$fig->coupling_and_evidence.

=cut

sub add_chr_clusters_and_pins {
    my($self,$peg,$hits) = @_;
    my(@clusters,@pins,$x,$sc,$neigh,$pairs,$y,@corr,@orgs,%projection);
    my($genome,$cluster,$pin,$peg2);

    if (@$hits > 0) {
        @clusters = ();
        @pins = ();
        my @pinMap = ($peg, map { $_->[1] } @$hits);
        push @clusters, \@pinMap;
        foreach $x (@$hits) {
            ($sc,$neigh,$pairs) = @$x;
            my @mapped = ($neigh, map { $_->[1] } @$pairs);
            push @pins, \@mapped;
            foreach $y (@$pairs) {
                $peg2 = $y->[0];
                if ($peg2 =~ /^fig\|(\d+\.\d+)/) {
                    $projection{$1}->{$peg2} = 1;
                }
            }
        }
        @corr = ();
        @orgs = keys(%projection);
        if (@orgs > 0) {
            foreach $genome (sort { $a <=> $b } @orgs) {
                push @corr, sort { &FIG::by_fig_id($a,$b) } keys(%{$projection{$genome}});
            }
            push(@pins,[$peg,@corr]);
        }

        foreach $cluster (@clusters) {
            $self->add_chromosomal_cluster($cluster);
        }

        foreach $pin (@pins) {
            $self->add_pch_pin($pin);
        }
    }
}

sub coupling_ev {
    my($self,$genome1,$sim1,$sim2,$bound) = @_;
    my($ev,$sc,$i,$j);

    $ev = [];

    $i = 0;
    $j = 0;
    while (($i < @$sim1) && ($j < @$sim2))
    {
        if ($sim1->[$i]->[0] < $sim2->[$j]->[0])
        {
            $i++;
        }
        elsif ($sim1->[$i]->[0] > $sim2->[$j]->[0])
        {
            $j++;
        }
        else
        {
            $self->accumulate_ev($genome1,$sim1->[$i]->[1],$sim2->[$j]->[1],$bound,$ev);
            $i++;
            $j++;
        }
    }
    my @mapped = map { $_->[0] } @$ev;
    return ($self->score(\@mapped),$ev);
}

sub accumulate_ev {
    my($self,$genome1,$feature_ids1,$feature_ids2,$bound,$ev) = @_;
    my($genome2,@locs1,@locs2,$i,$j,$x);

    if ((@$feature_ids1 == 0) || (@$feature_ids2 == 0)) { return 0 }

    $feature_ids1->[0] =~ /^fig\|(\d+\.\d+)/;
    $genome2 = $1;
    @locs1 = map { $x = $self->feature_location($_); $x ? [$self->boundaries_of($x)] : () } @$feature_ids1;
    @locs2 = map { $x = $self->feature_location($_); $x ? [$self->boundaries_of($x)] : () } @$feature_ids2;

    for ($i=0; ($i < @$feature_ids1); $i++)
    {
        for ($j=0; ($j < @$feature_ids2); $j++)
        {
            if (($feature_ids1->[$i] ne $feature_ids2->[$j]) &&
            &close_enough($locs1[$i],$locs2[$j],$bound))
            {
            push(@$ev,[$feature_ids1->[$i],$feature_ids2->[$j]]);
            }
        }
    }
}

sub close_enough {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($locs1,$locs2,$bound) = @_;

#   print STDERR &Dumper(["close enough",$locs1,$locs2]);
    return (($locs1->[0] eq $locs2->[0]) && (abs((($locs1->[1]+$locs1->[2])/2) - (($locs2->[1]+$locs2->[2])/2)) <= $bound));
}

sub acceptably_close {
    my($self,$feature_id,$sim_cutoff) = @_;
    my(%by_org,$id2,$genome,$sim);

    my($ans) = [];

    foreach $sim ($self->sims($feature_id,1000,$sim_cutoff,"fig"))
    {
        $id2 = $sim->id2;
        if ($id2 =~ /^fig\|(\d+\.\d+)/)
        {
            my $genome = $1;
            if (! $self->is_eukaryotic($genome))
            {
                push(@{$by_org{$genome}},$id2);
            }
        }
    }
    foreach $genome (sort { $a <=> $b } keys(%by_org))
    {
        push(@$ans,[$genome,$by_org{$genome}]);
    }
    return $ans;
}

################ Translations of PEGsand External Protein Sequences  ##########################


=head3 translatable

usage: $fig->translatable($prot_id)

The system takes any number of sources of protein sequences as input (and builds an nr
for the purpose of computing similarities).  For each of these input fasta files, it saves
(in the DB) a filename, seek address and length so that it can go get the translation if
needed.  This routine simply returns true iff info on the translation exists.

=cut

sub translatable {
    my($self,$prot) = @_;

    return &translation_length($self,$prot) ? 1 : 0;
}


=head3 translation_length

usage: $len = $fig->translation_length($prot_id)

The system takes any number of sources of protein sequences as input (and builds an nr
for the purpose of computing similarities).  For each of these input fasta files, it saves
(in the DB) a filename, seek address and length so that it can go get the translation if
needed.  This routine returns the length of a translation.  This does not require actually
retrieving the translation.

=cut

sub translation_length {
    my($self,$prot) = @_;

    if ($self->is_deleted_fid($prot)) { return undef }

    $prot =~ s/^([^\|]+\|[^\|]+)\|.*$/$1/;
    my $rdbH = $self->db_handle;
    my $relational_db_response = $rdbH->SQL("SELECT slen,seek FROM protein_sequence_seeks
                                             WHERE  id = \'$prot\' ");

    my @vals = sort { $b->[1] <=> $a->[1] } @$relational_db_response;
    return (@vals > 0) ? $vals[0]->[0] : undef;
}


=head3 get_translation

    my $translation = $fig->get_translation($prot_id);

The system takes any number of sources of protein sequences as input (and builds an nr
for the purpose of computing similarities).  For each of these input fasta files, it saves
(in the DB) a filename, seek address and length so that it can go get the translation if
needed.  This routine returns the stored protein sequence of the specified PEG feature.

=over 4

=item prot_id

ID of the feature (PEG) whose translation is desired.

=item RETURN

Returns the protein sequence string for the specified feature.

=back

=cut
#: Return Type $;
sub get_translation {
    my($self,$id) = @_;
    my($rdbH,$relational_db_response,$fileN,$file,$fh,$seek,$ln,$tran);

    if ($self->is_deleted_fid($id)) { return '' }

    if ($self->{memcache})
    {
	my $trans = $self->{memcache}->get("tr:$id");
	if (defined($trans))
	{
	    # print STDERR "found translation for $id\n";
	    return $trans;
	}
    }

    $rdbH = $self->db_handle;
    my $orig_id = $id;
    $id =~ s/^([^\|]+\|[^\|]+)\|.*$/$1/;

    $relational_db_response = $rdbH->SQL("SELECT fileno, seek, len FROM protein_sequence_seeks WHERE  id = \'$id\' ");

    if ((! ($relational_db_response && @$relational_db_response > 0)) &&
        ($id !~ /^fig\|/) &&
        ($id = $self->by_alias($id)))
    {
        $relational_db_response = $rdbH->SQL("SELECT fileno, seek, len FROM protein_sequence_seeks WHERE  id = \'$id\' ");
    }

    if ($relational_db_response && @$relational_db_response > 0)
    {
        my @vals = sort { $b->[1] <=> $a->[1] } @$relational_db_response;
        ($fileN,$seek,$ln) = @{$vals[0]};
        if (($fh = $self->openF($self->N2file($fileN))) &&
             ($ln > 10))
        {
            seek($fh,$seek,0);
            read($fh,$tran,$ln-1);
            $tran =~ s/\s//g;
	    $self->{memcache}->set("tr:$id", $tran) if $self->{memcache};
            return $tran;
        }
    }

    #
    # If it is a xxx identifier, try finding it in the indexed NR.
    #
    if ($orig_id =~ /^xxx/)
    {
	my $nr = "$FIG_Config::global/nr-std";
	if (-f "$nr.pal" or -f "$nr.phr")
	{
	    my $pipe;
	    if (!open($pipe, "-|", "$FIG_Config::ext_bin/fastacmd", "-d", $nr, "-s", "lcl|$orig_id"))
	    {
		return '';
	    }
	    my $h = <$pipe>;
	    local $/;
	    undef $/;
	    my $ret = <$pipe>;
	    close($pipe);
	    $ret =~ s/\s*//g;
	    $self->{memcache}->set("tr:$id", $ret) if $self->{memcache};
	    return $ret;
	}
    }

    return '';
}

sub get_translation_bulk {
    my($self,$ids) = @_;
    my($rdbH,$relational_db_response,$fileN,$file,$fh,$seek,$ln,$tran);

    my %out;
    my %need;
    $need{$_} = 1 foreach @$ids;
    if ($self->{memcache})
    {
	my $mcout = $self->{memcache}->get_multi(map { "tr:$_" } @$ids);
	map { my $k = $_; s/^tr://; $out{$_} = $mcout->{$k}; delete $need{$_} } keys %$mcout;
    }

    $rdbH = $self->db_handle;

    if (keys %need)
    {
	my $qs = join(", ", map { "?" } keys %need);
	
	my $qry = "SELECT s.id, ft.file, s.seek, s.len " .
	    "FROM protein_sequence_seeks s JOIN file_table ft ON s.fileno = ft.fileno " .
		"WHERE id IN ($qs)";
	# print STDERR "$qry\n";
	
	my $res = $rdbH->SQL($qry, undef, keys %need);
	# print STDERR Dumper(\%need, $res, $qs);
	
	my @update;
	for my $ent (@$res)
	{
	    my($id, $file, $seek, $ln) = @$ent;
	    my $fname = File::Spec->rel2abs($file, $FIG_Config::fig_disk);
	    my $fh = $self->openF($fname);

	    if ($fh)
	    {
		my $tran;
		seek($fh,$seek,0);
		read($fh,$tran,$ln-1);
		$tran =~ s/\s//g;
		$out{$id} = $tran;
		push(@update, ["tr:$id", $tran]) if $self->{memcache};
	    }
	    else
	    {
		warn "Cannot open $fname: $!";
	    }
	}
	if (@update)
	{
	    my $res = $self->{memcache}->set_multi(@update);
	    # print STDERR "get_translation_bulk update returns " . Dumper($res) . "\n";
	}
    }
    return \%out;
}

sub recast_ids {
    my($self,$pat,$ids) = @_;

    my($id,@to,%to_return,$x);
    foreach $id (@$ids)
    {
	@to = map { $_->[0] } $self->mapped_prot_ids($id);
	foreach $x (@to,$id)
	{
	    if ($x =~ /$pat/)
	    {
		$to_return{$x} = 1;
	    }
	}
    }
    return sort keys(%to_return);
}

=head3 mapped_prot_ids

usage: @mapped = $fig->mapped_prot_ids($prot)

This routine is at the heart of maintaining synonyms for protein sequences.  The system
determines which protein sequences are "essentially the same".  These may differ in length
(presumably due to miscalled starts), but the tails are identical (and the heads are not "too" extended).
Anyway, the set of synonyms is returned as a list of 2-tuples [Id,length] sorted
by length.

=cut

sub mapped_prot_ids {
    my($self,$id) = @_;
    my $rdbH = $self->db_handle;
    my $dbh = $rdbH->{_dbh};

    if ($id =~ /^fig\|/ && $self->is_deleted_fid($id)) { return () }

    #
    # Manage cached statement handles to accelerate multiple queries into the db.
    #
    my $query_cache = $self->cached("_mapped_prot_ids_cache");
    if (not exists($query_cache->{q1}))
    {
        $query_cache->{q1} = $dbh->prepare(qq(SELECT maps_to
                                              FROM peg_synonyms
                                              WHERE syn_id = ?));
    }
    if (not exists($query_cache->{q2}))
    {
        #
        # Select distinct to work around the duplicate-rows bug in build_nr.
	#
	# 2015-02-06 - don't do that, it incurs a performance penalty. Work around
	# it with %seen below; it is not even clear at this point if that is needed.
        #
        $query_cache->{q2} = $dbh->prepare(qq(SELECT syn_id,syn_ln,maps_to_ln
                                              FROM peg_synonyms
                                              WHERE maps_to = ?));
    }

    my $q1_sth = $query_cache->{q1};
    my $q2_sth = $query_cache->{q2};

    #
    # Determine the principal synonym for $id.
    #

    $q1_sth->execute($id);
    my $relational_db_response = $q1_sth->fetchall_arrayref();
#    my $relational_db_response = $rdbH->SQL("SELECT maps_to FROM peg_synonyms WHERE  syn_id = \'$id\' ");

    if ($relational_db_response && (@$relational_db_response))
    {
        $id = $relational_db_response->[0]->[0];
        #
        # if we have more than one, we have duplicate lines. Warn and let it still work.
        #
        if (@$relational_db_response > 1)
        {
            warn "Duplicates found in peg_synonyms for syn_id $id\n";
        }
    }

    #
    # Retrieve the list of synonyms for the principal synonym.
    #

    $q2_sth->execute($id);
    $relational_db_response = $q2_sth->fetchall_arrayref();

#    $relational_db_response = $rdbH->SQL("SELECT syn_id,syn_ln,maps_to_ln FROM peg_synonyms WHERE maps_to = \'$id\' ");
    my @good = ();   # we need to filter out deleted fids
    if ($relational_db_response && (@$relational_db_response > 0))
    {
	my %seen;
        foreach my $tuple (@$relational_db_response)
        {
	    my($syn_id, $syn_ln, $maps_to_ln) = @$tuple;
	    next if $seen{$syn_id, $syn_ln, $maps_to_ln};
	    $seen{$syn_id, $syn_ln, $maps_to_ln} = 1;
            if (($syn_id !~ /^fig\|/) || (! $self->is_deleted_fid($syn_id)))
            {
                push(@good,$tuple);
            }
        }
    }

    if ($relational_db_response && (@good > 0))
    {
        return ([$id,$good[0]->[2]],map { [$_->[0],$_->[1]] } @good);
    }
    else
    {
        #
        # If the sequence is a singleton, return it as such.
        #

        my $len = $self->translation_length($id);
        if ($len)
        {
            return ([$id,$len]);
        }
        else
        {
            return ();
        }
    }
}

sub maps_to_id {
    my($self,$id) = @_;

    my $rdbH = $self->db_handle;
    my $relational_db_response = $rdbH->SQL("SELECT maps_to FROM peg_synonyms WHERE  syn_id = \'$id\' ");
    return ($relational_db_response && (@$relational_db_response == 1)) ? $relational_db_response->[0]->[0] : $id;
}

#
# ID correspondence table utilities
#

=head3 get_corresponding_ids

    my @id_list = $fig->get_corresponding_ids($id, $with_type_info);

Return a list of the identifiers that correspond to the given identifier, based on
the  PIR id correspondence table.

=over 4

=item id

Identifer to look up.

=item with_type_info

Pass a true value here to return tuples [id, source-type, link-information] instead of identifiers.

=item RETURN

A list of identifiers if $with_type_info not true; a  list of tuples
[id, source-type, link-information] otherwise.

=back

=cut

sub get_corresponding_ids
{
    my($self, $id, $with_type_info) = @_;

    my $dbf = $self->db_handle();
    my $dbh = $dbf->{_dbh};

    if ($with_type_info)
    {
	my $res = $dbf->SQL(qq(SELECT i2.protein_id, t2.name
			       FROM id_correspondence i1 JOIN id_correspondence i2 ON i1.file_num = i2.file_num AND i1.set_id = i2.set_id
			       	JOIN id_correspondence_type t2 ON t2.id = i2.type
			       	JOIN id_correspondence_type t1 ON t1.id = i1.type
			       WHERE i1.protein_id = ? AND t1.searchable = 1
			       ), undef, $id);
	return @$res;
    }
    else
    {
	my $res = $dbh->selectcol_arrayref(qq(SELECT i2.protein_id
					      FROM id_correspondence i1 JOIN id_correspondence i2 ON i1.file_num = i2.file_num AND i1.set_id = i2.set_id
					      JOIN id_correspondence_type t1 ON t1.id = i1.type
					      WHERE i1.protein_id = ? AND t1.searchable = 1), undef, $id);
	return @$res;
    }
}

################ GFF3 utilities  ##########################

sub get_gff_writer
{
    my($self, %options) = @_;

    my $w = GFFWriter->new($self, %options);

    return $w;
}


################ Assignments of Function to PEGs  ##########################

# set to undef to unset user
#
sub set_user {
    my($self,$user) = @_;
    $self->{_user} = $user;
}

sub get_user {
    my($self) = @_;

    return $self->{_user};
}

=head3 function_of

    my $function = $fig->function_of($id, $user);

or

    my @functions = $fig->function_of($id);

In a scalar context, returns the most recently-determined functional
assignment of a specified feature by a particular user. In a list
context, returns a list of 2-tuples, each consisting of a user ID
followed by a functional assighment by that user. In this case,
the list contains all the functional assignments for the feature.

=over 4

=item id

ID of the relevant feature.

=item user

ID of the user whose assignment is desired (scalar context only)

=item RETURN

Returns the most recent functional assignment by the given user in scalar
context, and a list of functional assignments in list context. Each assignment
in the list context is a 2-tuple of the form [$user, $assignment].

=back

=cut

# Note that we do not return confidence.  I propose a separate function to get both
# function and confidence
#
sub function_of {
    my($self,$id,$user,$strip_comments) = @_;
    my($relational_db_response,@tmp,$entry,$i);
    my $wantarray = wantarray();
    my $rdbH = $self->db_handle;

    if ($self->is_deleted_fid($id)) { return $wantarray ? () : "" }

    if (($id =~ /^fig\|(\d+\.\d+\.peg\.\d+)/) && ($wantarray || $user))
    {
        if (($relational_db_response = $rdbH->SQL("SELECT made_by,assigned_function FROM assigned_functions WHERE ( prot = ? )", undef, $id)) &&
            (@$relational_db_response >= 1))
        {
            @tmp = sort { $a->[0] cmp $b->[0] } map { $_->[1] =~ s/^\s//; $_->[1] =~ s/(\t\S)?\s*$//; [$_->[0],$_->[1]] } @$relational_db_response;
	    #@tmp = grep { $_->[1] !~/\[SS\]/ } @tmp;
	    if ($strip_comments) { @tmp = map { $_->[1] =~ s/\s*\#.*$//; $_ } @tmp }
            for ($i=0; ($i < @tmp) && ($tmp[$i]->[0] ne "master"); $i++) {}
            if ($i < @tmp)
            {
                $entry = splice(@tmp,$i,1);
                unshift @tmp, ($entry);
            }

            my $val;
            if     ($wantarray)                                         { return @tmp }
            elsif  ($user && ($val  = &extract_by_who(\@tmp,$user)))    { return $val }
            elsif  ($user && ($val  = &extract_by_who(\@tmp,"master"))) { return $val }
            else                                                        { return ""   }
        }
    }
    else
    {
        if (($relational_db_response = $rdbH->SQL("SELECT assigned_function FROM assigned_functions WHERE ( prot = ? AND made_by = 'master' )", undef, $id)) &&
            (@$relational_db_response >= 1))
        {
	    @tmp = @$relational_db_response;
	    #@tmp = grep { $_->[0] !~/\[SS\]/ } @tmp;
	    if ($strip_comments) { @tmp = map { $_->[0] =~ s/\s*\#.*$//; $_ } @tmp }
            $tmp[0]->[0]  =~ s/^\s//; $tmp[0]->[0] =~ s/(\t\S)?\s*$//;
            return $wantarray ? (["master",$tmp[0]->[0]]) : $tmp[0]->[0];
        }
    }

    return $wantarray ? () : "";
}

sub function_of_quick {
    my($self,$id,$user) = @_;

    my $cache = $self->cached('_function_of');
    my $sth = $cache->{sth};
    if (!$sth)
    {
	$sth = $self->db_handle()->{_dbh}->prepare(qq(SELECT assigned_function
						      FROM assigned_functions
						      WHERE prot = ?));
	$cache->{sth} = $sth;
    }

    $sth->execute($id);
    my($fn) = $sth->fetchrow();
    return $fn;
}


=head3 function_of_bulk

    my $functionHash = $fig->function_of_bulk(\@fids, $no_del_check);

Return a hash mapping the specified proteins to their master functional assignments.

=over 4

=item fids

Reference to a list of feature IDs.

=item no_del_check

If TRUE, then deleted features B<will not> be removed from the list. The default
is FALSE, which means deleted feature B<will> be removed from the list.

=item RETURN

Returns a reference to a hash mapping feature IDs to their main functional assignments.

=back

=cut

# sub function_of_bulk {
#     my($self,$id_list, $no_del_check) = @_;
#     my($relational_db_response,@tmp,$entry,$i);
#     my $wantarray = wantarray();
#     my $rdbH = $self->db_handle;

#     my(@ids);

#     if ($no_del_check)
#     {
#         @ids = @$id_list;
#     }
#     else
#     {
#         @ids = grep { not $self->is_deleted_fid($_) } @$id_list;
#     }

#     unless (scalar(@ids)) {
# 	return {};
#     }

#     my $cond = join(" or ", map { "prot = '$_'" } @ids);

#     my $res = $rdbH->SQL(qq(SELECT prot, assigned_function
#                             FROM assigned_functions
#                             WHERE ( ( $cond ) AND made_by = 'master' )));

#     my $out = {};
#     map { $out->{$_->[0]} = $_->[1] } @$res;
#     return($out);
# }

sub function_of_bulk {
    my($self,$id_list, $no_del_check) = @_;
    my($relational_db_response,@tmp,$entry,$i);
    my $wantarray = wantarray();
    my $dbh = $self->db_handle->{_dbh};

    my(@ids);

    if ($no_del_check)
    {
        @ids = @$id_list;
    }
    else
    {
        @ids = grep { not $self->is_deleted_fid($_) } @$id_list;
    }

    unless (scalar(@ids)) {
	return {};
    }

    my $out = {};
    #
    # Search in memcache first.
    #
    if ($self->{memcache})
    {
	my $mcout  = $self->{memcache}->get_multi(map { "f:$_" } @ids);
	# print "memcache get " . Dumper($mcout);
	map { my $k = $_; s/^f://;  $out->{$_} = $mcout->{$k} } keys %$mcout;

	@ids = grep { !$mcout->{"f:$_"} } @ids;
    }

    #
    # Query for any remaining.
    #

    my @update;
#    my $sth = $dbh->prepare("SELECT prot, assigned_function FROM assigned_functions WHERE prot = ? AND made_by = 'master'");

    while (@ids)
    {
	my @batch = splice(@ids, 0, 1000);

	my $prots = join(", ", map { $dbh->quote($_) } @batch);
	my $r = $dbh->selectall_arrayref(qq(SELECT prot, assigned_function
					    FROM assigned_functions
					    WHERE prot IN ($prots)));
	for my $ent (@$r)
	{
	    my($id, $func) = @$ent;
	    $out->{$id} = $func;
	    push(@update, ["f:$id", $func, 12 * 60 * 60]);
	}
    }

    if ($self->{memcache} && @update)
    {
	my $res = $self->{memcache}->set_multi(@update);
	# print STDERR  "Updated with memcached: " . Dumper(\@update, $res);
    }

    return($out);
}

=head3 translated_function_of

usage: $function  = $fig->translated_function_of($peg,$user)

You get just the translated function.

=cut

sub translated_function_of {
    my($self,$id,$user) = @_;

    my $func = $self->function_of($id,$user);
    if ($func)
    {
        $func = $self->translate_function($func);
    }
    return $func;
}


sub extract_by_who {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($xL,$who) = @_;
    my($i);

    for ($i=0; ($i < @$xL) && ($xL->[$i]->[0] ne $who); $i++) {}
    return ($i < @$xL) ? $xL->[$i]->[1] : "";
}


=head3 translate_function

usage: $translated_func = $fig->translate_function($func)

Translates a function based on the function.synonyms table.

=cut

sub translate_function {
    my($self,$function) = @_;

    my ($tran,$from,$to,$line);
    if (! ($tran = $self->{_function_translation}))
    {
        $tran = {};
        if (open(TMP,"<$FIG_Config::global/function.synonyms"))
        {
            while (defined($line = <TMP>))
            {
                chomp $line;
                ($from,$to) = split(/\t/,$line);
                $tran->{$from} = $to;
            }
            close(TMP);
        }
        foreach $from (keys(%$tran))
        {
            $to = $tran->{$from};
            if ($tran->{$to})
            {
                delete $tran->{$from};
            }
        }
        $self->{_function_translation} = $tran;
    }

    while ($to = $tran->{$function})
    {
        $function = $to;
    }
    return $function;
}

sub who_set_function {
    my($self,$peg) = @_;

    my $func = $self->function_of($peg);

    my @annotations = grep { index($_->[3],$func) >= 0 }  sort { $b->[1] <=> $a->[1] } $self->feature_annotations($peg,1);
    return (@annotations > 0) ? $annotations[0]->[2] : "";
}


=head3 assign_function_bulk: Use assign_function( $target, $usage, $function, $options )

*** Use: assign_function( $target, $usage, $function, $options ) ***

usage:

    $status = $fig->assign_function_bulk( $fidL, $user, $func, $anno, $conf, $time, $no_expand )

    ( $n_succ, $n_ignore, $n_fail ) = $fig->assign_function_bulk( $fidL, $user, $func, $anno, $conf, $time, $no_expand )

Assigns a function and annotation to a list of features (@$fidL).
Confidence can (and should be if unusual) included.
Assignments are also logged in the annotation file, as is with the supplied
annotation.
The function is assigned to all identical protein sequences, a behavior that
can be overridden by setting the no_md5_expand argument.

In scalar context, return status is number of successful assignments (if any),
or minus the number ignored, or zero if all assignments failed.

In array context, return values are the number successful assignments, the
number of ignored assignments (that is, the function was already that desired),
and the number of failed assignment attempts.

=cut

sub assign_function_bulk
{
    my ( $self, $fidL, $user, $function, $annotation, $confidence, $timestamp, $no_md5_expand ) = @_;
    $self && $fidL && $user && defined $function
        or return wantarray ? () : 0;

    my $opts = {};
    $opts->{ annotation    } = $annotation  if $annotation;
    $opts->{ confidence    } = $confidence  if $confidence;
    $opts->{ timestamp     } = $timestamp   if $timestamp;
    $opts->{ no_md5_expand } =  1           if $no_md5_expand;

    $self->assign_function( $fidL, $user, $function, $opts );
}


=head3 assign_function

Assign a new function to one or more features.

Usage, new form:

    $return = $fig->assign_function( $target, $user, $function, \%options )
    @return = $fig->assign_function( $target, $user, $function, \%options )

Usage, classic form:

    $return = $fig->assign_function( $target, $user, $function, $confidence, $timestamp, $no_md5_expand )
    @return = $fig->assign_function( $target, $user, $function, $confidence, $timestamp, $no_md5_expand )

=over 4

=item target (required)

An id or a reference to a list of ids. An id can be either a fid or an md5.
When explicitly supplied as a target, an md5 is expanded to all matching fids,
regardless of the settings of the expand_md5/no_expand_md5 options.

=item user (required)

User name to be associated with the assignment(s).

=item function (required)

Function to be assigned to the target feature(s). If this will not change the
current assignment, then the assignment is not made, and no annotation is made
for the particular feature.


=item options (optional, new form of command only)

A reference to a hash of key value pairs that modify the behavior, or supply
optional data.

=item Option keys:

 annotation                => string  # annotation to be added for all features in
                                      #     target list whose function is changed. 

 confidence                => value   # confidence of assignment

 expand_md5                => bool    # force md5 expansion (default)

 expand_solid_rectangle    => bool    # force solid rectangle expansion

 no_expand_md5             => bool    # prevent md5 expansion; disables expand_solid_rectangle

 no_expand_solid_rectangle => bool    # prevent solid rectangle expansion

 return_value              => keyword # see Return values, below

 timestamp                 => int     # time stamp to be used (D = current time)

Although optional, it is highly recommended that an annotation be supplied that
explains the reason for annotation. This only applies to features that are
explicitly defined by the target; assignments that are induced transitively by
the md5 of a target fid, or expansion of a solid rectangle, are annotated with
the identity of the user-supplied fid (or md5), and the basis of the projection.

=item confidence (optional if no later parameters are used, classic form of command only)

Confidence can be (and, if unusual, should be) specified.

=item timestamp (optional if no later parameters are used, classic form of command only)

Integer representation of the time to be associated with the assignment.

=item no_md5_expand (optional, classic form of command only)

Do not expand the a fid to all fids of identical sequence.

=item Return values

Return values can be altered by the return_value option (see below).

Default scalar context return codes:

   $status

 where a $status value

    < 0 is minus the number of requested assignments, when none would change the function.
    = 0 indicates no successful assignments and either no valid features or at least one assignment failed.
    > 0 is the number of successful assignments.

Default array context return values:

   ( $n_changed, $n_failed, $n_moot )

     $n_changed is number of function changes.
     $n_failed  is number of attempted changes that failed.
     $n_moot    is number of requested changes than would not change function.

Explicit return_value option values and the corresponding data returned in scalar
and array contexts are:

 return_value => 'success_list'

   \@fids_changed
    @fids_changed

 return_value => 'fail_list'

   \@fids_failed
    @fids_failed

 return_value => 'all_lists'

   [ \@fids_changed, \@fids_failed, \@fids_moot ]
   ( \@fids_changed, \@fids_failed, \@fids_moot )

 return_value => 'success_count'

   $n_changed

 return_value => 'fail_count'

   $n_failed

 return_value => 'all_counts'

   [ $n_changed, $n_failed, $n_moot ]
   ( $n_changed, $n_failed, $n_moot )

 return_value => 'default'

   $status   (see above)   
   ( $n_changed, $n_failed, $n_moot )

=item Configuration

$FIG_Config::expand_solid_rectangle sets default behavior. To set expansion as default:

   package FIG_Config;
   $expand_solid_rectangle = 1;

=back

=cut

sub assign_function
{
    my ( $self, $target, $user, $function ) = splice @_, 0, 4;
    $self
        or print STDERR "FIG::assign_function() called without FIG object.\n"
           and return ();
    $target
        or print STDERR "FIG::assign_function() called without a target.\n"
           and return ();
    $user
        or print STDERR "FIG::assign_function() called without a user.\n"
           and return ();
    defined $function
        or print STDERR "FIG::assign_function() called without a function.\n"
           and return ();

    #  Ensure that $target is a list.

    $target = [ $target ] if ! ref( $target );

    #  Ensure that this looks like a new format call, with an options hash

    my $opts = {};
    #  Already with options
    if ( @_ && $_[0] && ref($_[0]) eq 'HASH' )
    {
        $opts = shift;
    }
    #  Old format: put arg values into hash
    else
    {
        my ( $confidence, $timestamp, $no_expand_md5 ) = @_;
        $opts->{ confidence }    = $confidence     if defined $confidence;
        $opts->{ timestamp }     = $timestamp      if $timestamp;
        $opts->{ no_expand_md5 } = $no_expand_md5  if $no_expand_md5;
    }

    my $confidence = $opts->{ confidence };
    my $timestamp  = $opts->{ timestamp };
    my $return_val = $opts->{ return_value } || 'default';

    #  List of [ $fid, $annotation ] pairs

    my @fid_anno;
    my @to_expand;

    #
    # We need to partition our target list into those that are propagation-locked
    # (in which case we don't expand), and those that are not and hence are expanded.
    #

    for my $fid (@$target)
    {
	if ($self->is_propagation_locked_fid($fid))
	{
	    push(@fid_anno, [$fid, '']);
	}
	else
	{
	    push(@to_expand, $fid);
	}
    }
    if (@to_expand)
    {
	my @more_fid_anno = $self->expand_assignment_target( \@to_expand, $opts );
	if ( ! @more_fid_anno )
	{
	    my $targ_str = join( ', ', @$target );
	    print STDERR "FIG::assign_function() failed to expand target [$targ_str].\n";
	}
	@more_fid_anno = grep { !$self->is_propagation_locked_fid($_->[0]) } @more_fid_anno;
	push(@fid_anno, @more_fid_anno);
    }

    #  Mark start of batch assignment. We have not yet verified that the
    #  assignment will be changed for a particular fid, but we have verified
    #  that they are real features, so the this annotation will be valid.

    my $min_fids_for_batch_assign_anno = 4;
    if ( @fid_anno >= $min_fids_for_batch_assign_anno )
    {
        my $fid = $fid_anno[0]->[0];
        $self->write_opening_batch_annotation( $fid, $user, $target, $function, $opts );
    }

    my ( @succ, @moot, @fail );
    foreach ( @fid_anno )
    {
        my ( $fid, $annotation ) = @$_;
        my $stat = $self->assign_function_0( $fid, $user, $function, $confidence, $timestamp );

        if ( ! $stat )
        {
            print STDERR "FIG::assign_function_0('$fid', '$user', '$function', ...) failed.\n";
        }
        # Record outcomes
        if    ( ! $stat )   { push @fail, $fid; next }
        elsif ( $stat > 0 ) { push @succ, $fid }
        else                { push @moot, $fid }

        # Add annotation if function was actually changed
        $self->add_annotation( $fid, $user, $annotation ) if length( $annotation) && $stat > 0;
    }

    #  Mark end of match assignment

    if ( @fid_anno >= $min_fids_for_batch_assign_anno )
    {
        my $fid = $fid_anno[0]->[0];
        $self->write_closing_batch_annotation( $fid, $user, $target, $function, $opts );
    }

    if ( wantarray )
    {
        return ( $return_val =~ m/^success_?list/i ) ? @succ
             : ( $return_val =~ m/^fail_?list/i )    ? @fail
             : ( $return_val =~ m/^all_?list/i )     ? ( \@succ, \@fail, \@moot )
             : ( $return_val =~ m/^n_?change/i )     ? scalar @succ
             : ( $return_val =~ m/^n_?fail/i )       ? scalar @fail
             : ( $return_val =~ m/^n_?all/i )        ? ( scalar @succ, scalar @fail, scalar @moot )
             :                                         ( scalar @succ, scalar @fail, scalar @moot );
    }
    else
    {
        return ( $return_val =~ m/^success_?list/i ) ? \@succ
             : ( $return_val =~ m/^fail_?list/i )    ? \@fail
             : ( $return_val =~ m/^all_?list/i )     ? [ \@succ, \@fail, \@moot ]
             : ( $return_val =~ m/^n_?change/i )     ? scalar @succ
             : ( $return_val =~ m/^n_?fail/i )       ? scalar @fail
             : ( $return_val =~ m/^n_?all/i )        ? [ scalar @succ, scalar @fail, scalar @moot ]
             :                                         ( @succ ? scalar @succ : @fail ? 0 : -(scalar @moot) );
    }
}


#
#   @fid_annotation = $fig->expand_assignment_target( $target, $options )
#  \@fid_annotation = $fig->expand_assignment_target( $target, $options )
#
#       @fid_annotation = ( [ $fid1, $anno1 ], [ $fid2, $anno2 ], ... )
#
#       $target is an id or a reference to a list of ids
#       an id is a fid or an md5
#
#  Options:
#
#     annotation                  # annotation to be added for all assignments
#                                 #     explicitly in target list
#     expand_md5                  # force md5 expansion (current default)
#     no_expand_md5               # prevent md5 expansion; also disables
#                                 #     expand_solid_rectangle, which requires md5s
#
#     expand_solid_rectangle      # force solid rectangle expansion
#     no_expand_solid_rectangle   # prevent solid rectangle expansion
#
#  $FIG_Config::expand_solid_rectangle sets default. To set expansion as default:
#
#   package FIG_Config;
#   $expand_solid_rectangle = 1;
#
sub expand_assignment_target
{
    my ( $self, $target, $opts ) = @_;
    $self && $target or return wantarray ? () : [];
    $opts = {} unless $opts && ref( $opts ) eq 'HASH';

    my $anno                   = $opts->{ annotation } || '';

    my $expand_md5             = ! $opts->{ no_expand_md5 };                 # default is expand md5

    my $expand_solid_rectangle = $expand_md5                                 # cannot do it without md5
                              && ( $opts->{ expand_solid_rectangle }         # explicit request to expand
                                || ( $FIG_Config::expand_solid_rectangle     # default from $FIG_Config
                                  && ! $opts->{ no_expand_solid_rectangle }  # unless explicit veto
                                   )
                                 );

    #  Build the sets of explicit fid assignments, and md5 values.
    #  The md5 data includes the md5, its annotation, and the fid that
    #  it came from (so that the expansions can be precisely annotated).

    my %fid2anno;  # Feature ids and their associated annotation
    my %md5data;   # md5s that need to be expanded; values are [ md5, anno, rootID ]
    my %done;      # remove duplicates (other hashes are not comprehensive lists)

    $target = [ $target ] if ! ref( $target );
    foreach my $id ( @$target )
    {
        next if ( ! $id ) || $done{ $id };
        $done{ $id } = 1;

        #  If this is an md5, not a fid, add it to the set of md5s to be
        #  processed. This overwrites an md5 that was induced by a previous
        #  fid in the target, as an explicitly listed md5 in the target will
        #  get a different annotation.

        if ( $id =~ /^[0-9A-Za-z]{32}$/ )
        {
            $md5data{ $id } = [ $id, $anno, "md5 $id" ];
            next;
        } 

        #  This is all there is to adding a fid, making sure that it is real.
        #  Regardless, we will try to expand it as an md5 (below).

        if ( $self->is_real_feature( $id ) )
        {
            $fid2anno{ $id } = $anno ;
        }
        else
        {
            print STDERR "FIG::expand_assignment_target: id '$id' failed is_real_feature() test\n";
        }

        #  Expansion by md5?

        if ( $expand_md5 )
        {
            my $md5 = $self->md5_of_peg( $id ) or next;
            #  Do not clobber an existing datum, it might have a stronger "reason".
            $md5data{ $md5 } ||= [ $md5, "Assignment projected from $id based on identical sequence", $id ];
        }
    }

    #  Expand the md5s by solid rectangles?

    if ( $expand_solid_rectangle )
    {
        #  project_md5s_by_solid_rectangles() takes and produces lists, not
        #  hashes, so this might look a little convoluted. If you think about
        #  it, a list is the natural exchange. The hashes are just an index
        #  into the list to avoid duplicates. It returns all existing values,
        #  and any new ones that it finds (which is why it overwrites the list).

        %md5data = map { $_->[0] => $_ }
                   $self->project_md5s_by_solid_rectangles( values %md5data );
    }

    #  Expand the md5s to fids:

    foreach ( values %md5data )
    {
        my ( $md5, $anno, undef ) = @$_;

        #  Other fids with identical sequence. pegs_with_md5() does not
        #  currently check for is_real_feature(), so we do it here.

        my @fids = grep { ( ! $done{ $_ }++ ) && $self->is_real_feature( $_ ) }
                   $self->pegs_with_md5( $md5 );

        foreach my $fid ( @fids ) { $fid2anno{ $fid } ||= $anno }
    }

    #  Return list of [ $fid, $annotation ] pairs

    my @fid_anno = map { [ $_, $fid2anno{$_} ] } keys %fid2anno;

    wantarray ? @fid_anno : \@fid_anno;
}


#
#   @md5_anno_rootID = $fig->project_md5s_by_solid_rectangles( @md5_anno_rootID )
#  \@md5_anno_rootID = $fig->project_md5s_by_solid_rectangles( @md5_anno_rootID )
#
#  Beware that each input item is a triple with the md5, the annotation associated
#  with assignments to that md5, and the id from which an assignment is propagated.
#  With this information we can construct informative annotations for the projections.
#  The output list is the same format, and it includes the original parameters.
#
#  We will cache the projection data, in case there are multiple independent calls
#  in a script.
#
sub project_md5s_by_solid_rectangles
{
    my $self = shift;

    my %md5data = map  { $_->[0] => $_ }
                  grep { ref($_) && $_->[0] =~ /^[0-9A-Fa-f]{32}$/ }
                  @_;

    my ( $md5_projections, $md5s_in_projection ) = $self->load_md5_projection_data();

    my %sets;
    foreach my $md5_anno_root ( values %md5data )
    {
        my ( $md5, undef, $root ) = @$md5_anno_root;
        my $anno = "Assignment projected from $root based on solid rectangle";

        my $setlist = $md5_projections->{ $md5 };
        if ( $setlist )
        {
            foreach ( @$setlist )
            {
                $sets{$_} ||= [ $md5, $anno, $root ];
            }
        }
    }

    foreach ( keys %sets )
    {
        my ( undef, $anno, $root ) = @{ $sets{$_} };
        my $md5list = $md5s_in_projection->{ $_ };
        foreach ( $md5list ? @$md5list : () )
        {
            $md5data{$_} ||= [ $_, $anno, $root ];
        }
    }

    wantarray ? values %md5data : [ values %md5data ];
}


#
#  Load hashes for going from md5s to solid rectangle projects sets, and back
#  again.
#
#   ( \%md5_projections, \%md5s_in_projection ) = $fig->load_md5_projection_data();
#
#  The structure of the call is intended to support extracting subsets (by
#  providing a list of relevant md5 values) at some point in the future.
#
sub load_md5_projection_data
{
    my $self = shift;

    my $md5_projections    = $self->cached( '_md5_projections' );
    my $md5s_in_projection = $self->cached( '_md5s_in_projection' );
    if ( ! keys %$md5_projections )
    {
        $md5_projections->{ '' } = [];    # we have tried to load the table
        # read the data
        if ( open( RECT_DATA, "<", "$FIG_Config::global/SolidRectangles.sets") )
        {
            while ( <RECT_DATA> )
            {
                chomp;
                my ( $setid, $md5 ) = split /\t/;
                push @{ $md5_projections->{ $md5 } }, $setid;
                push @{ $md5s_in_projection->{ $setid } }, $md5;
            }
            close( RECT_DATA );
        }
    }

    ( $md5_projections, $md5s_in_projection );
}


sub write_opening_batch_annotation
{
   my( $self, $fid, $user, $target, $function, $opts ) = @_;

   $opts = {} unless $opts;
   $target = [ $target ] unless ref( $target );
   my $fidstr = join( "", map { "\t$_\n" } @$target );
   my $optstr = join( "", map { my $esc = uri_escape_utf8($opts->{$_}); "\t$_:\t$esc\n" } keys %$opts );
   $function = uri_escape_utf8( $function );
   $self->add_annotation( $fid, $user, "batch assignment start for '$function'\n$fidstr\noptions\n$optstr" );
}


sub write_closing_batch_annotation
{
   my( $self, $fid, $user, $target, $function, $opts ) = @_;

   $function = uri_escape_utf8( $function );
   $self->add_annotation( $fid, $user, "batch assignment end for '$function'\n" );
}


#
#  This is now wrapped in a function that expands a peg to all pegs of identical
#  sequence. -- GJO 2012/10/19
#

=head3 assign_function_0

Usage:

    $status = $fig->assign_function_0( $fid, $user, $function, $confidence, $timestamp )

Assigns a new function to a feature. Assignments are also logged in the
annotation file associated with the genome. No action is taken if the request
would not change the current function.

This function is called by assign_function(); it is not intended to be called directly.

=over 4

=item fid

ID of the feature to be changed.

=item user

User name to be associated with the change.

=item function

New function to be assigned to the feature.

=item confidence

Optional. Confidence can (and should be if unusual) included.

=item timestamp

Optional. Integer representation of the time to be associated with the
change.

=item status

Return codes:

  -1  the assignment was not performed because the function would not be changed.
   0  the assignment failed.
   1  the assignment was successful.

=back

=cut

sub assign_function_0 {
    my($self,$fid,$user,$function,$confidence,$timestamp) = @_;
    my($role,$kvs,$kv,$k,$v);

    #  2011-12-14 -- GJO
    #  We need canonical form before testing identity to current:
    $function =~ s/\s+/ /sg;  # No multiple spaces
    $function =~ s/^\s+//;    # No space at begining
    $function =~ s/\s+$//;    # No space at end
    $function =~ s/ ; /; /g;  # No space before semicolon

    my $current_func = $self->function_of($fid,$user,0);
    if ($current_func eq $function) {
        # print STDERR "Not assigning: function already set $function\n";
        return -1;
    }  ### return status of -1 if the function is already current

    $user =~ s/^master://i;
    if (! $self->is_real_feature($fid))
    {
	print STDERR "assign_function: $fid is not real\n";
	return 0;
    }
    if (! $user)
    {
	print STDERR "assign_function: user is not set\n";
	return 0;
    }
    if ($self->is_locked_fid($fid))
    {
	print STDERR "assign_function: $fid is locked\n";
	$self->add_annotation($fid,$user,"attempted to alter assignment, but lock was set");
	return 0;
    }

    my $genome = $self->genome_of($fid);

    my $realuser = $user;  #  For annotation
    $user = 'master';      #  Actual assignments are treated as master assignments

    if ($function =~ /^(.*?)[\!](.*)/)
    {
        ($function,$kvs) = ($1,$2);
        if ($kvs)
        {
            $kvs =~ s/^\s+//;
            $kvs =~ s/\s+$//;
            foreach $kv (split(/\s+[\!\#]\s+/,$kvs))
            {
                if ($kv =~ /^([A-Za-z0-9._\-\+\%]+)\s+\^\s+(.*)$/)
                {
                    ($k,$v) = ($1,$2);
                    if ($v !~ /\S/)
                    {
                        &replace_peg_key_value($self,$fid,$k,"");
                    }
                    else
                    {
                        &replace_peg_key_value($self,$fid,$k,$v);
                    }
                }
                elsif ($kv =~ /^([A-Za-z0-9._\-\+\%]+)$/)
                {
                    &replace_peg_key_value($self,$fid,$1,1);
                }
            }
        }
    }

    $confidence = '' if ! defined $confidence;

    my $uuid = $self->get_uuid;

    $self->broker_log("assign_function", {
	    user => $realuser,
	    genome => $genome,
	    fid => $fid,
	    function => $function,
	    confidence => $confidence,
	    uuid => $uuid,
	});

    my $rdbH = $self->db_handle;

    #
    # Insert new func into the audit trail table first.
    #
    # Ugly repeated code here but it'd be uglier to factor and include conditionals on whether annotation_written
    # field is available.
    #

    if ($self->{have_assignment_auditing})
    {
	$rdbH->SQL(qq(INSERT INTO assigned_functions_log (prot, made_by, assigned_function, quality, org)
		      VALUES (?, ?, ?, ?, ?)), undef,
		   $fid, $realuser, $function, $confidence, $genome);

	#
	# Attempt an update first, since it is most likely to succeed.
	#
	my $rows = $rdbH->SQL(qq(UPDATE assigned_functions
				 SET prot = ?, made_by = ?, assigned_function = ?, quality = ?, org = ?, annotation_written = 'N'
				 WHERE prot = ? AND made_by = ?), undef,
			      $fid, $user, $function, $confidence, $genome,
			      $fid, $user);

	#
	# If it failed (did not update any values), insert a new row.
	#
	if ($rows == 0)
	{
	    my $r2 = $rdbH->SQL(qq(INSERT INTO assigned_functions ( prot, made_by, assigned_function, quality, org, annotation_written )
				   VALUES (?, ?, ?, ?, ?, 'N')), undef,
				$fid, $user, $function, $confidence, $genome);
	}
    }
    else
    {
	#
	# Attempt an update first, since it is most likely to succeed.
	#
	my $rows = $rdbH->SQL(qq(UPDATE assigned_functions
				 SET prot = ?, made_by = ?, assigned_function = ?, quality = ?, org = ?
				 WHERE prot = ? AND made_by = ?), undef,
			      $fid, $user, $function, $confidence, $genome,
			      $fid, $user);

	#
	# If UPDATE failed (did not update any values), insert a new row.
	#
	if ($rows == 0)
	{
	    my $r2 = $rdbH->SQL(qq(INSERT INTO assigned_functions ( prot, made_by, assigned_function, quality, org )
				   VALUES (?, ?, ?, ?, ?)), undef,
				$fid, $user, $function, $confidence, $genome);

	    # print STDERR "Insert returns '$r2'\n";
	}
    }

    $rdbH->SQL("DELETE FROM roles WHERE ( prot = ? AND made_by = ?)", undef,
	       $fid, $user);

    foreach $role (&roles_of_function($function))
    {
        $rdbH->SQL(qq(INSERT INTO roles ( prot, role, made_by, org )
		      VALUES (?, ?, ?, ?)), undef,
		   $fid, $role, $user,  $genome);
    }

    my $file;
    if ( $user eq "master" )
    {
        $file = "$FIG_Config::organisms/$genome/assigned_functions";
    }
    else
    {
        &verify_dir("$FIG_Config::organisms/$genome/UserModels/$user");
        $file = "$FIG_Config::organisms/$genome/UserModels/$user/assigned_functions";
    }

    my $status = 1;

    eval {

	if ( open( TMP, ">>$file" ) )
	{
	    flock(TMP,LOCK_EX) || confess "cannot lock assigned_functions: $!";
	    #  Is there a reason for the seek when the file was openned for append?
	    #  Does flock have a side effect?
	    seek(TMP,0,2)      || confess "failed to seek to the end of the file: $!";
	    print TMP "$fid\t$function\t$confidence\n" || confess "assign_function write failed: $!";
	    close(TMP) || confess "assign_function close failed: $!";
	    chmod(0777,$file);
	}
	else
	{
	    print STDERR "FAILED ASSIGNMENT: $fid\t$function\t$confidence. Could not open $file for append: $!\n";
	    $status = 0;
	    if ($self->{have_assignment_auditing})
	    {
		$rdbH->SQL(qq(INSERT INTO assigned_functions_log (prot, made_by, assigned_function, quality, org)
			      VALUES (?, ?, ?, ?, ?)), undef,
			   $fid, $realuser, "FAILED ASSIGNMENT for $function: could not open $file for append: $!", $confidence, $genome);
	    }
	}
    };
    if ($@)
    {
	if ($self->{have_assignment_auditing})
	{
	    $rdbH->SQL(qq(INSERT INTO assigned_functions_log (prot, made_by, assigned_function, quality, org)
			  VALUES (?, ?, ?, ?, ?)), undef,
		       $fid, $realuser, "FLATFILE WRITE ERRROR: $@\nfor write for function\n$function", $confidence, $genome);
	}
	confess $@;
    }

    #
    # Update the memcached if we are using it.
    #
    if ($self->{memcache})
    {
	$self->{memcache}->set("f:$fid", $function);
    }

    # Check for a Sapling cross-update.
    my $loader = $self->SaplingCheck();
    if ($loader) {
        $loader->UpdateFunction($fid, $function);
    }

    #  We are not getting annotations logged.  So, we will impose it here.

    #
    # if we were passed a valid timestamp, mark the annotation with that
    # timestamp. This is for annotation transfer.
    #
    if ($timestamp && $timestamp !~ /^\d+$/)
    {
	undef $timestamp;
    }
    $self->add_annotation( $fid, $realuser, "Set master function to\n$function\n", $timestamp);

    #print STDERR "Inserting $fid, $realuser, $function\n";
    #
    if ($FIG_Config::function_assigment_trail) {
        my @subs;
        foreach my $x ($self->subsystems_for_peg($fid))
        {
            my ($subsystem,$role) = @$x;
	        push(@subs, $subsystem);
	    }
        my $subsystems = join(":", @subs);

        #print "SUBSYSTEMS $subsystems\n";

        my $res = $rdbH->SQL("SELECT made_by, assigned_function FROM function_trail where prot = ? order by mod_time desc LIMIT 1", undef, $fid);

        my $was_made_by = $res->[0]->[0];
        my $was_assigned_function = $res->[0]->[1];

        #print STDERR "Inserting $fid, $was_made_by, $was_assigned_function, $realuser, $function, $subsystems\n";
        #
        $rdbH->SQL("INSERT INTO function_trail (prot, was_made_by, was_assigned_function, made_by, assigned_function, subsystem) VALUES (?,?,?,?,?,?)", undef, $fid,$was_made_by,$was_assigned_function, $realuser,$function, $subsystems);
    }

    if ($self->{have_assignment_auditing})
    {
	$rdbH->SQL(qq(UPDATE assigned_functions
		      SET annotation_written = 'Y'
		      WHERE prot = ? AND made_by = ?), undef,
		   $fid, $user);
    }

#   print STDERR "realuser = $realuser\n";

    return $status;
}

#######################################################################
# Changes one or more roles in matching assignments / subsystems etc. #
#######################################################################
# If more than one role is supplied, the roles and operators between
# them must all match, in order. So,
#
#    abc; def -> abc @ jkl
#
# will change
#
#    xyz / abc; def / ghi -> xyz / abc @ jkl / ghi
#
# but will not change
#
#    xyz / abc @ def / ghi
#
# or
#
#    xyz / def; abc / ghi
#
# The substitutions are global, so
#
#    abc -> def
#
# will change
#
#    xyz / abc; abc / ghi -> xyz / def; def / ghi
#
# The substitutions are non-overlapping, so
#
#    abc / abc -> def / abc
#
# will change
#
#    abc / abc / abc / abc -> def / abc / def / abc
#
sub change_funcrole {
    my ( $fig, $role, $newname, $seeduser, $synFlag ) = @_;

    my $comment = '';
    my $error   = '';

    if ( ! ( defined($role) && length($role) ) ) {
        $error .= "FIG::change_funcrole() called with bad old role.<BR />\n";
    }
    if ( ! ( defined($newname) && length($newname) ) ) {
        $error .= "FIG::change_funcrole() called with bad new role.<BR />\n";
    }
    if ( ! ( defined($seeduser) && length($seeduser) ) ) {
        $error .= "FIG::change_funcrole() called with bad seeduser.<BR />\n";
    }
    return ( $comment, $error ) if $error;

    #  We need canonical newname before testing identity to current:

    $newname =~ s/^\s+//;    # No space at begining
    $newname =~ s/\s+$//;    # No space at end
    $newname =~ s/\s+/ /g;   # No multiple spaces
    $newname =~ s/ ; /; /g;  # No space before semicolon
    if ( $role eq $newname ) {
        $comment .= "Old name and new name are the same in rename request.<BR />\n";
        return ( $comment, $error );
    }

    ##################################
    # get all pegs with the old role #
    ##################################
    #  Add handling of composite functions.
    #  Make a hash of the distinct roles in old 'role'

    my $composite = $role =~ /\s*;\s+|\s+[\@\/]\s+/;
    my %roles = map { $_ => 1 } split /\s*;\s+|\s+[\@\/]\s+/, $role;

    # Filter candidates for valid features; add support for composite functions.
    # Find candidate pegs by counting the distinct roles covered

    my %pegs;
    foreach ( keys %roles )
    {
        foreach ( $fig->is_real_feature_bulk( [ $fig->seqs_with_role( $_, "master" ) ] ) )
        {
            $pegs{$_}++;
        }
    }

    # for a peg to be a candidate for substitution, it must have all roles

    my $nrequired = keys %roles;
    my @pegs = grep { $pegs{$_} == $nrequired } keys %pegs;
    #print STDERR scalar(@pegs) . " pegs found that might have role.\n";

    ##################################
    # change annotation for each peg #
    ##################################
    #  Build a regular expression for the old role for s/$role_qr/$newfunc/g

    my $role_q  = quotemeta( $role );
    my $role_qr = qr((?:\A|(?<=.;\s|\s[@/]\s))$role_q(?=\s*;\s|\s+[@/]\s|\z));

    # track the changed pegs
    # my $pegcounter = 0;

    my ( %done, %fail );
    my $assign_opt = { return_value => 'all_lists' };
    foreach my $p ( @pegs ) {
        next if $done{ $p };
        my $function = $fig->function_of( $p );

        #  Strip comment (there really should be a space before the ! or #)
        my $rest = '';
	
        # if ( $function =~ /^(.*\S)(\s*[\!\#].*)$/ ) {
        #     $rest     = $2;
        #     $function = $1;
        # }
	($function, $rest) = SeedUtils::strip_func_comment($function);

        # Do a global substitution with lookbehind and lookahead for delimiters

        my $newfunc = $function;
        $newfunc =~ s/$role_qr/$newname/g;

        if ( $newfunc ne $function )
        {
            $assign_opt->{ annotation } = "Role changed from '$function' to '$newfunc'";
            $newfunc .= $rest;

            # now annotate the peg #
            # print STDERR "Assigning for $seeduser to $p: $newfunc\n";
            my ( $succL, $failL, $mootL ) = $fig->assign_function( $p, $seeduser, $newfunc, $assign_opt );
            if ( $succL && @$succL )
            {
                foreach my $p2 ( @$succL )
                {
                    next if $done{ $p2 };
                    $done{ $p2 } = 1;
                    # print STDERR "Annotating role change for $p2.\n";
                    # $fig->add_annotation($p2, $seeduser, "Role changed from '$function' to '$newfunc'");
                    # print STDERR "Annotating new assignment for $p2.\n";
                    # my $time = time;
                    # if ( ! $fig->add_annotation( $p2, $seeduser, "Master function set by $seeduser at $time\n$newfunc\n" ) ) {
                    #     $error .= "Could not set annotation of $p2 to $newfunc<BR />\n";
                    # }
                }
            }
            if ( $failL && @$failL )
            {
                foreach my $p2 ( @$failL )
                {
                    next if $done{ $p2 } || $fail{ $p2 };
                    $fail{ $p2 } = 1;
                    $error .= "Failed in assignment to $p2.<BR />\n";
                }
            }
        }
        elsif ( ! $composite ) {
            $error .= "Could not find $role in $function<BR />\n";
        }
    }

    my $pegcounter = keys %done;
    $comment .= "Changed Role for $pegcounter features based on " . scalar( @pegs ) . " candidates<BR />\n";

    #################################################
    # Collect Model SEED stuff together in the code #
    #################################################
#     my $fba;
#     eval {
#         require ModelSEED::FBAMODELserver;
#         $fba = FBAMODELserver->new();
#     };
#     if ($@) { warn "Error creating fba: $@"; }
    
#     #print STDERR "Model seed object created.\n";

#     if ($pegcounter && $fba) {
#         $fba->changeModelRole({ oldRole => $role, newRole => $newname,
#                                           user => $seeduser, syntaxOnly => $synFlag });
#     }

    ###################################
    # change subsystems for this role #
    ###################################
    #print STDERR "Changing role in subsystems.\n";
    my @subsystems = $fig->function_to_subsystems( $role );

    foreach my $ssn ( @subsystems ) {
        print STDERR "Processing subsystem $ssn.\n";
        # create a new subsystem object #
        my $subsystem = new Subsystem( $ssn, $fig, 0 );
        if ( !defined( $subsystem ) ) {
            print STDERR "Could not get Subsystem Object for $ssn\n";
            next;
        }


        # print STDERR "Change role in $ssn\n";
        my ( $succ, $error ) = $subsystem->change_role( $role, $newname );

        my $thiscomment = '';

        if ( !defined( $error ) ) {
            $thiscomment = "Changed role $role to $newname in subsystem $ssn<BR>\n";
        }

        # here we really edit the files in the subsystem directory #
        $subsystem->incr_version();
        #print STDERR "Updating database.\n";
        $subsystem->db_sync();
        #print STDERR "Writing subsystem.\n";
        $subsystem->write_subsystem();
        #print STDERR "Subsystem $ssn updated.\n";

        #
        # Reload subsystem to ensure consistent state, and resync index.
        #
        my $new_subsystem = new Subsystem( $ssn, $fig, 0 );
        $new_subsystem->db_sync();
        #print STDERR "Subsystem resynced.\n";
        $comment .= $thiscomment;
    }

    ############################
    # print entry into logfile #
    ############################
    my $logfile = "$FIG_Config::data/Logs/functionalroles.rewrite";

    if ( open( LOG, ">>$logfile" ) ) {
        print LOG "Role $role was replaced by $newname\n";
    }
    else {
        $error .= "Logfile could not be opened<BR />\n";
    }

    return ( $comment, $error );
}


sub screwed_up {
    my($self,$peg) = @_;

    my $func = $self->function_of($peg);
    return $func && ($func =~ /(frameshift)|(frame shift)|(interrupt)|(fragment)|(truncate)/);
}

sub hypo {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my $x = (@_ == 1) ? $_[0] : $_[1];

    if (! $x)                             { return 1 }
    if ($x =~ /lmo\d+ protein/i)          { return 1 }
    if ($x =~ /hypoth/i)                  { return 1 }
    if ($x =~ /conserved protein/i)       { return 1 }
    if ($x =~ /gene product/i)            { return 1 }
    if ($x =~ /interpro/i)                { return 1 }
    if ($x =~ /B[sl][lr]\d/i)             { return 1 }
    if ($x =~ /^U\d/)                     { return 1 }
    if ($x =~ /^orf[^_]/i)                { return 1 }
    if ($x =~ /uncharacterized/i)         { return 1 }
    if ($x =~ /pseudogene/i)              { return 1 }
    if ($x =~ /^predicted/i)              { return 1 }
    if ($x =~ /AGR_/)                     { return 1 }
    if ($x =~ /similar to/i)              { return 1 }
    if ($x =~ /similarity/i)              { return 1 }
    if ($x =~ /glimmer/i)                 { return 1 }
    if ($x =~ /unknown/i)                 { return 1 }
    if (($x =~ /domain/i) ||
        ($x =~ /^y[a-z]{2,4}\b/i) ||
        ($x =~ /complete/i) ||
        ($x =~ /ensang/i) ||
        ($x =~ /unnamed/i) ||
        ($x =~ /EG:/) ||
        ($x =~ /orf\d+/i) ||
        ($x =~ /RIKEN/) ||
        ($x =~ /Expressed/i) ||
        ($x =~ /[a-zA-Z]{2,3}\|/) ||
        ($x =~ /predicted by Psort/) ||
        ($x =~ /^bh\d+/i) ||
        ($x =~ /cds_/i) ||
        ($x =~ /^[a-z]{2,3}\d+[^:\+\-0-9]/i) ||
        ($x =~ /similar to/i) ||
        ($x =~ / identi/i) ||
        ($x =~ /ortholog of/i) ||
        (index($x, "Phage protein") == 0) ||
        ($x =~ /structural feature/i))    { return 1 }
    return 0;
}


############################  Similarities ###############################


=head3 nsims

New sims code.

This code takes advantage of a network similarity server if it is available.

We gather sims in the following manner:

    If a local sims directory exists, gather the raw sims for our peg.
    If dynamic sims are available, gather the raw sims from there as well.

    Do an initial pruning of these raw sims, based on the conditions
    passed in to the sims call.

    Locally expand these sims.

    If we are using network sims, retrieve them now, and add to the local sims set.

    Do a final pruning of this set of sims, and sort.

=cut

sub nsims
{
    my ( $self, $id, $maxN, $maxP, $select, $max_expand, $filters ) = @_;

    my $filter_func =  $self->create_sim_filter($maxP, $filters);

    $max_expand = defined( $max_expand ) ? $max_expand : 10000;
    return () if $self->is_deleted_fid( $id );

    my @raw_sims;

    #@raw_sims = $self->get_local_sims($id, $filter_func);


    my %seen = map { $_->[1] => 1 } @raw_sims;

    my @exp_sims;
    if ($select eq 'raw')
    {
        @exp_sims = @raw_sims;
    }
    else
    {
        @exp_sims = $self->expand_local_sims(\@raw_sims, \%seen, $select, $filters);
    }

    #
    # Retrieve network sims if we don't have a sims directory.
    #

    my $want_net_sims = ! -e "$FIG_Config::data/Sims";
    $want_net_sims = 1;

    if ($want_net_sims)
    {
        my @net_sims = $self->get_network_sims($id, \%seen, $maxN, $maxP, $select, $max_expand, $filters);
        push(@exp_sims, @net_sims);
    }

    #
    # If we had no sims for a particular genome in our list,
    # see if we have local sims (from a RAST job) for that genome.
    #

    my %inp_genomes;
    my $ids = ref($id) ? $id : [$id];
    push(@{$inp_genomes{$_->[0]}}, $_->[1]) for map { /fig\|(\d+\.\d+)/ ? [$1, $_] : () } @$ids;

    delete $inp_genomes{$_} for map { $_->id1 =~ /fig\|(\d+\.\d+)/ ? $1 : () } @exp_sims;

    my @need = keys %inp_genomes;

    my($prefix, $file);

    if ($select eq 'raw')
    {
	$file = "similarities";
	$prefix = "_raw_sims";
    }
    else
    {
	$file = "expanded_similarities";
	$prefix = "_sims";
    }

    for my $g (@need)
    {
	my $dir = "$FIG_Config::organisms/$g";
	next unless -f "$dir/$file";
	my $figv = FIGV->new($dir);

	for my $peg (@{$inp_genomes{$g}})
	{
	    my @s = $figv->retrieve_sims($peg, $prefix, $maxP, $select);
	    push(@exp_sims, @s);
	}
    }

    #
    # Do a final filtering for dups.
    #


    #
    # And sort.
    #

    my @sims = $self->sort_sims(\@exp_sims, $filters);

#    print STDERR "Returning for $id: ", Dumper(\@sims);


    return @sims;
}

#
# Create a sim filter-function from the parameters passed.
# Returns true if the sim passed as an argument meets all the requirements.
#
sub create_sim_filter
{
    my($self, $maxP, $filters) = @_;

    my $txt = "sub { \$_ = shift;\n";

    my ( $show_env, $min_sim, $sim_meas, $min_q_cov, $min_s_cov, $sort_by );
    if ( $filters && ref( $filters ) eq "HASH" )
    {
        defined( $filters->{ maxP }      ) and $maxP      = $filters->{ maxP };
        defined( $filters->{ show_env }  ) and $show_env  = $filters->{ show_env };
        defined( $filters->{ min_sim }   ) and $min_sim   = $filters->{ min_sim };
        defined( $filters->{ sim_meas }  ) and $sim_meas  = $filters->{ sim_meas };
        defined( $filters->{ min_q_cov } ) and $min_q_cov = $filters->{ min_q_cov };
        defined( $filters->{ min_s_cov } ) and $min_s_cov = $filters->{ min_s_cov };
        defined( $filters->{ sort_by }   ) and $sort_by   = $filters->{ sort_by };
    }
    defined( $maxP )      or $maxP       =    10;
    defined( $show_env )  or $show_env   =     1;
    defined( $min_sim )   or $min_sim    =     0;
    defined( $sim_meas )  or $sim_meas   =   'id';
    defined( $min_q_cov ) or $min_q_cov  =     0;
    defined( $min_s_cov ) or $min_s_cov  =     0;
    defined( $sort_by )   or $sort_by = 'bits';

    #
    # Initial filter
    #

    $txt .= "return unless \$_->[10] <= $maxP;\n";

    if ($min_sim > 0)
    {
        if ($sim_meas eq 'id')
        {
            $txt .= "return unless \$_->[2] >= $min_sim;\n";
        }
        elsif ($sim_meas eq 'bpp')
        {
            $txt .= "return unless \$_->[2] >= $min_sim * ( \$_->[7] - \$_->[6] + 1); \n";
        }
    }
    #  Query coverage filter

    if ( $min_q_cov > 0 )
    {
        my $thresh = 0.01 * $min_q_cov;
        $txt .= "return unless ( abs( \$_->[7] - \$_->[6] ) + 1 ) >= ( $thresh * \$_->[12] ); \n";
    }

    #  Subject coverage filter

    if ( $min_s_cov > 0 )
    {
        my $thresh = 0.01 * $min_s_cov;
        $txt .= "return unless ( abs( \$_->[9] - \$_->[8] ) + 1 ) >= ( $thresh * \$_->[13] ); \n";
    }

    $txt .= " return 1; }\n";

    #print STDERR "Filter text: $txt\n";

    my $initial_filter = eval $txt;

    return $initial_filter;
}

=head3 osims

usage: @sims = $fig->osims($peg,$maxN,$maxP,$select,$max_expand, $filters)

Returns a list of similarities for $peg such that

    there will be at most $maxN similarities,

    each similarity will have a P-score <= $maxP, and

    $select gives processing instructions:

        "raw" means that the similarities will not be expanded (by far fastest option)
        "fig" means return only similarities to fig genes
        "all" means that you want all the expanded similarities.
        "figx" means exapand until the maximum number of fig sims

By "expanded", we refer to taking a "raw similarity" against an entry in the non-redundant
protein collection, and converting it to a set of similarities (one for each of the
proteins that are essentially identical to the representative in the nr).

Each entry in @sims is a refence to an array. These are the values in each array position:

 0.  The query peg
 1.  The similar peg
 2.  The percent id
 3.  Alignment length
 4.  Mismatches
 5.  Gap openings
 6.  The start of the match in the query peg
 7.  The end of the match in the query peg
 8.  The start of the match in the similar peg
 9.  The end of the match in the similar peg
10.  E value
11.  Bit score
12.  Length of query peg
13.  Length of similar peg
14.  Method

=cut

sub osims {
    my ( $self, $id, $maxN, $maxP, $select, $max_expand, $filters ) = @_;

    my @res;

    if (ref($id) eq 'ARRAY')
    {
	for my $one_id (@$id)
	{
	    push(@res, &osims_one($self, $one_id, $maxN, $maxP, $select, $max_expand, $filters));
	}
	return @res;
    }
    else
    {
	return &osims_one;
    }
}

sub osims_one {
    my ( $self, $id, $maxN, $maxP, $select, $max_expand, $filters ) = @_;
    my( $sim );

    $max_expand = defined( $max_expand ) ? $max_expand : 10000;

    return () if $self->is_deleted_fid( $id );

    #
    # Retrieve the list of synonyms for this peg. The first in the list
    # is the principal synonym.
    #
    my @maps_to = $self->mapped_prot_ids( $id );
    ( @maps_to > 0 ) or return ();

    my $rep_id = $maps_to[0]->[0];
    if ( ! defined( $maps_to[0]->[1] ) )
    {
        print STDERR &Dumper( \@maps_to );
        confess "bad";
    }

    #
    # Find my entry in the list.
    #
    my @entry = grep { $_->[0] eq $id } @maps_to;
    ( @entry == 1 ) and defined( $entry[0]->[1] ) or return ();

    #
    #  Get the similarities. They are based on the principal synonym.
    #

    my @raw_sims = get_raw_sims( $self, $rep_id, $maxP, $filters );

    #  If the query is not the representative, make sims look like it is
    #  by replacing id1 and fixing match coordinates if lengths differ.

    my $delta = $maps_to[0]->[1] - $entry[0]->[1];
    if ( $id ne $rep_id )
    {
        foreach $sim ( @raw_sims )
        {
            $sim->[0]  = $id;
            $sim->[6] -= $delta;
            $sim->[7] -= $delta;
        }
    }

    #  The query must be present for expanding matches to identical sequences.

    if ( ( $max_expand > 0 ) && ( $select ne "raw" ) )
    {
        unshift( @raw_sims, bless( [ $id,
                                     $rep_id,
                                     "100.00",
                                     $entry[0]->[1],
                                     0,
                                     0,
                                     1,        $entry[0]->[1],
                                     $delta+1, $maps_to[0]->[1],
                                     0.0,
                                     2 * $entry[0]->[1],
                                     $entry[0]->[1],
                                     $maps_to[0]->[1],
                                     "blastp"
                                   ], 'Sim'
                                 )
               );
        $max_expand++;
    }

    # print STDERR "\n\n"; for ( @raw_sims ) { print STDERR join( ", ", @{ $_ } ), "\n" }

    #  expand_raw_sims now handles sanity checks on id1 eq id2 and id2
    #  is not deleted.  This lets it keep count of the actual number of
    #  sims reported!

    return expand_raw_sims( $self, \@raw_sims, $maxN, $maxP, $select, 1, $max_expand, $filters );
}


#
# Choose the old sims code.
#

sub sims
{
    my @sims;
    my $which;
    if ($FIG_Config::try_sim_server)
    {
        #
        # Choose the new sims code.
        #
        @sims = &nsims;
        $which = 'new';
    }
    else
    {
        @sims = &osims;
        $which = 'old';
    }
    #open(SIMLOG, ">>$FIG_Config::temp/simlog");
    #print SIMLOG join("\t", $which, @_), "\n";
    #for my $s (@sims)
    #{
        #print SIMLOG join("\t", @$s), "\n";
    #}
    #print SIMLOG "//\n";
    #close(SIMLOG);
    return @sims;
}


sub get_local_sims {
    my ($self, $id, $filter_func) = @_;
    my( $sim );

    return () if $self->is_deleted_fid( $id );

    #
    # Retrieve the list of synonyms for this peg. The first in the list
    # is the principal synonym.
    #
    my @maps_to = $self->mapped_prot_ids( $id );
    ( @maps_to > 0 ) or return ();

    my $rep_id = $maps_to[0]->[0];
    if ( ! defined( $maps_to[0]->[1] ) )
    {
        print STDERR &Dumper( \@maps_to );
        confess "bad";
    }

    #
    # Find my entry in the list.
    #
    my @entry = grep { $_->[0] eq $id } @maps_to;
    ( @entry == 1 ) and defined( $entry[0]->[1] ) or return ();

    #
    #  Get the similarities. They are based on the principal synonym.
    #

    my @raw_sims = get_raw_sims_new( $self, $rep_id, $filter_func);

    
    #  If the query is not the representative, make sims look like it is
    #  by replacing id1 and fixing match coordinates if lengths differ.

    my $delta = $maps_to[0]->[1] - $entry[0]->[1];
    if ( $id ne $rep_id )
    {
        foreach $sim ( @raw_sims )
        {
            $sim->[0]  = $id;
            $sim->[6] -= $delta;
            $sim->[7] -= $delta;
        }
    }

    #  The query must be present for expanding matches to identical sequences.

    unshift( @raw_sims, bless( [ $id,
                                $rep_id,
                                "100.00",
                                $entry[0]->[1],
                                0,
                                0,
                                1,        $entry[0]->[1],
                                $delta+1, $maps_to[0]->[1],
                                0.0,
                                2 * $entry[0]->[1],
                                $entry[0]->[1],
                                $maps_to[0]->[1],
                                "blastp"
                                ], 'Sim'
                             )
           );

    return @raw_sims;
}

sub get_network_sims
{
    my($self, $id, $seen, $maxN, $maxP, $select, $max_expand, $filters) = @_;
    # Get the similarities.
    my $retVal = FIGRules::GetNetworkSims($self, $id, $seen, $maxN, $maxP, $select, $max_expand, $filters);
    # If an error occurred, return an empty list instead of C<undef>.
    if (! defined($retVal)) {
        $retVal = [];
    }
    return @{$retVal};
}

sub sort_sims
{
    my($self, $sims, $filters) = @_;
    my @sorted;

    my $sort_by;
    if ( $filters && ref( $filters ) eq "HASH" )
    {
        defined( $filters->{ sort_by }   ) and $sort_by   = $filters->{ sort_by };
    }
    defined( $sort_by )   or $sort_by = 'bits';

    if    ( $sort_by eq 'id' )                        # Percent identity
    {
        @sorted = sort { $a->[0] cmp $b->[0] or $b->[2] <=> $a->[2] } @$sims;
    }

    elsif ( $sort_by eq 'id2' )                       # Percent identity adjusted
    {
        #  Lower percent identity by 2 standard deviations to prevent random
        #  fluctuation in short sequences from moving them up so often.

        my ( $p, $len, $sigma );
        @sorted = map  { $_->[0] }
                 sort { $a->[0]->[0] cmp $b->[0]->[0] or $b->[1] <=> $a->[1] }
                 map  { $p = 0.01 * $_->[2];                 # fraction identity
                        $len = abs( $_->[7] - $_->[6] ) + 1; # seq len
                        $sigma = sqrt( $p * ( 1 - $p ) / $len ); # binomial sigma
                        [ $_, $_->[2] - 200 * $sigma ]
                      }
                 @$sims;
    }

    elsif ( $sort_by eq 'bpp' )                       # Bits per position
    {
        @sorted = map  { $_->[0] }
                 sort { $a->[0]->[0] cmp $b->[0]->[0] or $b->[1] <=> $a->[1] }
                 map  { [ $_, $_->[11] / abs( $_->[7] - $_->[6] ) ] }
                 @$sims;
    }

    elsif ( $sort_by eq 'bpp2' )                      # Bits per position adjusted
    {
        #  Lower score by 2 standard deviations to prevent random
        #  fluctuation in short sequences from moving them up so often.

        my ( $bpp, $len, $sigma );
        @sorted = map  { $_->[0] }
                 sort { $a->[0]->[0] cmp $b->[0]->[0] or $b->[1] <=> $a->[1] }
                 map  { $len = abs( $_->[7] - $_->[6] ) + 1; # seq len
                        $bpp = $_->[11] / $len;              # bit per pos
                        $sigma = 2.5 * sqrt( 1 / $len );  # simple estimate
                        [ $_, $bpp - 2 * $sigma ]
                      }
                 @$sims;
    }

    else                                              # Bit score (bits)
    {
        @sorted = sort { $a->[0] cmp $b->[0] or $b->[11] <=> $a->[11] } @$sims;
    }

    return @sorted;
}

sub expand_local_sims {
    my( $self, $raw_sims, $seen, $select, $filters) = @_;
    my( $sim, $id1, $id2, %others, $x );

    my $show_env;
    if ( $filters && ref( $filters ) eq "HASH" )
    {
        defined( $filters->{ show_env }   ) and $show_env   = $filters->{ show_env };
    }
    defined( $show_env )   or $show_env   =       1;   # Show environmental by default

    my @sims = ();
    foreach $sim ( @$raw_sims )
    {
        $id2 = $sim->id2;
        $id1 = $sim->id1;

        next if ( $id1 eq $id2 ) || $self->is_deleted_fid( $id2 );

        my @relevant = ();

        #
        # If we are expanding, determine the set of proteins that
        # are equivalent to the protein that we are similar to.
        #
        # Depending on the options passed in, we filter the
        # equivalent proteins found.
        #

        my @maps_to = $self->mapped_prot_ids( $id2 );
        my $ref_len = $maps_to[0]->[1];

        @maps_to = grep { $_->[0] !~ /^xxx\d+/ } grep { $_->[0] !~ /^gnl\|md5/ } @maps_to;

        if ( $select =~ /^figx?$/ )          # Only fig
        {
            @relevant = grep { $_->[0] =~ /^fig/ } @maps_to;
        }
        elsif ( $select =~ /^figx?_?pref/ )  # FIG preferred
        {
            @relevant = grep { $_->[0] =~ /^fig/ } @maps_to;
            #
            # If this id doesn't map to any fig ids, and id2 isn't an xxx id,
            # go ahead and include this sim (and don't bother expanding).
            #
            if ( ! @relevant and $id2 !~ /^xxx\d+$/)
            {
                if (not $seen->{$id2})
                {
                    push @sims, $sim;
                    $seen->{$id2}++;
                }
                next;
            }
        }
        elsif ( $select =~ /^ext/i )         # Not fig
        {
            @relevant = grep { $_->[0] !~ /^fig/ } @maps_to;
        }
        else                                 # All
        {
            @relevant = @maps_to;
        }

        #
        # Include the relevant sims.
        #

        foreach $x ( @relevant )
        {
            my ( $x_id, $x_ln ) = @$x;

            next if $seen->{$x_id};

            $seen->{$x_id} = 1;

            defined( $x_ln ) || confess "x_ln id2='$id2' x_id='$x_id'";
            #next if ( ! $show_env && ( $x_id =~ /^fig\|9999999/ ) );
	    next if ( ! $show_env && ( $self->is_environmental($self->genome_of($x_id)) ) ); # a more inclusive is environmental flag
            next if ( $id1 eq $x_id ) || $self->is_deleted_fid( $x_id );

            defined( $ref_len ) || confess "maps_to";
            my $delta2  = $ref_len - $x_ln;   # Coordinate shift
            my $sim1    = [ @$sim ];                  # Make a copy
            $sim1->[1]  = $x_id;
            $sim1->[8] -= $delta2;
            $sim1->[9] -= $delta2;
            bless( $sim1, "Sim" );
            push( @sims, $sim1 );
        }
    }

    return @sims;
}

sub expand_raw_sims {
    my( $self, $raw_sims, $maxN, $maxP, $select, $dups, $max_expand, $filters ) = @_;
    my( $sim, $id1, $id2, %others, $x );

    #  Set up behavior defaults (pretty wide open):

    my ( $show_env );
    if ( $filters && ref( $filters ) eq "HASH" )
    {
        defined( $filters->{ maxN }       ) and $maxN       = $filters->{ maxN };
        defined( $filters->{ select }     ) and $select     = $filters->{ select };
        defined( $filters->{ max_expand } ) and $max_expand = $filters->{ max_expand };
        defined( $filters->{ show_env }   ) and $show_env   = $filters->{ show_env };
        defined( $filters->{ dups }       ) and $dups       = $filters->{ dups };
    }
    defined( $maxN )       or $maxN       = 1000000;   # Unlimited sims
    defined( $select )     or $select     =    'all';  # Show all expansions
    defined( $max_expand ) or $max_expand =       0;   # But none by default
    defined( $show_env )   or $show_env   =       1;   # Show environmental by default

    $max_expand = 1000000000 if ( $select =~ /^figx/ ); # figx forces unlimited expand

    my @sims = ();
    foreach $sim ( @$raw_sims )
    {
        $id2 = $sim->id2;
        if ( ! $dups )
        {
            next if $others{ $id2 };
            $others{ $id2 } = 1;
        }

        $id1 = $sim->id1;
        if ( ( $select eq "raw" ) || ( $max_expand <= 0 ) )
        {
            #next if ( ! $show_env && ( $id2 =~ /^fig\|9999999/ ) );
	    next if ( ! $show_env && ( $self->is_environmental($self->genome_of($id2)) ) ); # a more inclusive is environmental flag
            next if ( $id1 eq $id2 ) || $self->is_deleted_fid( $id2 );
            push( @sims, $sim );
            return @sims if ( @sims >= $maxN );
        }
        else
        {
            my @relevant = ();
            $max_expand--;

            #
            # If we are expanding, determine the set of proteins that
            # are equivalent to the protein that we are similar to.
            #
            # Depending on the options passed in, we filter the
            # equivalent proteins found.
            #

            my @maps_to = grep { $_->[0] !~ /^xxx\d+/ } grep { $_->[0] !~ /^gnl\|md5/ }  $self->mapped_prot_ids( $id2 );
            if ( $select =~ /^figx?$/ )          # Only fig
            {
                @relevant = grep { $_->[0] =~ /^fig/ } @maps_to;
            }
            elsif ( $select =~ /^figx?_?pref/ )  # FIG preferred
            {
                @relevant = grep { $_->[0] =~ /^fig/ } @maps_to;
                #
                # If this id doesn't map to any fig ids, and id2 isn't an xxx id,
                # go ahead and include this sim.
                #
                if ( ! @relevant and $id2 !~ /^xxx\d+$/)
                {
                    push @sims, $sim;
                    return @sims if ( @sims >= $maxN );
                    next;
                }
            }
            elsif ( $select =~ /^ext/i )         # Not fig
            {
                @relevant = grep { $_->[0] !~ /^fig/ } @maps_to;
            }
            else                                 # All
            {
                @relevant = @maps_to;
            }

            #
            # Include the relevant sims.
            #

            foreach $x ( @relevant )
            {
                my ( $x_id, $x_ln ) = @$x;
                defined( $x_ln ) || confess "x_ln id2='$id2' x_id='$x_id'";
                #next if ( ! $show_env && ( $x_id =~ /^fig\|9999999/ ) );
		next if ( ! $show_env && ( $self->is_environmental($self->genome_of($x_id)) ) ); # a more inclusive is environmental flag
                next if ( $id1 eq $x_id ) || $self->is_deleted_fid( $x_id );

                defined( $maps_to[0]->[1] ) || confess "maps_to";
                my $delta2  = $maps_to[0]->[1] - $x_ln;   # Coordinate shift
                my $sim1    = [ @$sim ];                  # Make a copy
                $sim1->[1]  = $x_id;
                $sim1->[8] -= $delta2;
                $sim1->[9] -= $delta2;
                bless( $sim1, "Sim" );
                push( @sims, $sim1 );
                return @sims if ( @sims >= $maxN );
            }
        }
    }

    return @sims;
}


sub get_raw_sims_new {
    my ( $self, $rep_id, $filter_func) = @_;
    my ( $sim_chunk, $seek, $fileN, $ln, $fh, $file, @lines, $sim );


    my $rdbH = $self->db_handle;
    my $relational_db_response = $rdbH->SQL("SELECT seek, fileN, len FROM sim_seeks WHERE id = \'$rep_id\' ");

    #  Gather all of the acceptable "lines" from the sim chunks

    foreach $sim_chunk ( @$relational_db_response )
    {
        ( $seek, $fileN, $ln ) = @$sim_chunk;
        $file = $self->N2file( $fileN );
        $fh   = $self->openF( $file );
        $fh or confess "could not open sims for $file";

        #  Read file, parse lines, sanity check values, and filter E-value
        #   0.  The query peg
        #   1.  The similar peg
        #   2.  The percent id
        #   3.  Alignment length
        #   4.  Mismatches
        #   5.  Gap openings
        #   6.  The start of the match in the query peg
        #   7.  The end of the match in the query peg
        #   8.  The start of the match in the similar peg
        #   9.  The end of the match in the similar peg
        #  10.  E-value
        #  11.  Bit score
        #  12.  Length of query peg
        #  13.  Length of similar peg
        #  14.  Method

        push @lines, grep { ( @$_ >= 15 ) &&
                            ( $_->[10] =~ /^[0-9.e-]+$/ ) &&  # E-value
                            ( $_->[11] =~ /^[0-9.]+$/ ) &&    # bit score
                            ( $_->[12] =~ /^\d+$/ ) &&        # query len
                            ( $_->[13] =~ /^\d+$/ ) &&        # subj len
                            ( $_->[6]  =~ /^\d+$/ ) &&        # q-match start
                            ( $_->[7]  =~ /^\d+$/ ) &&        # q-match end
                            ( $_->[8]  =~ /^\d+$/ ) &&        # s-match start
                            ( $_->[9]  =~ /^\d+$/ ) &&        # s-match end
                            ( $_->[2]  =~ /^[0-9.]+$/ ) &&    # percent id
                            ( &$filter_func($_) )             # compiled sim filter
                          }
                     map  { [ split( /\t/, $_ ), "blastp" ] }
                     @{ read_block( $fh, $seek, $ln-1 ) };
    }

    push(@lines,     grep { ( @$_ >= 15 ) &&
                            ( $_->[10] =~ /^[0-9.e-]+$/ ) &&  # E-value
                            ( $_->[11] =~ /^[0-9.]+$/ ) &&    # bit score
                            ( $_->[12] =~ /^\d+$/ ) &&        # query len
                            ( $_->[13] =~ /^\d+$/ ) &&        # subj len
                            ( $_->[6]  =~ /^\d+$/ ) &&        # q-match start
                            ( $_->[7]  =~ /^\d+$/ ) &&        # q-match end
                            ( $_->[8]  =~ /^\d+$/ ) &&        # s-match start
                            ( $_->[9]  =~ /^\d+$/ ) &&        # s-match end
                            ( $_->[2]  =~ /^[0-9.]+$/ ) &&    # percent id
                            ( &$filter_func($_) )             # compiled sim filter
                          }
                     &get_dynamic_sims($self,$rep_id));



    #  Bless the raw sims:

    return map { bless( $_, 'Sim' ); $_ } @lines;
}

sub get_raw_sims {
    my ( $self, $rep_id, $maxP, $filters ) = @_;
    my ( $sim_chunk, $seek, $fileN, $ln, $fh, $file, @lines, $sim );

    #  Set up behavior defaults (pretty wide open):

    my ( $show_env, $min_sim, $sim_meas, $min_q_cov, $min_s_cov, $sort_by );
    if ( $filters && ref( $filters ) eq "HASH" )
    {
        defined( $filters->{ maxP }      ) and $maxP      = $filters->{ maxP };
        defined( $filters->{ show_env }  ) and $show_env  = $filters->{ show_env };
        defined( $filters->{ min_sim }   ) and $min_sim   = $filters->{ min_sim };
        defined( $filters->{ sim_meas }  ) and $sim_meas  = $filters->{ sim_meas };
        defined( $filters->{ min_q_cov } ) and $min_q_cov = $filters->{ min_q_cov };
        defined( $filters->{ min_s_cov } ) and $min_s_cov = $filters->{ min_s_cov };
        defined( $filters->{ sort_by }   ) and $sort_by   = $filters->{ sort_by };
    }
    defined( $maxP )      or $maxP       =    10;
    defined( $show_env )  or $show_env   =     1;
    defined( $min_sim )   or $min_sim    =     0;
    defined( $sim_meas )  or $sim_meas   =   'id';
    defined( $min_q_cov ) or $min_q_cov  =     0;
    defined( $min_s_cov ) or $min_s_cov  =     0;
    defined( $sort_by )   or $sort_by = 'bits';

    my $rdbH = $self->db_handle;
    my $relational_db_response = $rdbH->SQL("SELECT seek, fileN, len FROM sim_seeks WHERE id = \'$rep_id\' ");

    #  Gather all of the acceptable "lines" from the sim chunks

    foreach $sim_chunk ( @$relational_db_response )
    {
        ( $seek, $fileN, $ln ) = @$sim_chunk;
        $file = $self->N2file( $fileN );
        $fh   = $self->openF( $file );
        $fh or confess "could not open sims for $file";

        #  Read file, parse lines, sanity check values, and filter E-value
        #   0.  The query peg
        #   1.  The similar peg
        #   2.  The percent id
        #   3.  Alignment length
        #   4.  Mismatches
        #   5.  Gap openings
        #   6.  The start of the match in the query peg
        #   7.  The end of the match in the query peg
        #   8.  The start of the match in the similar peg
        #   9.  The end of the match in the similar peg
        #  10.  E-value
        #  11.  Bit score
        #  12.  Length of query peg
        #  13.  Length of similar peg
        #  14.  Method

	my @raw = @{ read_block( $fh, $seek, $ln-1 ) };
        push @lines, grep { ( @$_ >= 15 ) &&
                            ( $_->[10] =~ /^[0-9.e-]+$/ ) &&  # E-value
                            ( $_->[10] <= $maxP )   &&        # E-value test
                            ( $_->[11] =~ /^[0-9.]+$/ ) &&    # bit score
                            ( $_->[12] =~ /^\d+$/ ) &&        # query len
                            ( $_->[13] =~ /^\d+$/ ) &&        # subj len
                            ( $_->[6]  =~ /^\d+$/ ) &&        # q-match start
                            ( $_->[7]  =~ /^\d+$/ ) &&        # q-match end
                            ( $_->[8]  =~ /^\d+$/ ) &&        # s-match start
                            ( $_->[9]  =~ /^\d+$/ ) &&        # s-match end
                            ( $_->[2]  =~ /^[0-9.]+$/ )       # percent id
                          }
                     map  { [ split( /\t/, $_ ), "blastp" ] }
                     @raw;
    }

    push(@lines,     grep { ( @$_ >= 15 ) &&
                            ( $_->[10] =~ /^[0-9.e-]+$/ ) &&  # E-value
                            ( $_->[10] <= $maxP )   &&        # E-value test
                            ( $_->[11] =~ /^[0-9.]+$/ ) &&    # bit score
                            ( $_->[12] =~ /^\d+$/ ) &&        # query len
                            ( $_->[13] =~ /^\d+$/ ) &&        # subj len
                            ( $_->[6]  =~ /^\d+$/ ) &&        # q-match start
                            ( $_->[7]  =~ /^\d+$/ ) &&        # q-match end
                            ( $_->[8]  =~ /^\d+$/ ) &&        # s-match start
                            ( $_->[9]  =~ /^\d+$/ ) &&        # s-match end
                            ( $_->[2]  =~ /^[0-9.]+$/ )       # percent id
                          }
                     &get_dynamic_sims($self,$rep_id));



    my @linesS = sort { $a->[10] <=> $b->[10] } @lines;  # now sort and remove duplicates
    @lines = ();
    foreach $_ (@linesS)
    {
        if ((@lines == 0) || ($lines[$#lines]->[0] ne $_->[0]) || ($lines[$#lines]->[1] ne $_->[1]))
        {
            push(@lines,$_);
        }
    }

    #  Similarity filter

    if ( $min_sim > 0 )
    {
        if    ( $sim_meas eq 'id'  )
        {
            @lines = grep { $_->[2] >= $min_sim } @lines;
        }
        elsif ( $sim_meas eq 'bpp' )
        {
            @lines = grep { $_->[11] >= $min_sim * ( $_->[7] - $_->[6] + 1 ) } @lines;
        }
    }

    #  Query coverage filter

    if ( $min_q_cov > 0 )
    {
        my $thresh = 0.01 * $min_q_cov;
        @lines = grep { ( abs( $_->[7] - $_->[6] ) + 1 ) >= ( $thresh * $_->[12] ) } @lines;
    }

    #  Subject coverage filter

    if ( $min_s_cov > 0 )
    {
        my $thresh = 0.01 * $min_s_cov;
        @lines = grep { ( abs( $_->[9] - $_->[8] ) + 1 ) >= ( $thresh * $_->[13] ) } @lines;
    }

    #  Order the surviving raw sims by requested criterion:

    if    ( $sort_by eq 'id' )                        # Percent identity
    {
        @lines = sort { $b->[2] <=> $a->[2] } @lines;
    }

    elsif ( $sort_by eq 'id2' )                       # Percent identity adjusted
    {
        #  Lower percent identity by 2 standard deviations to prevent random
        #  fluctuation in short sequences from moving them up so often.

        my ( $p, $len, $sigma );
        @lines = map  { $_->[0] }
                 sort { $b->[1] <=> $a->[1] }
                 map  { $p = 0.01 * $_->[2];                 # fraction identity
                        $len = abs( $_->[7] - $_->[6] ) + 1; # seq len
                        $sigma = sqrt( $p * ( 1 - $p ) / $len ); # binomial sigma
                        [ $_, $_->[2] - 200 * $sigma ]
                      }
                 @lines;
    }

    elsif ( $sort_by eq 'bpp' )                       # Bits per position
    {
        @lines = map  { $_->[0] }
                 sort { $b->[1] <=> $a->[1] }
                 map  { [ $_, $_->[11] / abs( $_->[7] - $_->[6] ) ] }
                 @lines;
    }

    elsif ( $sort_by eq 'bpp2' )                      # Bits per position adjusted
    {
        #  Lower score by 2 standard deviations to prevent random
        #  fluctuation in short sequences from moving them up so often.

        my ( $bpp, $len, $sigma );
        @lines = map  { $_->[0] }
                 sort { $b->[1] <=> $a->[1] }
                 map  { $len = abs( $_->[7] - $_->[6] ) + 1; # seq len
                        $bpp = $_->[11] / $len;              # bit per pos
                        $sigma = 2.5 * sqrt( 1 / $len );  # simple estimate
                        [ $_, $bpp - 2 * $sigma ]
                      }
                 @lines;
    }

    else                                              # Bit score (bits)
    {
        @lines = sort { $b->[11] <=> $a->[11] } @lines;
    }

    #  Bless the raw sims:
    return map { bless( $_, 'Sim' ); $_ } @lines;
}


sub get_dynamic_sims {
    my($self,$prot_id) = @_;
    my $tuples;

    my $rdbH = $self->db_handle;

    if ($rdbH->table_exists('dynamic_sims') &&
        ($tuples = $rdbH->SQL("SELECT id1,id2,iden,ali_ln,mismatches,gap_openings,b1,e1,b2,e2,psc,bit_sc,ln1,ln2 FROM dynamic_sims WHERE id1 = '$prot_id'")) &&
        (@$tuples > 0))
    {
        my @tuples = ();
        foreach $_ (@$tuples)
        {
            push(@$_,"blastp");
            push(@tuples,$_);
        }
        return @tuples;
    }
    return ();
}

sub insert_dynamic_sims {
    my($self,$sims) = @_;
    my($sim);
    my $rdbH = $self->db_handle;

    if (! $rdbH->table_exists('dynamic_sims'))
    {
        $rdbH->create_table( tbl => 'dynamic_sims',
                             flds => 'id1 varchar(32),id2 varchar(32), iden float, ali_ln integer, mismatches float,' .
                                     'gap_openings float, b1 integer, e1 integer, b2 integer, e2 integer, ' .
                                     'psc float, bit_sc float, ln1 integer, ln2 integer');
        $rdbH->create_index( tbl => 'dynamic_sims', idx => 'dynamic_sims_idx_id1', flds => 'id1');
        $rdbH->create_index( tbl => 'dynamic_sims', idx => 'dynamic_sims_idx_id2', flds => 'id2');
    }

    my $rc = 1;
    foreach $sim (@$sims)
    {
        my($id1,$id2,$iden,$ali_ln,$mismatches,$gap_openings,$b1,$e1,$b2,$e2,$psc,$bit_sc,$ln1,$ln2) = @$sim;
        if (! ($rdbH->SQL("INSERT INTO dynamic_sims
                           (id1,id2,iden,ali_ln,mismatches,gap_openings,b1,e1,b2,e2,psc,bit_sc,ln1,ln2)
                           VALUES ('$id1','$id2',$iden,$ali_ln,$mismatches,$gap_openings,$b1,$e1,$b2,$e2,$psc,$bit_sc,$ln1,$ln2)") &&
               $rdbH->SQL("INSERT INTO dynamic_sims
                           (id1,id2,iden,ali_ln,mismatches,gap_openings,b1,e1,b2,e2,psc,bit_sc,ln1,ln2)
                           VALUES ('$id2','$id1',$iden,$ali_ln,$mismatches,$gap_openings,$b2,$e2,$b1,$e1,$psc,$bit_sc,$ln2,$ln1)")))

        {
            $rc = 0;
        }
    }
    return $rc;
}

sub insert_dynamic_sims_file {
    my($self,$sims_file) = @_;
    my($sim);
    my $rdbH = $self->db_handle;

    if (! $rdbH->table_exists('dynamic_sims'))
    {
        $rdbH->create_table( tbl => 'dynamic_sims',
                             flds => 'id1 varchar(32),id2 varchar(32), iden float, ali_ln integer, mismatches float,' .
                                     'gap_openings float, b1 integer, e1 integer, b2 integer, e2 integer, ' .
                                     'psc float, bit_sc float, ln1 integer, ln2 integer');
    }

    #
    # If we're using postgres we can optimize by opening a pipe
    # to a COPY table FROM STDIN
    #
    if ($rdbH->{_dbms} eq "Pg")
    {
        print STDERR "Using pg optimized insert\n";
        $rdbH->drop_index( tbl => 'dynamic_sims', idx => 'dynamic_sims_idx_id1');
        $rdbH->drop_index( tbl => 'dynamic_sims', idx => 'dynamic_sims_idx_id2');
        my $rc= $self->insert_dynamic_sims_pg($sims_file);
        $rdbH->create_index( tbl => 'dynamic_sims', idx => 'dynamic_sims_idx_id1', flds => 'id1');
        $rdbH->create_index( tbl => 'dynamic_sims', idx => 'dynamic_sims_idx_id2', flds => 'id2');
        return $rc;
    }

    my $rc = 1;

    my $sth = $rdbH->{_dbh}->prepare(
                      qq(INSERT INTO dynamic_sims
                         (id1,id2,iden,ali_ln,mismatches,gap_openings,b1,e1,b2,e2,psc,bit_sc,ln1,ln2)
                         VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)));


    if (!open(SIMS, "<$sims_file"))
    {
        warn "Cannot open $sims_file: $!\n";
        return 0;
    }

    while (<SIMS>)
    {
        chomp;
        my($id1,$id2,$iden,$ali_ln,$mismatches,$gap_openings,$b1,$e1,$b2,$e2,$psc,$bit_sc,$ln1,$ln2) = split(/\t/);

        if (!$sth->execute($id1,$id2,$iden,$ali_ln,$mismatches,$gap_openings,$b1,$e1,$b2,$e2,$psc,$bit_sc,$ln1,$ln2))
        {
            warn "SQL error: " . $rdbH->{_dbh}->errstr;
            return 0;
        }
        if (!$sth->execute($id2,$id1,$iden,$ali_ln,$mismatches,$gap_openings,$b2,$e2,$b1,$e1,$psc,$bit_sc,$ln2,$ln1))
        {
            warn "SQL error: " . $rdbH->{_dbh}->errstr;
            return 0;
        }
    }
    return $rc;
}

sub insert_dynamic_sims_pg {
    my($self,$sims_file) = @_;
    my($sim);
    my $rdbH = $self->db_handle;
    my $db = $rdbH->{_dbh};

    $db->do("copy dynamic_sims from stdin");

    open(S, "<$sims_file") or die "Cannot open sims $sims_file: $!\n";

    my $num_per_copy = 5000;
    my $count = 0;
    while (<S>)
    {
        chomp;
        my($id1,$id2,$iden,$ali_ln,$mismatches,$gap_openings,$b1,$e1,$b2,$e2,$psc,$bit_sc,$ln1,$ln2) = split(/\t/);

        $db->func(join("\t", $id1,$id2,$iden,$ali_ln,$mismatches,
                       $gap_openings,$b1,$e1,$b2,$e2,$psc,$bit_sc,$ln1,$ln2) . "\n", 'putline');
        $db->func(join("\t", $id2,$id1,$iden,$ali_ln,$mismatches,
                       $gap_openings,$b2,$e2,$b1,$e1,$psc,$bit_sc,$ln2,$ln1) . "\n", 'putline');

        if ($count++ >= $num_per_copy)
        {
            $db->func("\\.\n", 'putline');
            $db->func("endcopy");
            print "Write $.\n";
            $db->do("copy dynamic_sims from stdin");
            $count = 0;
        }
    }
    close(S);
    $db->func("\\.\n", 'putline');
    $db->func("endcopy");
}


sub read_block {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($fh,$seek,$ln) = @_;
    my($piece,$readN);

    seek($fh,$seek,0);
    my @lines = ();

        $readN = read($fh,$piece,$ln);
        ($readN == $ln)
            || confess "could not read the block of sims at $seek for $ln characters; $readN actually read";
    return [ split( /\n/, $piece ) ];
}

=head3 bbhs

    my @bbhList = $fig->bbhs($peg, $cutoff);

Return a list of the bi-directional best hits relevant to the specified PEG.

=over 4

=item peg

ID of the feature whose bidirectional best hits are desired.

=item cutoff

Similarity cutoff. If omitted, 1e-10 is used.

=item RETURN

Returns a list of 3-tuples. The first element of the list is the best-hit PEG; the second element is the score. A lower score indicates a better match. The third element is the normalized bit score for the pair, and is normalized to the length of the protein.

=back

=cut
#: Return Type @@;
sub bbhs {
    my($self,$peg,$cutoff) = @_;
    my($sim,$peg2,$genome2,$i,@sims2,%seen);

    if ($self->is_deleted_fid($peg)) { return () }

    $cutoff = defined($cutoff) ? $cutoff : 1.0e-10;

    if ($FIG_Config::use_bbh_server)
    {
	return $self->net_bbhs($peg, $cutoff);
    }

    my @bbhs = ();
    my $rdbH = $self->db_handle;
    my $relational_db_response = $rdbH->SQL("SELECT * FROM bbh WHERE peg1 = \'$peg\' ");

    return sort { $a->[1] <=> $b->[1] }
           grep { $_->[1] <= $cutoff }
           map { [$_->[1],$_->[2],$_->[3]] }
           @{$relational_db_response};
}

=head3 bbh_list

    my $bbhHash = $fig->bbh_list($genomeID, \@featureList);

Return a hash mapping the features in a specified list to their bidirectional best hits
on a specified target genome.

(Modeled after the Sprout call of the same name.)

=over 4

=item genomeID

ID of the genome from which the best hits should be taken.

=item featureList

List of the features whose best hits are desired.

=item RETURN

Returns a reference to a hash that maps the IDs of the incoming features to the best hits
on the target genome.

=back

=cut
#: Return Type %;
sub bbh_list {
    my($self, $genome, $features) = @_;

    my $cutoff = 1.0e-10;

    my $out = {};
    for my $feature (@$features) {
        my @bbhs = $self->bbhs($feature, $cutoff);
        my @featureList = grep { /fig\|$genome\.peg/ } map { $_->[0] } @bbhs;
        $out->{$feature} = \@featureList;
    }
    return $out;
}

=head3 get_figfams_data

usage: $dir = $fig->get_figfams_data($mydir)
usage: $dir = &FIG::get_figfams_data($mydir)

Returns the Figfams data directory to use. If $mydir is passed, use that value. Otherwise
see if $FIG_Config::FigfamsData is defined, and use that. Otherwise default to
$FIG_Config::data/FigfamsData.

=cut

sub get_figfams_data
{
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);

    my ($dir) = @_;

    if (defined($dir))
    {
	return $dir;
    }
    if ($FIG_Config::FigfamsData ne '')
    {
	return $FIG_Config::FigfamsData;
    }
    return "$FIG_Config::data/FigfamsData";
}

=head3 dsims

usage: @sims = $fig->dsims($id,$seq,$maxN,$min_nbsc)

Returns a list of similarities for $seq against PEGs from FIGfams such that

    there will be at most $maxN similarities, and

    each similarity will have a normalized bit-score >= $min_nbsc

The "dsims" or "dynamic sims" are not precomputed.  They are computed using a heuristic which
is much faster than blast, but misses some similarities.  Essentially, you have an "index" or
representative sequences, a quick blast is done against it, and if there are any hits these are
used to indicate which sub-databases to blast against.  This implies that the p-scores are fairly
meaningless; use the normalized bit-scores ($sim->nbsc)

=cut

sub dsims {
    my($self,$id,$seq,$maxN,$min_nbsc,$figfams_data,$blast_parms) = @_;
    my($sim,$partition,%hits);

    $figfams_data = $self->get_figfams_data($figfams_data);

    my $reps_db = "$figfams_data/repdb";

    (-s $reps_db) || return ();

    my @index1 = &blastitP('query',$seq,$reps_db,1.0e-5);
    my %indexH;
    foreach $_ (@index1)
    {
	my $id2 = $_->id2;
	if ($id2 =~ /^(FIG\d+)/)
	{
	    my $nbsc = $_->nbsc;
	    if ((! $indexH{$id2}) || ($indexH{$id2} < $nbsc))
	    {
		$indexH{$id2} = $nbsc;
	    }
	}
    }
    my @index = sort {$indexH{$b} <=> $indexH{$a} } keys(%indexH);
    if (@index > 20) { $#index = 19 }

#   print STDERR "index contains ",scalar @index, " entries\n";
    foreach my $id2 (@index)
    {
	if ($id2 =~ /^(FIG\d+)/)
	{
	    my $fam_id = $1;
#	    print STDERR "\ttrying $fam_id\n";
	    my $fam_dir = &FF::fam_dir($figfams_data,$fam_id);
	    if ((-s "$fam_dir/blast.partition") && open(PARTITION,"<$fam_dir/blast.partition"))
	    {
		if (defined($_ = <PARTITION>) && ($_ =~ /^(\d+)/))
		{
		    my $partition = $1;
		    my $partitionF = "$figfams_data/Partitions/" . ($partition % 1000) . "/$partition/fasta";
		    foreach $sim (&blastitP('query',$seq,$partitionF,1.0e-3))
		    {
			if ($sim->nbsc >= $min_nbsc)
			{
			    my $hit = $sim->id2;
			    if ((! $hits{$hit}) || ($sim->nbsc > $hits{$hit}->nbsc))
			    {
				$sim->[0] = $id;
				$hits{$hit} = $sim;
			    }
			}
		    }
		}
		close(PARTITION);
	    }
	}
    }

    return sort { $b->nbsc <=> $a->nbsc } map { $hits{$_} } keys(%hits);
}

sub blastit {
    my($id,$seq,$db,$maxP,$parms) = @_;
    return &blastitP($id,$seq,$db,$maxP,$parms);
}

sub blastitP {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($id,$seq,$db,$maxP,$parms) = @_;
    $parms ||= '';  # avoid append to undef
    return &blast($id,$seq,$db,$maxP,$parms . " -p blastp");
}

sub blast {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($id,$seq,$db,$maxP,$parms) = @_;

    # Move sanity checks before altering eol

    my $ln1 = length($seq);
    (-s $db) && $ln1 or return ();

    $id = 'unidentified' unless defined($id) && $id =~ /\S/;
    $id = s/^\s+//;
    $id = s/\s.*$//;
    $maxP  ||= 1.0e-5;
    $parms ||= '';      # Avoid concat with undef

    my $old_eol = $/;
    $/ = "\n";

    my $tmpF = "$FIG_Config::temp/blastin.$$.fasta";

    open(TMP,">$tmpF") || die "could not open $tmpF";
    print TMP ">$id\n$seq\n";
    close(TMP);

    my $tool = ($parms =~ /-p (t?blast\s+)/) ? $1 : "blastp";
#   print STDERR "blastall -i $tmpF -d $db -m 8 $parms -e $maxP\n";
    my @blastout = map { chomp; [split(/\s+/,$_),$ln1,undef,$tool] }
                   `$FIG_Config::ext_bin/blastall -i $tmpF -d $db -m 8 $parms -e $maxP`;
    unlink $tmpF;

    my %to_find = map { $_->[1] => 1 } @blastout;
    open(DB,"<$db") || die "could not open $db";

    my %ln2;
    $/ = "\n>";
    while (defined($_ = <DB>))
    {
	chomp;
	if (($_ =~ /^>?(\S+)[^\n]*\n(.*)/sg) && $to_find{$1})
	{
	    my $id2  = $1;
	    my $seq2 = $2;
	    $seq2 =~ s/\s+//g;    # Slight efficiency gain of one-or-more
	    $ln2{$id2} = length($seq2);
	}
    }
    close(DB);

    $/ = $old_eol;
    return map { $_->[13] = $ln2{$_->[1]}; bless($_,'Sim') } @blastout;
}

sub related_by_func_sim {
    my($self,$peg,$user) = @_;
    my($func,$sim,$id2,%related);

    if ($self->is_deleted_fid($peg)) { return () }

    if (($func = $self->function_of($peg,$user)) && (! &FIG::hypo($func)))
    {
        foreach $sim ($self->sims($peg,500,1,"fig",500))
        {
            $id2 = $sim->id2;
            if ($func eq $self->function_of($id2,$user))
            {
                $related{$id2} = 1;
            }
        }
    }
    return keys(%related);
}

################################# chromosomal clusters ####################################

=head3 in_cluster_with

usage: @pegs = $fig->in_cluster_with($peg)

Returns the set of pegs that are thought to be clustered with $peg (on the
chromosome).

=cut

sub in_cluster_with {
    my($self,$peg) = @_;
    my($set,$id,%in);

    return $self->in_set_with($peg,"chromosomal_clusters","cluster_id");
}

=head3 add_chromosomal_clusters

usage: $fig->add_chromosomal_clusters($file)

The given file is supposed to contain one predicted chromosomal cluster per line (either
comma or tab separated pegs).  These will be added (to the extent they are new) to those
already in $FIG_Config::global/chromosomal_clusters.

=cut


sub add_chromosomal_clusters {
    my($self,$file) = @_;
    my($set,$added);

    open(TMPCLUST,"<$file")
        || die "aborted";
    while (defined($set = <TMPCLUST>))
    {
        print STDERR ".";
        chomp $set;
        $added += $self->add_chromosomal_cluster([split(/[\t,]+/,$set)]);
    }
    close(TMPCLUST);

    if ($added)
    {
        my $rdbH = $self->db_handle;
        $self->export_set("chromosomal_clusters","cluster_id","$FIG_Config::global/chromosomal_clusters");
        return 1;
    }
    return 0;
}

#=pod
#
#=head3 export_chromosomal_clusters
#
#usage: $fig->export_chromosomal_clusters
#
#Invoking this routine writes the set of chromosomal clusters as known in the
#relational DB back to $FIG_Config::global/chromosomal_clusters.
#
#=cut
#
sub export_chromosomal_clusters {
    my($self) = @_;

    $self->export_set("chromosomal_clusters","cluster_id","$FIG_Config::global/chromosomal_clusters");
}

sub add_chromosomal_cluster {
    my($self,$ids) = @_;
    my($id,$set,%existing,%in,$new,$existing,$new_id);

#   print STDERR "adding cluster ",join(",",@$ids),"\n";
    foreach $id (@$ids)
    {
        foreach $set ($self->in_sets($id,"chromosomal_clusters","cluster_id"))
        {
            $existing{$set} = 1;
            foreach $id ($self->ids_in_set($set,"chromosomal_clusters","cluster_id"))
            {
                $in{$id} = 1;
            }
        }
    }
#   print &Dumper(\%existing,\%in);

    $new = 0;
    foreach $id (@$ids)
    {
        if (! $in{$id})
        {
            $in{$id} = 1;
            $new++;
        }
    }
#   print STDERR "$new new ids\n";
    if ($new)
    {
        foreach $existing (keys(%existing))
        {
            $self->delete_set($existing,"chromosomal_clusters","cluster_id");
        }
        $new_id = $self->next_set("chromosomal_clusters","cluster_id");
#       print STDERR "adding new cluster $new_id\n";
        $self->insert_set($new_id,[keys(%in)],"chromosomal_clusters","cluster_id");
        return 1;
    }
    return 0;
}

################################# PCH pins  ####################################

=head3 in_pch_pin_with

usage: $fig->in_pch_pin_with($peg)

Returns the set of pegs that are believed to be "pinned" to $peg (in the
sense that PCHs occur containing these pegs over significant phylogenetic
distances).

=cut

sub in_pch_pin_with_old {
    my($self,$peg) = @_;
    my($set,$id,%in);

    return $self->in_set_with($peg,"pch_pins","pin");
}

sub in_pch_pin_with
{
    my($self, $peg1, $diverse) = @_;

    my @all = $self->in_pch_pin_with_and_evidence($peg1);
#    warn "Got all=" . Dumper(\@all);
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

    if ($FIG_Config::use_pch_server)
    {
	warn "returning net pch\n";
	return $self->net_in_pch_pin_with_and_evidence($peg1);
    }

    my $rdbH = $self->db_handle;
    if ($rdbH->table_exists('pchs') && $self->is_complete(&FIG::genome_of($peg1)))
    {
	my $res = $rdbH->SQL(qq(SELECT peg3, max(rep)
				FROM pchs
				WHERE peg1 = ?
				GROUP BY peg3
			       ), undef, $peg1);
	return @$res;
    }
    else
    {
	return ();
    }
}



=head3 add_pch_pins

usage: $fig->add_pch_pins($file)

The given file is supposed to contain one set of pinned pegs per line (either
comma or tab seprated pegs).  These will be added (to the extent they are new) to those
already in $FIG_Config::global/pch_pins.

=cut

sub add_pch_pins {
    my($self,$file) = @_;
    my($set,$added);

    open(TMPCLUST,"<$file")
        || die "aborted";
    while (defined($set = <TMPCLUST>))
    {
        print STDERR ".";
        chomp $set;
        my @tmp = split(/[\t,]+/,$set);
        if (@tmp < 200)
        {
            $added += $self->add_pch_pin([@tmp]);
        }
    }
    close(TMPCLUST);

    if ($added)
    {
        my $rdbH = $self->db_handle;
        $self->export_set("pch_pins","pin","$FIG_Config::global/pch_pins");
        return 1;
    }
    return 0;
}

sub export_pch_pins {
    my($self) = @_;

    $self->export_set("pch_pins","pin","$FIG_Config::global/pch_pins");
}

sub add_pch_pin {
    my($self,$ids) = @_;
    my($id,$set,%existing,%in,$new,$existing,$new_id);

#   print STDERR "adding cluster ",join(",",@$ids),"\n";
    foreach $id (@$ids)
    {
        foreach $set ($self->in_sets($id,"pch_pins","pin"))
        {
            $existing{$set} = 1;
            foreach $id ($self->ids_in_set($set,"pch_pins","pin"))
            {
                $in{$id} = 1;
            }
        }
    }
#   print &Dumper(\%existing,\%in);

    $new = 0;
    foreach $id (@$ids)
    {
        if (! $in{$id})
        {
            $in{$id} = 1;
            $new++;
        }
    }

    if ($new)
    {
        if (keys(%in) < 300)
        {
            foreach $existing (keys(%existing))
            {
                $self->delete_set($existing,"pch_pins","pin");
            }
            $new_id = $self->next_set("pch_pins","pin");
#           print STDERR "adding new pin $new_id\n";
            $self->insert_set($new_id,[keys(%in)],"pch_pins","pin");
        }
        else
        {
            $new_id = $self->next_set("pch_pins","pin");
#           print STDERR "adding new pin $new_id\n";
            $self->insert_set($new_id,$ids,"pch_pins","pin");
        }
        return 1;
    }
    return 0;
}

################################# Annotations  ####################################

=head3 add_annotation

    my $okFlag = $fig->add_annotation($fid, $user, $annotation, $time_made);

Add an annotation to a feature.

=over 4

=item fid

ID of the feature to be annotated.

=item user

Name of the user making the annotation.

=item annotation

Text of the annotation.

=item time_made (optional)

Time of the annotation, in seconds since the epoch. If omitted, the
current time is used.

=item RETURN

Returns 1 if successful, 0 if any of the parameters are invalid or an
error occurs.

=back

=cut

sub add_annotation {
    my($self,$feature_id,$user,$annotation, $time_made) = @_;
    my($genome);

    $time_made = time unless defined( $time_made ) && $time_made =~ /^\d+$/;

    if ($self->is_deleted_fid($feature_id)) { return 0 }

#   print STDERR "add: fid=$feature_id user=$user annotation=$annotation\n";
    if ($genome = $self->genome_of($feature_id))
    {
        my $file = "$FIG_Config::organisms/$genome/annotations";
        my $fileno = $self->file2N($file);
        my $ma   = ($annotation =~ /^Set master function to/);

	$self->broker_log("add_annotation", {
	    user => $user,
	    genome => $genome,
	    fid => $feature_id,
	    annotation => $annotation,
	    time_made => $time_made,
	});

	# Check for a Sapling cross-update.
	my $loader = $self->SaplingCheck();
	if ($loader) {
	    $loader->MakeAnnotation($feature_id, $annotation, $user, $time_made);
	}

        if (open(TMP,">>$file"))
        {
            flock(TMP,LOCK_EX) || confess "cannot lock annotations";
            seek(TMP,0,2)      || confess "failed to seek to the end of the file";

            # Tweaked this section for Windows compatability. The size on disk of
            # "\n" is not constant.
            my $seek1 = tell TMP;
            my $dataLine = "$feature_id\n$time_made\n$user\n$annotation" . ((substr($annotation,-1) eq "\n") ? "" : "\n");
            print TMP $dataLine . "//\n";
            close(TMP);
            chmod 0777, $file;
            my $ln = length($dataLine);
            my $rdbH = $self->db_handle;
            if ($rdbH->SQL("INSERT INTO annotation_seeks ( fid, dateof, who, ma, fileno, seek, len ) VALUES ( \'$feature_id\', $time_made, \'$user\', \'$ma\', $fileno, $seek1, $ln )"))
            {
                return 1;
            }
        }
    }
    return 0;
}

=head3 add_annotation_batch

    my ($n_added, $badList) = $fig->add_annotation_batch($file);

Install a batch of annotations.

=over 4

=item file

File containing annotations.

=item RETURN

Returns the number of annotations successfully added in $n_added. If annotations failed,
they are returned in $badList as a tuple [$peg, $error_msg, $entry].

=back

=cut

#
# This method exists because it is hugely slow to add a large number
# of annotations with add_annotation (it opens and closes the annotation
# file for each individual annotation, and uses individual INSERT statements
# to update the database). This method batches updates to the files and creates
# a load file for the database update.
#
# if the annotations are sorted by genome, so much the better: it will
# do a single file open for the annotation file for that genome.
#

sub add_annotation_batch
{
    my($self, $file) = @_;

    my $anno_fh = new FileHandle("<$file");

    if (not $anno_fh)
    {
        confess "Cannot open $file for reading: $!\n";
    }

    my $dbtmp = "$FIG_Config::temp/add_anno_db.$$";

    my $dbfh = new FileHandle(">$dbtmp");
    if (not $dbfh)
    {
        confess "Cannot write database tmpfile $dbtmp for writing: $!\n";
    }

    local $/ = "///\n";
    my $count = 0;

    my $last_file;
    my $anno_out_fh;
    my $errors = [];

    while (my $anno = <$anno_fh>)
    {
        chomp $anno;

        my ($feature_id, $time_made, $user, $annotation) = split(/\n/, $anno, 4);

        if ($feature_id eq '' or $time_made eq '' or $user eq '' or $annotation eq '')
        {
            push(@$errors, [$feature_id, "Empty fields in annotation", $anno]);
            next;
        }

        next if $self->is_deleted_fid($feature_id);

        my $genome = $self->genome_of($feature_id);
        if (not $genome)
        {
            push(@$errors, [$feature_id, "no genome found for fid '$feature_id'", $anno]);
            next;
        }

        my $file = "$FIG_Config::organisms/$genome/annotations";
        my $fileno = $self->file2N($file);
        my $ma   = ($annotation =~ /^Set master function to/) ? 1 : 0;

        #
        # if this is the first time through or if we have a new file, close and reopen.
        #
        if (not $last_file or $file ne $last_file)
        {
            close($anno_out_fh) if $anno_out_fh;
            chmod 0777, $last_file;
            print "Close $last_file, open $file\n";
            $anno_out_fh = new FileHandle(">>$file");
            if (not $anno_out_fh)
            {
                push(@$errors, [$feature_id, "cannot open annotation file $file: $!", $anno]);
                next;
            }
            $last_file = $file;
            flock($anno_out_fh, LOCK_EX)  or confess "cannot lock assigned_functions $file: $!";
            seek($anno_out_fh, 0, 2)      or confess "failed to seek to the end of the file $file: $!";
        }

        # Tweaked this section for Windows compatability. The size on disk of
        # "\n" is not constant.
        my $seek1 = tell $anno_out_fh;

        my $dataLine = "$feature_id\n$time_made\n$user\n$annotation" . ((substr($annotation,-1) eq "\n") ? "" : "\n");
        print $anno_out_fh $dataLine . "//\n";
        my $ln = length($dataLine);

        print $dbfh join("\t", $feature_id, $time_made, $user, $ma, $fileno, $seek1, $ln), "\n";
        $count++;
    }
    close($anno_out_fh);
    chmod 0777, $last_file;
    print "Loading $count annotations into database from $dbtmp\n";

    close($dbfh);

    my $rows = $self->db_handle()->load_table(file => $dbtmp,
                                              tbl => 'annotation_seeks');
    print "Loaded $rows rows\n";
    return $count, $errors;
}

=head3 merged_related_annotations

usage: @annotations = $fig->merged_related_annotations($fids)

The set of annotations of a set of PEGs ($fids) is returned as a list of 4-tuples.
Each entry in the list is of the form [$fid,$timestamp,$user,$annotation].

=cut

sub merged_related_annotations {
    my($self,$fids) = @_;
    my($fid);
    my(@ann) = ();

    foreach $fid (@$fids)
    {
        push(@ann,$self->feature_annotations1($fid));
    }
    return map { $_->[1] = localtime($_->[1]); $_ } sort { $a->[1] <=> $b->[1] } @ann;
}

=head3 feature_annotations

    my @annotations = $fig->feature_annotations($fid, $rawtime);

Return a list of the specified feature's annotations. Each entry in the list
returned is a 4-tuple containing the feature ID, time stamp, user ID, and
annotation text. These are exactly the values needed to add the annotation
using L</add_annotation>, though in a different order.

=over 4

=item fid

ID of the features whose annotations are to be listed.

=item rawtime (optional)

If TRUE, the times will be returned as PERL times (seconds since the epoch);
otherwise, they will be returned as formatted time strings.

=item RETURN

Returns a list of 4-tuples, one per annotation. Each tuple is of the form
I<($fid, $timeStamp, $user, $annotation)> where I<$fid> is the feature ID,
I<$timeStamp> is the time the annotation was made, I<$user> is the name of
the user who made the annotation, and I<$annotation> is the text of the
annotation.

=back

=cut

sub feature_annotations {
    my($self,$feature_id,$rawtime) = @_;

    if ($self->is_deleted_fid($feature_id)) { return () }

    if ($rawtime)
    {
        return $self->feature_annotations1($feature_id);
    }
    else
    {
        return map {  $_->[1] = localtime($_->[1]); $_ } $self->feature_annotations1($feature_id);
    }
}

sub feature_annotations1 {
    my($self,$feature_id) = @_;
    my($tuple,$fileN,$seek,$ln,$annotation,$feature_idQ);
    my($file,$fh);

    if ($self->is_deleted_fid($feature_id)) { return () }

    my $rdbH = $self->db_handle;
    my $relational_db_response = $rdbH->SQL("SELECT fileno, seek, len  FROM annotation_seeks WHERE  fid = \'$feature_id\' ");
    my @annotations = ();

    foreach $tuple (@$relational_db_response)
    {
        ($fileN,$seek,$ln) = @$tuple;
        $annotation = $self->read_annotation($fileN,$seek,$ln);
        $feature_idQ = quotemeta $feature_id;
        if ($annotation =~ /^$feature_idQ\n(\d+)\n([^\n]+)\n(.*)/s)
        {
            push(@annotations,[$feature_id,$1,$2,$3]);
        }
        else
        {
            print STDERR "malformed annotation for $feature_idQ\n$annotation\n";
        }
    }
    return sort { $a->[1] <=> $b->[1] } @annotations;
}

sub read_annotation {
    my($self,$fileN,$seek,$ln) = @_;
    my($readN,$readC);

    my $file = $self->N2file($fileN);
    if (! $file) { return "" }

    my $fh   = $self->openF($file);
    if (! $fh)
    {
        confess "could not open annotations for $file";
    }

    #
    # See if the seek address is after the end of the file. If it is,
    # we're likely looking at an annotation file that is older than the
    # database entry. This can come from instantaneous database replication
    # with file replication happening at a slower rate.
    #

    if ($seek > -s $fh)
    {
        warn "Attempting to seek past the end of $file\n";
        return "";
    }

    seek($fh,$seek,0);
    $readN = read($fh,$readC,$ln);
    my $len2 = length $readC;
    ($readN == $ln)
        || confess "could not read the block of annotations at $seek for $ln characters; $readN actually read from $file\n$readC";
    return $readC;
}

sub epoch_to_readable {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($epoch) = @_;

    my($sec,$min,$hr,$dd,$mm,$yr) = localtime($epoch);
    $mm++;
    $yr += 1900;
    return "$mm-$dd-$yr:$hr:$min:$sec";
}

=head3 read_all_annotations

    my @annotations = $fig->read_all_annotations($genomeID);

Return a list of the specified genome's annotations. Each entry in the list
returned is a 4-tuple containing the feature ID, time stamp, user ID, and
annotation text. The values are read directly from the annotation flat
file without resorting to the database.

=over 4

=item genomeID

ID of the genome whose annotations are to be read.

=item RETURN

Returns a list of 4-tuples, one per annotation. Each tuple is of the form
I<($fid, $timeStamp, $user, $annotation)> where I<$fid> is the feature ID,
I<$timeStamp> is the time the annotation was made, I<$user> is the name of
the user who made the annotation, and I<$annotation> is the text of the
annotation.

=back

=cut
#: Return Type ;
sub read_all_annotations {
    # Get the parameters.
    my ($self, $genomeID) = @_;
    # Declare the return variable.
    my @retVal = ();
    # Locate the genome's annotation file.
    my $annoFileName = "$FIG_Config::organisms/$genomeID/annotations";
    # Only proceed if it exists. If it doesn't, we have no annotations.
    if (-e $annoFileName) {
        # Open the file.
        Open(\*ANNOTATIONS, "<$annoFileName");
        # Loop through the file.
        while (my $record = read_annotation_record(\*ANNOTATIONS)) {
            # Clear the trailing newline.
            chomp $record;
            # Split out the parts.
            my ($featureID, $time, $user, @data) = split /\s*\n/, $record;
            # Rejoin the data records.
            my $data = join("\n", @data);
            # Verify the feature ID.
            if (! $self->is_deleted_fid($featureID)) {
                push @retVal, [$featureID, $time, $user, $data];
            }
        }
    }
    # Return the result.
    return @retVal;
}

=head3 read_annotation_record

    my $annoString = FIG::read_annotation_record($fileHandle);

Read an annotation record from the specified file handle. Will return the
annotation record if successful, and C<undef> if end-of-file is read. An
annotation record consists of multiple lines of text separated by a
line containing a double-slash C<//>.

=over 4

=item fileHandle

The file handle from which to read the record.

=item RETURN

Returns either the entire annotation record (without the double-slash) or
C<undef>, indicating end-of-file. Null records will not be returned.

=back

=cut
#: Return Type ;
sub read_annotation_record {
    # Get the parameters.
    my ($fileHandle) = @_;
    # Declare the return variable.
    my $retVal = "";
    # Loop until we find a non-null record or end-of-file.
    while (defined($retVal) && $retVal eq "") {
        # Loop through the file records, stuffing them into the return
        # variable.
        my $line = <$fileHandle>;
        while (defined($line) && $line !~ m!^//!) {
            $retVal .= $line;
            $line = <$fileHandle>;
        }
        # Check for the end-of-file possibility.
        if (!defined($line)) {
            $retVal = undef;
        }
    }
    # Return the result.
    return $retVal;
}

=head3 parse_date

usage: $date = $fig->parse_date(date-string)

Parse a date string, returning seconds-since-the-epoch, or undef if the date did not parse.

Accepted formats include an integer, which is assumed to be seconds-since-the-epoch an
is just returned; MM/DD/YYYY;  or a date that can be parsed by the routines in
the Date::Parse module.

=cut

sub parse_date
{
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);

    my($date) = @_;

    $date or return undef;

    my $epoch_date;

    if ($date =~ /^(\d{1,2})\/(\d{1,2})\/(\d{4})$/)
    {
        my($mm,$dd,$yyyy) = ($1,$2,$3);
        $epoch_date = &Time::Local::timelocal(0,0,0,$dd,$mm-1,$yyyy-1900,0,0,0);
    }
    elsif ($date =~ /^\d+$/)
    {
        $epoch_date = $date;
    }
    elsif ($haveDateParse)
    {
        $epoch_date = str2time($date);
    }
    return $epoch_date;
}

#
# This now calls assignments_made_full and remaps the output.
#
sub assignments_made
{
    my($self,$genomes,$who,$date) = @_;

    my @a = $self->assignments_made_full($genomes, $who, $date);

    return map { [ @{$_}[0,1]] } @a;
}

#
# Looks up and returns assignments made; return is a list of
# tuples [peg, assignment, date, who]
#

sub assignments_made_full {
    my($self,$genomes,$who,$date) = @_;
    my($relational_db_response,$entry,$fid,$fileno,$seek,$len,$ann);
    my($epoch_date,$when,%sofar,$x);

    if (! defined($genomes)) { $genomes = [$self->genomes] }

    my %genomes = map { $_ => 1 } @$genomes;

    $epoch_date = $self->parse_date($date);

    $epoch_date = defined($epoch_date) ? $epoch_date-1 : 0;

    my @assignments = ();
    my $rdbH = $self->db_handle;
    if ($who eq "master")
    {
        $relational_db_response = $rdbH->SQL("SELECT fid, dateof, fileno, seek, len  FROM annotation_seeks WHERE ((ma = \'1\') AND (dateof > $epoch_date))");
    }
    else
    {
        $relational_db_response = $rdbH->SQL("SELECT fid, dateof, fileno, seek, len  FROM annotation_seeks WHERE (( who = \'$who\' ) AND (dateof > $epoch_date))");
    }

    if ($relational_db_response && (@$relational_db_response > 0))
    {
        foreach $entry (@$relational_db_response)
        {
            ($fid,$when,$fileno,$seek,$len) = @$entry;
            if (($fid =~ /^fig\|(\d+\.\d+)/) && $genomes{$1} && (! $self->is_deleted_fid($fid)))
            {
                if ($len < 4)
                {
                    print STDERR "BAD: fid=$fid when=$when fileno=$fileno seek=$seek len=$len\n";
                    next;
                }
                $ann = $self->read_annotation($fileno,$seek,$len);

                if (($ann =~ /^(fig\|\d+\.\d+\.peg\.\d+)\n(\d+)\n(\S+)\nSet ([^\n]*)function[^\n]*\n(\S[^\n]+\S)/s) &&
                    (($who eq $3) || (($4 eq "master ") && ($who eq "master"))) &&
                    ($2 >= $epoch_date))
                {
                    if ((! $sofar{$1}) || (($x = $sofar{$1}) && ($when > $x->[0])))
                    {
                        $sofar{$1} = [$when, $5, $3];
                    }
                }
            }
        }
    }
    @assignments = map { $x = $sofar{$_}; [$_,$x->[1], $x->[0], $x->[2]] } keys(%sofar);
    return @assignments;
}

=head3 extract_assignments_from_annotations

Extract a list of assignments from an annotations package as created by
annotations_made_fast. Assumes that the user and date filtering was
done by the annotations extraction, so all this has to do is to
sort the lists of annotations by date and grab the latest one.

Return value is a list of tuples [$peg, $assignment, $date, $who].

=cut

sub extract_assignments_from_annotations
{
    my($self, $annos) = @_;

    #
    # $annos is a list of pairs [$genome, $genomeannos]
    # $genomeannos is a hash keyed on peg. value is a list of lists [$peg, $time, $who, $anno].
    #

    #
    # Sort on genome.
    #
    my @annos = sort { &FIG::by_genome_id($a->[0], $b->[0]) } @$annos;

    my @out;
    for my $gent (@annos)
    {
        my($genome, $genome_anno_list) = @$gent;

        #
        # Sort on peg id.
        for my $peg (sort { &FIG::by_fig_id($a, $b) } keys %$genome_anno_list)
        {
            my $anno_list = $genome_anno_list->{$peg};

            #
            # Pull assignment annotations.
            #

            my @a = grep { $_->is_assignment() } @$anno_list;

            next unless @a > 0;

            #
            # and sort by date, descending.
            #

            @a = sort { $b->anno_time() <=> $a->anno_time() } @a;

            my $winner = $a[0];

            $winner->fid() eq $peg or confess "KEY mismatch in annotations_made_fast output";

            push(@out, $winner);
        }
    }
    return @out;
}

sub assignments_made_for_protein {
    my($self, $fid) = @_;
    my($relational_db_response,$entry,$fileno,$seek,$len,$ann);
    my($epoch_date,$when,%sofar,$x);

    if ($self->is_deleted_fid($fid)) { return () }

    my @assignments = ();
    my $rdbH = $self->db_handle;

    $relational_db_response = $rdbH->SQL("SELECT fid, dateof, fileno, seek, len  FROM annotation_seeks WHERE (fid = '$fid')");

    if ($relational_db_response && (@$relational_db_response > 0))
    {
        foreach $entry (@$relational_db_response)
        {
            ($fid,$when,$fileno,$seek,$len) = @$entry;
            if ($len < 4)
            {
                print STDERR "BAD: fid=$fid when=$when fileno=$fileno seek=$seek len=$len\n";
                next;
            }
            $ann = $self->read_annotation($fileno,$seek,$len);

            if (my ($peg, $when, $who, $what, $func) =
                $ann =~ /^(fig\|\d+\.\d+\.peg\.\d+)\n(\d+)\n(\S+)\nSet ([^\n]*)function[^\n]*\n(\S[^\n]+\S)/s)
            {
                push(@assignments, [$peg, $when, $who, $what, $func]);
            }
        }
    }
    return @assignments;
}

=head3 annotations_made

usage: @annotations = $fig->annotations_made($genomes, $who, $date)

Return the list of annotations on the genomes in @$genomes  made by $who
after $date.

Each returned annotation is of the form [$fid,$timestamp,$user,$annotation].

=cut

sub annotations_made {
    my($self,$genomes,$who,$date) = @_;
    my($relational_db_response,$entry,$fid,$fileno,$seek,$len,$ann);
    my($epoch_date,$when,@annotations);

    if (! defined($genomes)) { $genomes = [$self->genomes] }

    my %genomes = map { $_ => 1 } @$genomes;
    if ($date =~ /^(\d{1,2})\/(\d{1,2})\/(\d{4})$/)
    {
        my($mm,$dd,$yyyy) = ($1,$2,$3);
        $epoch_date = &Time::Local::timelocal(0,0,0,$dd,$mm-1,$yyyy-1900,0,0,0);
    }
    elsif ($date =~ /^\d+$/)
    {
        $epoch_date = $date;
    }
    else
    {
        $epoch_date = 0;
    }
    $epoch_date = defined($epoch_date) ? $epoch_date-1 : 0;
    @annotations = ();
    my $rdbH = $self->db_handle;
    if ($who eq "master")
    {
        $relational_db_response = $rdbH->SQL("SELECT fid, dateof, fileno, seek, len  FROM annotation_seeks WHERE ((ma = \'1\') AND (dateof > $epoch_date))");
    }
    else
    {
        $relational_db_response = $rdbH->SQL("SELECT fid, dateof, fileno, seek, len  FROM annotation_seeks WHERE (( who = \'$who\' ) AND (dateof > $epoch_date))");
    }

    if ($relational_db_response && (@$relational_db_response > 0))
    {
        foreach $entry (@$relational_db_response)
        {
            ($fid,$when,$fileno,$seek,$len) = @$entry;
            if (($fid =~ /^fig\|(\d+\.\d+)/) && $genomes{$1} && (! $self->is_deleted_fid($fid)))
            {
                $ann = $self->read_annotation($fileno,$seek,$len);

                if ($ann =~ /^(fig\|\d+\.\d+\.peg\.\d+)\n(\d+)\n(\S+)\n(.*\S)/s)
                {
                    push(@annotations,[$1,$2,$3,$4]);
                }
            }
        }
    }
    return @annotations;
}

sub annotations_made_fast
{
    my($self, $genomes, $start_time, $end_time, $anno_by, $replace_master_with_group) = @_;

    if (!defined($anno_by))
    {
        $anno_by = 'master';
    }

    if (!defined($genomes))
    {
        $genomes = [$self->genomes()];
    }

    my $group = $FIG_Config::group;

    $group = 'FIG' unless $group;

    my $annos;
    my $pegs = {};

    if ($start_time !~ /^\d+$/)
    {
        my $st = parse_date($start_time);
        if (!defined($st))
        {
            confess "annotations_made_fast: unparsable start time '$start_time'";
        }
        $start_time = $st;
    }
    if (defined($end_time))
    {
        if ($end_time !~ /^\d+$/)
        {
            my $et = parse_date($end_time);
            if (!defined($et))
            {
                confess "annotations_made_fast: unparsable end time '$end_time'";
            }
            $end_time = $et;
        }
    }
    else
    {
        $end_time = time + 60;
    }

    #
    # We originally used a query to get the PEGs that needed to have annotations
    # sent. Unfortunately, this performed very poorly due to all of the resultant
    # seeking around in the annotations files.
    #
    # The code below just runs through all of the anno files looking for annos.
    #
    # A better way to do this would be to do a query to retrieve the genome id's for
    # genomes that have updates. The problem here is that the annotation_seeks
    # table doesn't have an explicit genome field.
    #
    # Surprisingly, to me anyway, the following query appers to run quickly, in both
    # postgres and mysql:
    #
    # SELECT distinct(substring(fid from 5 for position('.peg.' in fid) - 5))
    # FROM annotation_seeks
    # WHERE dateof > some-date.
    #
    # The output of that can be parsed to get the genome id and just those
    # annotations files searched.
    #

    my $master_anno = $anno_by eq 'master';

    for my $genome (@$genomes)
    {
        my $genome_dir = "$FIG_Config::organisms/$genome";
        next unless -d $genome_dir;
        my $gpegs = {};

        my $afh = new FileHandle("<$genome_dir/annotations");
        if ($afh)
        {
            my($fid, $anno_time, $who, $anno_text,$anno_who, @rest);

            while (not $afh->eof())
            {
                chomp($fid = <$afh>);
                next if $fid eq "//";
                chomp($anno_time = <$afh>);
                next if $anno_time eq "//";
                chomp($who = <$afh>);
                next if $who eq "//";
                @rest = ();

                while (<$afh>)
                {
                    chomp;
                    last if $_ eq "//";
                    push(@rest, $_);
                }

                #
                # Validate.
                #

                if ($fid !~ /^fig\|\d+\.\d+\.peg\.\d+$/)
                {
                    #warn "Invalid fid '$fid' in annotations ($genome_dir/annotations line $.)\n";
                    next;
                }
                elsif ($anno_time !~ /^\d+$/)
                {
                    warn "Invalid annotation time '$anno_time' in annotations ($genome_dir/annotations line $.)\n";
                    next;
                }

                #
                # Filter deleted fids.
                #

                next if $self->is_deleted_fid($fid);

                #
                # Filter on date.
                #

                next if $anno_time < $start_time or $anno_time > $end_time;

                my $aobj = new Annotation($fid, $anno_time, $who, @rest);

                if ($aobj->is_assignment())
                {
                    my $anno_who = $aobj->assignment_who();

                    #
                    # Filter on annotator.
                    #
                    if ($anno_by eq 'all' or
                        ($master_anno ?
                         ($anno_who eq 'FIG' or $anno_who eq 'master') :
                         ($who eq $anno_by)))
                    {
                        if ($replace_master_with_group)
                        {
                            $aobj->set_assignment_who($group);
                        }
                    }
                    else
                    {
                        next;
                    }
                }
                #
                # Non-assignment annotations are filtered such that:
                # If master annotations are requested, we take all non-assignment annotations.
                # Otherwise, we take only those where $who eq $anno_by.
                #
                elsif (not($master_anno or
                           $anno_by eq 'all' or
                           $anno_by eq $who))
                {
                    next;
                }

                #
                # Fall through: save this anno. Note that we do not save the final newline.
                #

                push(@{$gpegs->{$fid}}, $aobj);
            }

#           while (my $ann = <$afh>)
#           {
#               chomp $ann;

#               if ((($fid, $anno_time, $who, $anno_text, $anno_who) =
#                    ($ann =~ /^(fig\|\d+\.\d+\.peg\.\d+)\n(\d+)\n(\S+)\n(Set\s+(\S+)\s+function\s+to.*\S)/s)) and
#                   not $self->is_deleted_fid($fid) and
#                   $anno_time >= $start_time and
#                   $anno_time <= $end_time and
#                   ($anno_by eq 'all' or ($master_anno ? ($anno_who eq 'FIG' or $anno_who eq 'master') : ($who eq $anno_by))))
#               {
#                   #
#                   # Update users list.
#                   #
#                   {
#                       my $d =  $self->is_deleted_fid($fid);
#                   }

#                   if ($replace_master_with_group)
#                   {
#                       $anno_text =~ s/Set master function to/Set $group function to/;
#                   }

#                   my $anno = [$fid, $anno_time, $who, $anno_text];

#                   push(@{$gpegs->{$fid}}, $anno);
#               }
#           }

        }
        push(@$annos, [$genome, $gpegs]);
    }
    return $annos;
}

################################### ATTRIBUTES

=head2 Attributes

The attribute system automatically detects whether you are using a local
attribute database, a remote attribute server, or the SEED data store. For
details on the new attribute system see the documentation for the
B<CustomAttributes> module.

Because of the enormous number of attributes in the system (1.5 million and
growing), the old system, which combined a database table and flat file data
stores, has become too slow for live SEEDs. It is maintained for small test
SEEDs, such as what you might have running on a local PC. Be aware, however,
that not all functions of the old system work in the new system, and vice versa.

To activate the new attribute system, define C<$FIG_Config::attrHost> to be the
name of the server containing the attribute database (currently the annotator
SEED).

=head3 The SEED Data Store Interface

There are several base attribute methods:

 get_attributes
 add_attribute
 delete_attribute
 change_attribute

There are also methods for more complex things:

 get_keys
 get_values
 guess_value_format

By default all keys are case sensitive, and all keys have leading and trailing white space removed. Keys can not contain anything but [a-zA-Z0-9_] (or things matched by \w)

Attributes are not on a 1:1 correlation, so a single key can have several values.

Most attributes files are stored in the genome specific directories. These are in Organisms/nnnnn.n/Attributes for the organisms, and Organisms/nnnnn.n/Feaures/*/Attributes for the features. Attributes can also be stored in Global/Attributes where they will be loaded, but users are discouraged from doing this since there will be no packaging and sharing of those attibutes. Global should be reserved for those attributes that are calculated on a database-wide instance. There are several "special" files that we are using:

1. Definition files

These are the raw text files stored in the appropriate locations (Organisms/nnnnn.n/Attributes, Organisms/nnnnn.n/Feaures/*/Attributes, and Global/Attributes). The files should consist of ONLY feature, key, value, and optional URL. Any other columns will be ignored and not loaded into the database.

2. Global/Attributes/attribute_keys

This contains the definition of the attribute keys. There are currently 3 defined columns although others may be added and this file can contain lines of an arbitrary length.

3. Global/Attributes/transaction_log, Organisms/nnnnnn.n/Attributes/transaction_log, and Organisms/nnnnnn.n/Features/*/Attributes/transaction_log

These are the transaction logs that contain any modifications to the data. In general the data is loaded from a single definition file this is not modified by the software. Any changes to the attributes are made in the Database and then written to the transaction log. The transaction log has the following columns

1. command. This can be one of ADD/DELETE/CHANGE
2. feature. The feature id to be modified
3. key. The key to be modified
4. old value. The original value of the key
5. old url. The original URL
6. new value. The new value for the key. Ignored if the key is deleted.
7. new url. The new value for the URL. Ignored if the key is deleted.

Note that the old value and old url are optional. If they are not provided ALL instances of the key will be affected.

Notice also that the old file assigned_attributes is no longer used. This is replaced by the transaction log.

Finally, in the parsing of all files any line beginning with a pound sign is ignored as a comment.

A method, read_attribute_transaction_log, is provided to read the transaction_logs and implement the changes therein. In each of the methods add_attribute, delete_attribute, and change_attribute there is an optional boolean that can be set to prevent writing of the transaction_log. The read_attribute_transaction_log reads the log and then adds/changes/deletes the records as appropriate. Without this boolean there is a circular reference.

Get attributes requires one of four keys:
fid (which can be genome, peg, rna, or other id, or a reference to a list of ids),
key,
value,
url

It will find any attribute that has the characteristics that you request, and if any values match it will return a four-ple of:
[fid, key, value, url]

You can request an E. coli key like this
$fig->get_attributes('83333.1');

You can request a peg id like this:
$fig->get_attributes($peg);
$fig->get_attributes("fig|833333.1.peg.4");

You can request any structure key like this
$fig->get_attributes(undef, 'structure');

You can request any url like this
$fig->get_attributes(undef, undef, undef, 'http://pir.georgetown.edu/sfcs-cgi/new/pirclassif.pl?id=SF001547');

NOTE: If there are no attributes an empty array will be returned. You need to check for this and not assume that it will be undef.

=head3 get_attributes

    my @attributeList = $fig->get_attributes($objectID, $key, @values);

In the database, attribute values are sectioned into pieces using a splitter
value specified in the constructor (L</new>). This is not a requirement of
the attribute system as a whole, merely a convenience for the purpose of
these methods. If a value has multiple sections, each section
is matched against the corresponding criterion in the I<@valuePatterns> list.

This method returns a series of tuples that match the specified criteria. Each tuple
will contain an object ID, a key, and one or more values. The parameters to this
method therefore correspond structurally to the values expected in each tuple. In
addition, you can ask for a generic search by suffixing a percent sign (C<%>) to any
of the parameters. So, for example,

    my @attributeList = $attrDB->GetAttributes('fig|100226.1.peg.1004', 'structure%', 1, 2);

would return something like

    ['fig|100226.1.peg.1004', 'structure', 1, 2]
    ['fig|100226.1.peg.1004', 'structure1', 1, 2]
    ['fig|100226.1.peg.1004', 'structure2', 1, 2]
    ['fig|100226.1.peg.1004', 'structureA', 1, 2]

Use of C<undef> in any position acts as a wild card (all values). You can also specify
a list reference in the ID column. Thus,

    my @attributeList = $attrDB->GetAttributes(['100226.1', 'fig|100226.1.%'], 'PUBMED');

would get the PUBMED attribute data for Streptomyces coelicolor A3(2) and all its
features.

In addition to values in multiple sections, a single attribute key can have multiple
values, so even

    my @attributeList = $attrDB->GetAttributes($peg, 'virulent');

which has no wildcard in the key or the object ID, may return multiple tuples.

Value matching in this system works very poorly, because of the way multiple values are
stored. For the object ID and key name, we create queries that filter for the desired
results. For the values, we do a comparison after the attributes are retrieved from the
database. As a result, queries in which filter only on value end up reading the entire
attribute table to find the desired results.

=over 4

=item objectID

ID of object whose attributes are desired. If the attributes are desired for multiple
objects, this parameter can be specified as a list reference. If the attributes are
desired for all objects, specify C<undef> or an empty string. Finally, you can specify
attributes for a range of object IDs by putting a percent sign (C<%>) at the end.

=item key

Attribute key name. A value of C<undef> or an empty string will match all
attribute keys. If the values are desired for multiple keys, this parameter can be
specified as a list reference. Finally, you can specify attributes for a range of
keys by putting a percent sign (C<%>) at the end.

=item values

List of the desired attribute values, section by section. If C<undef>
or an empty string is specified, all values in that section will match. A
generic match can be requested by placing a percent sign (C<%>) at the end.
In that case, all values that match up to and not including the percent sign
will match. You may also specify a regular expression enclosed
in slashes. All values that match the regular expression will be returned. For
performance reasons, only values have this extra capability.

=item RETURN

Returns a list of tuples. The first element in the tuple is an object ID, the
second is an attribute key, and the remaining elements are the sections of
the attribute value. All of the tuples will match the criteria set forth in
the parameter list.

=back

=cut

sub get_attributes {
    my($self,@request) = @_;

    # @request = (  $id,  $key, $value, $url )
    # or         ( \@ids, $key, $value, $url )


    # Create the return list.
    my @attr;

    if (exists $self->{_ca}) {
        # Here we can use the new system.
        @attr = $self->{_ca}->GetAttributes(@request);
    } else {
        # Make sure that the id is a reference to a list of ids
        if (! defined($request[0]) || ref($request[0]) ne 'ARRAY' ) {
            $request[0] = [$request[0]];
        }

        ##  Should we still be doing this?  My guess is not.  It excludes
        ##  wildcard searches. Valid key names should be the responsibility
        ##  of the caller.
        # clean the keys if there are any
        # $request[1] && ($request[1] = $self->clean_attribute_key($request[1]));
        ####

        my $rdbH = $self->db_handle;
        return () unless ($rdbH);

        # An error check to make sure that we are operating on the new version of attributes
        # If we are not, we will print an error and then return. Otherwise continue
        eval {$rdbH->SQL("SELECT genome,ftype,id,tag,val,url FROM attribute LIMIT 1")};
        if ( $@ ) { return () }  # This was returning reference, not list
        #if ($@) {print STDERR "Please rerun load_attributes to install the newest set of attributes\n"; return []}

        # Loop through the ids.
        my ( $ids, $tag, $val, $url ) = @request;
        foreach my $fid ( @$ids ) {
            # Table columns are now genome, ftype, id, key, val, url
            # here we generate the select statement based on what is in the request. Only add those fields we need.
            # we add the where conditional to the @where array and the value for that conditional to the @values array
            # and then join the @where into the select statement. The @values is provided to the SQL statement to merge
            my @where;
            my @values;
            # Fix this to handle genome ids, and protein sequence md5 ids
            if ( $fid ) {
                my ( $genome, $ftype, $id ) = $self->split_attribute_oid( $fid );
                if ( $genome ) {
                    push @where,  $genome =~ /%/ ? 'genome LIKE ?' : 'genome = ?';
                    push @values, $genome;
                }
                if ( $ftype ) {
                    push @where,  $ftype =~ /%/ ? 'ftype LIKE ?' : 'ftype = ?';
                    push @values, $ftype;
                }
                if ( $id ) {
                    push @where,  $id =~ /%/ ? 'id LIKE ?' : 'id = ?';
                    push @values, $id;
                }
            }
            if ( defined $tag ) {
                push @where,  $tag =~ /%/ ? 'tag LIKE ?' : 'tag = ?';
                push @values, $tag;
            }
            if ( defined $val ) {
                push @where,  $val =~ /%/ ? 'val LIKE ?' : 'val = ?';
                push @values, $val;
            }
            if ( $url ) {
                push @where,  $url =~ /%/ ? 'url LIKE ?' : 'url = ?';
                push @values, $url;
            }

            my $select = "SELECT genome,ftype,id,tag,val,url FROM attribute";
            $select   .= " WHERE (" . join(' AND ', @where) . ")" if @where;
            #print STDERR "TRYING: $select and ?=", join(" ?=", @values), "\n";

            Trace("Where list for attributes is (" . join(", ", @where) . ")") if T(4);
            Trace("Value list for attributes is (" . join(", ", @values) . ")") if T(4);
            Trace("Select for attributes is: $select") if T(4);
	    #print STDERR "$select\n", join("\n", @values), "\n";
            my $res=$rdbH->SQL($select, undef, @values);

            # the following line takes the first 3 elements from each array and puts them back
            # to be a feature or genome using join_attribute_oid and then puts them back in the array.
	    #print STDERR Dumper($res);
            map { unshift @$_, $self->join_attribute_oid( splice(@$_, 0, 3) ) } @$res;
            push @attr, @{$res};
        }
    }

    # Evidence codes of proteins are keyed by sequence, not by fid. So, if
    # the user is asking for all attributes or evidence codes of a feature,
    # then we need to do a search keyed by the associated sequence.
    my ( $ids, $key ) = @request;
    $ids = [ $ids ] if $ids && ref($ids) ne 'ARRAY';

    if ( $ids && ! defined $key || $key eq 'evidence_code' ) {
        push @attr, map  { Dlits::protein_evidence_codes($self,$_) }
                    grep { $_ && /^fig\|\d+\.\d+\.peg\.\d+$/ }
                    @$ids;
    }
    return @attr;
}

=head3 query_attributes

    my @attributeData = $ca->query_attributes($filter, $filterParms);

Return the attribute data based on an SQL filter clause. In the filter clause,
the name C<$object> should be used for the object ID, C<$key> should be used for
the key name, C<$subkey> for the subkey value, and C<$value> for the value field.

=over 4

=item filter

Filter clause in the standard ERDB format, except that the field names are C<$object> for
the object ID field, C<$key> for the key name field, C<$subkey> for the subkey field,
and C<$value> for the value field. This abstraction enables us to hide the details of
the database construction from the user.

=item filterParms

Parameters for the filter clause.

=item RETURN

Returns a list of tuples. Each tuple consists of an object ID, a key (with optional subkey), and
one or more attribute values.

=back

=cut

sub query_attributes {
    my ($self, $filter, $filterParms) = @_;
    return $self->{_ca}->QueryAttributes($filter, $filterParms);
}


=head3 get_cv_attributes

A simple wrapper around get_attriubtes to return only those attributes
that have meta_data indicating that the key is a controlled vocabulary.

### DEPRECATED ### The controlled vocabulary feature was never used in the old
system, and in the new system, ALL the keys are controlled vocabulary.

=cut

sub get_cv_attributes {
    return get_attributes(@_);
}

=head3 add_attribute

Add a new key/value pair to something. Something can be a genome id, a peg, an rna, prophage, whatever.

Arguments:

        feature id, this can be a peg, genome, etc,
        key name. This is case sensitive and has the leading and trailing white space removed
        value
        optional URL to add
        boolean to prevent writing to the transaction log. See above

=cut

sub add_attribute {
    my ($self, @request) = @_;
    if (exists $self->{_ca}) {
        # Here we can use the new system.
        return $self->{_ca}->AddAttribute(@request);
    } else {
        my($peg,$k,$v, $url, $notl) = @request;
        return unless ($peg && $k); # we must have at least a peg and a tag to add (though other things can be undef)
        $k =  $self->clean_attribute_key($k);
        my $rdbH = $self->db_handle;

        # split the peg/feature/genome into pieces and parts
        $rdbH->SQL("INSERT INTO attribute ( genome,ftype,id,tag,val,url ) VALUES ( ?,?,?,?,?,? )", undef, $self->split_attribute_oid($peg),$k,$v,$url);
        my $location=$self->attribute_location($peg);
        &verify_dir("$location");
        if (!$notl && open(TMPATTR,">>$location/transaction_log"))
        {
            print TMPATTR "ADD\t$peg\t$k\t$v\t$url\n";
            close(TMPATTR);
        }
        return 1;
    }
}

=head3 delete_attribute

    $fig->delete_attribute($objectID, $key, @values);

Delete the specified attribute key/value combination from the database.

=over 4

=item objectID

ID of the object whose attribute is to be deleted.

=item key

Attribute key name.

=item values

One or more values associated with the key. If no values are specified, then all values
will be deleted. Otherwise, only a matching value will be deleted.

=back

=cut

sub delete_attribute {
    my ($self, @request) = @_;
    if (exists $self->{_ca}) {
        # Here we can use the new system.
        return $self->{_ca}->DeleteAttribute(@request);
    } else {
        my($peg,$k, $oldval, $oldurl, $notl) = @request;

        # we need a feature and a key to delete
        return unless ($peg && $k);

        # clean the key
        $k =  $self->clean_attribute_key($k);
        # get the transaction log
        my $location=$self->attribute_location($peg);
        &verify_dir("$location");
        if (!$notl && open(TMPATTR,">>$location/transaction_log"))
        {
            print TMPATTR "DELETE\t$peg\t$k\t$oldval\n";
            close(TMPATTR);
        }
        return $self->change_attribute($peg, $k, $oldval, $oldurl, undef, undef);
    }
}

=head3 parse_oid

    my ($type, $id) = FIG::parse_oid($idValue);

Convert an attribute object ID to an object type and an ID applicable to that type.
This information can be used to convert an ID string obtained from the L</get_attributes>
method to an object name and ID suitable for plugging into the C<GetEntity> method
of an B<ERDB> database.

=over 4

=item idValue

ID string from the attribute database.

=item RETURN

Returns a two-element list consisting of the object type and its individual ID.

=back

=cut

sub parse_oid {
    my ($idValue) = @_;
    my @retVal = CustomAttributes::ParseID($idValue);
    return @retVal;
}

=head3 form_oid

    my $idValue = FIG::form_oid($type, $id);

Convert an object type and ID into an ID string for the attribute database.

=over 4

=item type

Object type. This should usually correspond to an entity name in a database. It can
only contain letters. This means no digits, spaces, or even underscores.

=item id

Individual object ID.

=item RETURN

Returns the string used to represent the object in the attribute database.

=back

=cut

sub form_oid {
    # Get the parameters.
    my ($type, $id) = @_;
    my $retVal = CustomAttributes::FormID($type, $id);
    return $retVal;
}

=head3 delete_matching_attributes

    my @attributeList = $fig->delete_matching_attributes($objectID, $key, @values);

This method works identically to L</get_attributes>, except that the attributes are
deleted as they are retrieved.

=cut

sub delete_matching_attributes {
    # Get the parameters.
    my ($self, $objectID, $key, @values) = @_;
    my @retVal;
    # Declare the return variable.
    if (exists $self->{_ca}) {
        # Here we can use the new system.
        @retVal = $self->{_ca}->DeleteMatchingAttributes($objectID, $key, @values);
    } else {
        # Confess("delete_matching_attributes not supported in old code.");
	if ($self->delete_attribute($objectID, $key, @values)) {
		return ($objectID, $key, @values);
	}
    }
    # Return the result.
    return @retVal;
}

=head3 change_attribute

    $fig->change_attribute($objectID, $key, \@oldValues, \@newValues);

Change the value of an attribute key/value pair for an object. This is
implemented as a delete followed by an insert.

=over 4

=item objectID

ID of the genome or feature to which the attribute is to be changed. In general, an ID that
starts with C<fig|> is treated as a feature ID, and an ID that is all digits and periods
is treated as a genome ID. For IDs of other types, this parameter should be a reference
to a 2-tuple consisting of the entity type name followed by the object ID.

=item key

Attribute key name. This corresponds to the name of a field in the database.

=item oldValues

One or more values identifying the key/value pair to change.

=item newValues

One or more values to be put in place of the old values.

=back

=cut

sub change_attribute {
    my ($self, @request) = @_;
    if (exists $self->{_ca}) {
        # Here we can use the new system.
        return $self->{_ca}->ChangeAttribute(@request);
    } else {
        my($peg,$k,$oldval, $oldurl, $newval, $newurl, $notl) = @request;
        return (0) unless ($peg && $k); # we must have at least a peg and a key.
        $k =  $self->clean_attribute_key($k);
        my $rdbH = $self->db_handle;

        # An error check to make sure that we are operating on the new version of attributes
        # If we are not, we will print an error and then return. Otherwise continue
        eval {$rdbH->SQL("SELECT genome,ftype,id,tag,val,url FROM attribute LIMIT 1")};
        if ($@) {return []}
        #if ($@) {print STDERR "Please rerun load_attributes to install the newest set of attributes\n"; return []}

        # Build the delete statement "@boundValues" will be the values replacing the
        # parameter marks.
        my $exc="DELETE FROM attribute WHERE ";
        my @boundValues;
        my ($delgenome, $delftype, $delid, $deltag)=($self->split_attribute_oid($peg), $k);
        $delgenome && ($exc .= "genome = ? and ") && (push @boundValues, $delgenome);
        $delftype && ($exc .= "ftype = ? and ") && (push @boundValues, $delftype);
        $delid && ($exc .= "id = ? and ") && (push @boundValues, $delid);
        $deltag && ($exc .= "tag = ? and ") && (push @boundValues, $deltag);
        $exc =~ s/and\s+$//;

        if ($oldval) {
            $exc .= " and val = ?";
            push @boundValues, $oldval;
            if ($oldurl) {
                $exc .= " AND url = ?";
                push @boundValues, $oldurl;
            }
        }
        $rdbH->SQL($exc, undef, @boundValues);
        if (defined $newval) {
            $exc = "INSERT INTO attribute (  genome,ftype,id,tag,val,url ) VALUES ( ?,?,?,?,?,? )";
            $rdbH->SQL($exc, undef, $self->split_attribute_oid($peg), $k, $newval, $newurl);
            # write to the transaction log if we add a new value (writing deletes is handled above)
            my $location = $self->attribute_location($peg);
            &verify_dir("$location");
            if (!$notl && open(TMPATTR,">>$location/transaction_log"))
            {
                print TMPATTR "CHANGE\t$peg\t$oldval\t$oldurl\t$newval\t$newurl\n";
                close(TMPATTR);
            }
        }
        return 1;
    }
}

=head3 clean_attribute_key()

## DEPRECATED ## This process is no longer required in the new system.

use $key=$fig->clean_attribute_key($key)

Keys for attributes are used as filenames in the code, and there are limitations on the characters that can be used in the key name. We provide an extended explanation of each key, so the key does not necessarily need to be person-readable.

Keys are not allowed to contain any non-word character (i.e. they must only contain [a-zA-Z0-9] and _

This method will remove these.

=cut

sub clean_attribute_key {
 my ($self, $key)=@_;
 #$key =~ s/[\s\n\t\$\@\/\\\Q!#%^&*()`~{}[]|:;"'<>?,.\E]//g; # the \Q .. \E just allows not escaping all the intermediate metacharacters
 my $old = $key;
 $key =~ s/\s+/\_/g;
 $key =~ s/\-/\_/g;
 $key =~ s/\W//g;
 $key =~ s/\_+/\_/g;
 return $key;
}

=head3 essential

    my $flag = $fig->essential($fid);

Return TRUE if a feature is considered essential and FALSE otherwise. This method
provides a uniform method for determining essentiality that will remain consistent
during the various overhauls of essentiality. Currently a feature is essential
if it has an attribute with the value C<essential> or C<potential_essential>.

=over 4

=item fid

ID of the feature to check for essentiality.

=item RETURN

Returns TRUE if the feature is considered essential, else FALSE.

=back

=cut

sub essential {
    # Get the parameters.
    my ($self, $fid) = @_;
    # Declare the return variable. We assume FALSE until proven otherwise.
    my $retVal = 0;
    # Check for essentiality.
    my @essentials = $self->get_attributes($fid, undef, 'essential');
    if (@essentials) {
        $retVal = 1;
    } else {
        # Check for potential essentiality.
        my @potentials = $self->get_attributes($fid, undef, 'potential_essential');
        if (@potentials) {
            $retVal = 1;
        }
    }
    # Return the result.
    return $retVal;
}

=head3 virulent

    my $flag = $fig->virulent($fid);

Return TRUE if a feature is considered virulent and FALSE otherwise. This method
provides a uniform method for determining virulence that will remain consistent
during the various overhauls of virulence attributes. Currently a feature is virulent
if it has an attribute whose key begins with C<virulence_associated>.

=over 4

=item fid

ID of the feature to check for essentiality.

=item RETURN

Returns TRUE if the feature is considered essential, else FALSE.

=back

=cut

sub virulent {
    # Get the parameters.
    my ($self, $fid) = @_;
    # Declare the return variable. We assume FALSE until proven otherwise.
    my $retVal = 0;
    # Get all the attributes and filter for virulence.
    my @attributes = $self->get_attributes($fid);
    # We loop until we prove virulence or run out of attributes.
    while (! $retVal && scalar(@attributes) > 0) {
        my $attributeThing = pop @attributes;
        # Each attribute entry is a 4-tuple. The key name is the second element.
        if ($attributeThing->[1] =~ /^virulence_associated/i) {
            $retVal = 1;
        }
    }
    # Return the result.
    return $retVal;
}

=head2 Splitting and Joining Attributes "oids"

There was a big problem with attributes being very slow to recover, and having to recover all attributes just to get those for a peg or a genome. The current implementation splits the original ID (oid) into three columns, genome, ftype, and id. The ftype is peg, rna, pp, etc. The id is the feature number. The genome is the genome number.

Hence:
fig|83333.1.peg.1345 becomes 83333.1, peg, and 1345
83333.1 becomes 83333.1, '', and ''

To split an oid into an array with three parts:
        $self->split_attribute_oid($peg);

To join the three parts of a series of results:
map {unshift @$_, $self->join_attribute_oid(splice(@$_, 0, 3))} @$res;

This code splices the first three elements of the the array, joins them, and then unshifts the result of that join back into the start of the array. Cool, eh?

=head3 split_attribute_oid()

use my ($genome, $type, $id)=split_attribute_feature($id);

splits an id into genome, type, and id if it is a feature, or just genome and '', '' if it is a genome, and just the id and undef undef if it is not known

=cut

sub split_attribute_oid {
 my($self, $id)=@_;
 if ($id =~ /^\d+\.\d+$/)
 {
  # it appears to be a genome id
  return ($id, "", "");
 }
 elsif ($id =~ /^fig\|(\d+\.\d+)\.(\w+)\.(\d+)/)
 {
  # it appears to be a feature
  return ($1, $2, $3);
 }
 else
 {
  # not sure what it is
  return ($id, "", "");
 }
}

=head3  join_attribute_oid()

use my $id=join_attribute_oid($genome, $feature, $id);

Joins an attribute back together after it has been pulled from the mysql database

=cut

sub join_attribute_oid {
 my ($self, @parts)=@_;
 if (($parts[0] && $parts[0] =~ /^\d+\.\d+$/) && ($parts[1] && $parts[1] =~ /^\w+$/) && ($parts[2] && $parts[2] =~ /^\d+$/))
 {
  # it is a feature ID
  return "fig|$parts[0].$parts[1].$parts[2]";
 }
 elsif ($parts[0] =~ /^\d+\.\d+$/ && !($parts[1] && $parts[2]))
 {
  # it is a genome
  return $parts[0];
 }
 else
 {
  return join("", @parts);
 }
}

=head3 read_attribute_transaction_log

use: $fig->read_attribute_transaction_log($logfile);

This method reads the transaction_log described in $logfile and enacts the changes described therein. The changes must be one of add, delete, or change.

=cut

sub read_attribute_transaction_log {
 my ($self, $file)=@_;
 return unless (-e $file);
 open(IN, $file) || die "Can't open $file";
 while (<IN>) {
  next if (/^\s*\#/ || /^\s*$/);
  chomp;
  my @line=split /\t/;
  next if ($line[2] eq "evidence_code");
  my $type=shift @line;
  if (uc($type) eq "DELETE")
  {
   $line[4]=1;
   $self->delete_attribute(@line);
  }
  elsif (uc($type) eq "ADD")
  {
   # some of the adds are lik this
   #print TMPATTR "ADD\t$peg\t$k\t$v\t$url\n";
   # and some are like this;
   #print TMPATTR "ADD\t$peg\t$k\t\t\t$v\t$url\n";
   # the second is the correct format
   if ($#line >= 4 && !($line[2]) && !($line[3])) {splice(@line, 2, 2)}
   $line[4]=1;
   $self->add_attribute(@line);
  }
  elsif (uc($type) eq "CHANGE")
  {
   $line[7]=1;
   $self->change_attribute(@line);
  }
  else
  {
   print STDERR "Do not understand this line from $file. It doesn't appear to be a transaction record:\n$_\n";
  }
 }
}

=head3 erase_attribute_entirely

This method will remove any notion of the attribute that you give it. It is different from delete as that just removes a single attribute associated with a peg. This will remove the files and uninstall the attributes from the database so there is no memory of that type of attribute. All of the attribute files are moved to FIG_Tmp/Attributes/deleted_attributes, and so you can recover the data for a while. Still, you should probably use this carefully!

I use this to clean out old PIR superfamily attributes immediately before installing the new correspondence table.

e.g. my $status=$fig->erase_attribute_entirely("structure");

This will return the number of files that were moved to the new location

=cut

sub erase_attribute_entirely {
 my ($self, $attr)=@_;
 return 0 unless ($attr);
 if (exists $self->{_ca}) {
    # Here we can use the new system.
    return $self->{_ca}->EraseAttribute($attr);
 } else {
    my %path_to_files; # this hash has the path as the key and the genome id as the value

    # first, find all the features with our attribute
    foreach my $attributes ($self->get_attributes(undef, $attr))
    {
      unless ($attributes->[1] eq $attr)
      {
       print STDERR "Warning : expected to erase $attr but we retrieved ", $attributes->[1], "\n";
       next;
      }
      #print STDERR "deleting attr: ", join(" ", @$attributes), "\n";
      $self->delete_attribute(@$attributes, 1);
      $path_to_files{$self->attribute_location($attributes->[0])}=$self->genome_of($attributes->[0]);
    }

    # now we need to check that we have the files to delete
    # we are going to see if there are files to delete, and then we will make temp dirs and move them. If there are no files
    # to do we don't need to make the dir
    my @files;
    foreach my $path (keys %path_to_files)
    {
     if (-e "$path/$attr") {push @files, $path}
    }

    return 1 unless (scalar @files); # don't continue if there are no files to move

    $self->verify_dir("$FIG_Config::temp/Attributes/deleted_attributes");

    foreach my $path (@files)
    {
     my $genome=$path_to_files{$path};
     unless ($genome) {$genome='unknown'}
     my $dest="$FIG_Config::temp/Attributes/deleted_attributes/$genome";
     mkdir "$FIG_Config::temp/Attributes/deleted_attributes/$genome", 0755 unless (-e "$FIG_Config::temp/Attributes/deleted_attributes/$genome");
     $dest .= "/".$attr;
     if (-e $dest)
     {
      # don't overwrite the file
      my $count=1;
      while (-e "$dest.$count") {$count++}
      $dest .= ".$count";
     }

     system("mv $path/$attr $dest");
    }

    return scalar @files;
 }
}

=head3 get_group_keys

    my @keys = $fig->get_group_keys($groupName);

Return all the attribute keys in the named group.

=over 4

=item groupName

Name of the group whose keys are desired.

=item RETURN

Returns a list of the attribute keys in the named group.

=back

=cut

sub get_group_keys {
    # Get the parameters.
    my ($self, $groupName) = @_;
    # Declare the return variable.
    my @retVal = $self->{_ca}->GetAttributeKeys($groupName);
    # Return the results.
    return @retVal;
}

=head3 get_group_key_info

    my %keys = $fig->get_group_key_info($groupName);

Return the descriptive data for all the attribute keys in the named group.

=over 4

=item groupName

Name of the group whose keys are desired. If omitted, then all keys will be returned. This
could be expensive, but when it's necessary, it's necessary.

=item RETURN

Returns a hash mapping each relevant attribute key to an n-tuple containing the the attribute relation name,
the description, and the 0 or more group names.

=back

=cut

sub get_group_key_info {
    # Get the parameters.
    my ($self, $groupName) = @_;
    # Declare the return variable.
    my %retVal;
    # Check the parameter.
    if (defined $groupName) {
        # Get all keys in the group.
        %retVal = $self->{_ca}->GetAttributeData(group => $groupName);
    } else {
        # Get all the keys.
        %retVal = $self->{_ca}->GetAttributeData(name => '');
    }
    # Return the results.
    return %retVal;
}

=head3 get_genome_keys

Get all the keys that apply to genomes and only genomes.
This method takes no arguments and returns an array.

=cut

sub get_genome_keys {
 my($self)=@_;
 if (exists $self->{_ca}) {
    return $self->{_ca}->GetAttributeKeys('Genome');
 } else {
    my $rdbH = $self->db_handle;

    # An error check to make sure that we are operating on the new version of attributes
    # If we are not, we will print an error and then return. Otherwise continue
    eval {$rdbH->SQL("SELECT genome,ftype,id,tag,val,url FROM attribute LIMIT 1")};
    if ($@) {return []}
    #if ($@) {print STDERR "Please rerun load_attributes to install the newest set of attributes\n"; return []}

    my $res=$rdbH->SQL("SELECT DISTINCT tag from attribute where (genome is not null and ftype = '' and id = '')");
    return map {$_->[0]} @$res;
 }
}

=head3 get_peg_keys

Get all the keys that apply just to pegs.
This method takes no arguments and returns an array.

=cut

sub get_peg_keys {
    my( $self ) = @_;
    if (exists $self->{_ca})
    {
        return $self->{_ca}->GetAttributeKeys('peg');
    }

    #  Add caching.  This saves time when called twice in one script.  This is
    #  done in a manner that immediately extends to other Feature types. -- GJO

    my $keys_hash = $self->cached( '_attribute_keys' );
    my $ans = $keys_hash->{peg};  #  Feature type

    if ( ! $ans )
    {
        my $rdbH = $self->db_handle;

        # An error check to make sure that we are operating on the new version of attributes
        # If we are not, we will print an error and then return. Otherwise continue

        eval { $rdbH->SQL("SELECT genome,ftype,id,tag,val,url FROM attribute LIMIT 1") };
        if ( $@ )
        {
            # print STDERR "Please rerun load_attributes to install the newest set of attributes\n";
            return [];
        }

        my $res = $rdbH->SQL("SELECT DISTINCT tag FROM attribute WHERE (ftype = 'peg')");
        $keys_hash->{peg} = $ans = [ map { $_->[0] } @$res ];
    }

    return @$ans;
}

=head3 get_peg_keys_for_genome

Get all the keys that apply just to pegs from a specified genome.
This method takes a genome id as an argument and returns an array.

=cut

sub get_peg_keys_for_genome {
    my ($self, $genome)=@_;
    if (exists $self->{_ca}) {
        my @list1 = $self->{_ca}->GetAttributes($genome);
        my @list2 = $self->{_ca}->GetAttributes(['Feature', "fig|$genome.%"]);
        return @list1, @list2;
    } else {
        my $rdbH = $self->db_handle;

        # An error check to make sure that we are operating on the new version of attributes
        #  # If we are not, we will print an error and then return. Otherwise continue
        eval {$rdbH->SQL("SELECT genome,ftype,id,tag,val,url FROM attribute LIMIT 1")};
        if ($@) {return []}
        #if ($@) {print STDERR "Please rerun load_attributes to install the newest set of attributes\n"; return []}

        my $res=$rdbH->SQL("SELECT genome,ftype,id,tag,val,url from attribute where (genome = '$genome' and ftype = 'peg')");

        # the following line takes the first 3 elements from each array and puts them back
        # # to be a feature or genome using join_attribute_oid and then puts them back in the array.
        map {unshift @$_, $self->join_attribute_oid(splice(@$_, 0, 3))} @$res;
        return @{$res};
    }
}

=head3 get_genomes_with_attribute

Get a list of all genomes that have a specified attribute. This will search for all genomes that have some attribute.

This will also accept partial matches. Hence to find all genomes that have essentiality data you can do this:

my @genomes=$fig->get_genomes_with_attribute("essential");

This will find Essential_Gene_Sets_Bacterial, essential, etc

=cut

sub get_genomes_with_attribute {
    my ($self, $attr)=@_;
    if (exists $self->{_ca}) {
        my @attributes = $self->{_ca}->GetAttributes(undef, "%$attr%");
        my %retVal = ();
        for my $attribute (@attributes) {
            if ($attribute->[0] =~ /^fig\|(\d+\.\d+)/ ||
                $attribute->[0] =~ /^(\d+\.\d+)/) {
                $retVal{$1} = 1;
            }
        }
        return sort keys %retVal;
    } else {
        my $rdbH = $self->db_handle;
        # An error check to make sure that we are operating on the new version of attributes
        # do we still need this? Probably.
        eval {$rdbH->SQL("SELECT genome,ftype,id,tag,val,url FROM attribute LIMIT 1")};
        if ($@) {return []}

        my $res=$rdbH->SQL("SELECT distinct genome from attribute where (tag like '\%$attr\%')");
        return map {$_->[0]} @$res;
    }
}



=head3 key_info

DEPRECATED: in actual fact, no attribute metadata was ever put into the system.

Access a hash of key information. The data that are returned are currently:

hash key name           what is it                      data type
single                                                  [boolean]
description             Explanation of key              [free text]
readonly                whether to allow read/write     [boolean]
is_cv                   attribute is a cv term          [boolean]

Single is a boolean, if it is true only the last value returned should be used. Note that the other methods willl still return all the values, it is upto the implementer to ensure that only the last value is used.

Explanation is a user-derived explanation that can be free text

If a reference to a hash is provided, along with the key, those values will be set to the attribute_keys file

Returns an empty hash if the key is not provieded or doesn't exist

e.g.
$fig->key_info($key, \%data); # set the data
$data=$fig->key_info($key); # get the data

This data is stored in a file called $FIG_Config::global/Attributes/attribute_metadata and in a database called attribute_metadata. The data is strictly on a last in last out basis, so that if a datapoint is changed, the last datapoint in the database or file is returned. At the moment I am not coding the ability to edit data.

The method takes the following arguments

=over 4

=item key

The key to look for or add data to.

=item $data

A reference to a hash containing the new data to add to the database. If provided this will cause the database to be updated

=item $nowrite

Do not write the new data to the attributes_metadata file. This is mainly used by load_attributes to prevent a circular read/write condition.

=back

=cut

sub key_info {
 my ($self, $key, $data, $nowrite)=@_;
 return ();
 #
 #return unless ($key);
 #$key =  $self->clean_attribute_key($key);
 #my $rdbH = $self->db_handle;
 #
 ## An error check to make sure that we are operating on the new version of attributes
 ## If we are not, we will print an error and then return. Otherwise continue
 #eval {$rdbH->SQL("SELECT genome,ftype,id,tag,val,url FROM attribute LIMIT 1")};
 #if ($@) {return []}
 ##if ($@) {print STDERR "Please rerun load_attributes to install the newest set of attributes\n"; return []}
 #
 #unless ($data)
 #{
 # # we can just return the info right away
 # return $self->{'key_info'}->{$key} if ($self->{'key_info'}->{$key});
 # my $res=$rdbH->SQL("SELECT  metakey, metaval from attribute_metadata where attrkey = ?", undef, $key);
 # foreach my $result (@$res)
 # {
 #  $self->{'key_info'}->{$key}->{$result->[0]}=$result->[1];
 # }
 # return $self->{'key_info'}->{$key};
 #}
 #
 ## there is new data to add
 ## first, check if we have an old style attributes file and update it. eventually we should be able to delete this line.
 #if (-e "$FIG_Config::global/Attributes/attribute_keys") {$self->update_attributes_metadata("$FIG_Config::global/Attributes/attribute_keys")}
 #
 ## now append the new data to the attributes_metadata file
 #unless ($nowrite)
 #{
 # open (OUT, ">>$FIG_Config::global/Attributes/attribute_metadata") || die "Can't append to $FIG_Config::global/Attributes/attribute_metadata";
 #}
 #foreach my $datum (keys %$data)
 #{
 # unless (defined $data->{$datum}) {$data->{$datum}='true'} # just make it true so that it exists
 # unless ($nowrite) {print OUT "$key\t$datum\t", $data->{$datum}, "\n"}
 #
 # $rdbH->SQL("INSERT INTO attribute_metadata (attrkey, metakey, metaval) VALUES (?,?,?) ", undef, $key, $datum, $data->{$datum});
 #}
 #unless ($nowrite) {close OUT}
 #my $res=$rdbH->SQL("SELECT  metakey, metaval from attribute_metadata where attrkey = ?", undef, $key);
 #foreach my $result (@$res)
 #{
 # $self->{'key_info'}->{$key}->{$result->[0]}=$result->[1];
 #}
 #return $self->{'key_info'}->{$key};
}

=head3 update_attributes_metadata()

This method exists solely to update the attributes metadata file and make sure that it is in the right format.
This method can probably be deleted in a while, but it needs to be run on all machines with attributes data before then!

It is only called if an old attributes metadata file is found.

The method returns the filename where the data is now stored.

=cut


sub update_attributes_metadata {
 my ($self, $file)=@_;
 my $version=1;
 my $attr;
 open(IN, $file) || die "Can't open $file for reading";
 while (<IN>) {
  if (/^\#\s*Version\s*(\d+)/) {$version=$1}
  next if (/^\s*\#/);
  chomp;
  next unless ($_);
  my @a=split /\t/;
  # fix old versions of attribute_keys
  if ($version==1) {$attr->{$a[0]}->{'single'}=$a[1]; $attr->{$a[0]}->{'description'}=$a[2]; next}
  $attr->{$a[0]}->{$a[1]}=$a[2];
 }
 close IN;
 unlink($file);

 my $rdbH = $self->db_handle;

 # An error check to make sure that we are operating on the new version of attributes
 # If we are not, we will print an error and then return. Otherwise continue
 eval {$rdbH->SQL("SELECT genome,ftype,id,tag,val,url FROM attribute LIMIT 1")};
 if ($@) {return []}
 #if ($@) {print STDERR "Please rerun load_attributes to install the newest set of attributes\n"; return []}

 unless (-e "$FIG_Config::global/Attributes/attribute_metadata")
 {
  open (OUT, ">$FIG_Config::global/Attributes/attribute_metadata") || die "Can't open $FIG_Config::global/Attributes/attribute_metadata";
  print OUT "# Version 2\n# This file contains information about the attribute keys in this database. The columns are:\n";
  print OUT "# attribute key\n# tag associated for that key\n# value of that tag\n";
  print OUT "# Each attribute key can have as many of these as you want. The last one in the file will be used. This is used to store data applicable to\n";
  print OUT "# every key in the attributes\n";
  close OUT;
 }

 open (OUT, ">>$FIG_Config::global/Attributes/attribute_metadata") || die "Can't open $FIG_Config::global/Attributes/attribute_metadata";
 foreach my $keyName (keys %$attr) {
  foreach my $attrName (keys %{$attr->{$keyName}} ) {
   unless (defined $attr->{$keyName}->{$attrName}) {$attr->{$keyName}->{$attrName}=1}
   print OUT "$keyName\t$attrName\t", $attr->{$keyName}->{$attrName}, "\n";
   my $res=$rdbH->SQL("INSERT INTO attribute_metadata (attrkey, metakey, metaval) VALUES (?,?,?)", undef, $keyName, $attrName, $attr->{$keyName}->{$attrName});
  }
 }
 close OUT;
 return "$FIG_Config::global/Attributes/attribute_metadata";
}

=head3 get_values

Get all the values that we know about

Without any arguments:

Returns a reference to a hash, where the key is the type of feature (peg, genome, rna, prophage, etc), and the value is a reference to a hash where the key is the value and the value is the number of occurences

e.g. print "There are  " , {$fig->get_values}->{'peg'}->{'100'}, " keys with the value 100 in  the database\n";

With a single argument:

The argument is assumed to be the type (rna, peg, genome, etc).

With two arguments:

The first argument is the type (rna, peg, genome, etc), and the second argument is the key.

In each case it will return a reference to a hash.
E.g.

        $fig->get_values(); # will get all values

        $fig->get_values('peg'); # will get all values for pegs

        $fig->get_values('peg', 'structure'); # will get all values for pegs with attribute structure

        $fig->get_values(undef, 'structure'); # will get all values for anything with that attribute

=cut


sub get_values {
 my ($self, $want, $tag)=@_;
 unless ($want) {$want="all"}
 my $rdbH = $self->db_handle;
 $tag =~ s/^\s+//; $tag =~ s/\s+$//; $tag=uc($tag);

 my $sql="SELECT genome, ftype, id ,val FROM attribute";
 if ($tag) {$sql .= " WHERE tag = \'$tag\'"}

 my $relational_db_response=$rdbH->SQL($sql);

 my $tags;
 foreach my $res (@$relational_db_response) {
  my ($fid, $val)=($self->join_attribute_oid(splice(@$res,0,3)), $res->[$#$res]);
  my $type=$self->ftype($fid);
  if ($type && ($want eq $type || $want eq "all")) {
   $tags->{$type}->{$val}++;
  } elsif (($fid =~ /^\d+\.\d+$/) && (lc($want) eq "genome" || $want eq "all")) {
   $tags->{'genome'}->{$val}++;
  }
 }
 if ($want eq "all") {return $tags} else {return $tags->{$want}}
}




=head3 guess_value_format

There are occassions where I want to know what a value is for a key. I have three scenarios right now:

 1. strings
 2. numbers
 3. percentiles ( a type of number, I know)

In these cases, I may want to know something about them and do something interesting with them. This will try and guess what the values are for a given key so that you can try and limit what people add. At the moment this is pure guess work, although I suppose we could put some restrictions on t/v pairs I don't feel like.

This method will return a reference to an array. If the element is a string there will only be one element in that array, the word "string". If the value is a number, there will be three elements, the word "float" in position 0, and then the minimum and maximum values. You can figure out if it is a percent :-)

=cut

sub guess_value_format {
 my ($self, $tag)=@_;
 return unless ($tag);

 # I am using self->{'value_format'} to save the format so if this is called multiple times it is not recalculated each time
 return $self->{'value_format'}->{$tag} if (defined $self->{'value_format'}->{$tag});

 my $hash = $self->get_values(undef, $tag);
 return if (!$hash || !scalar keys %$hash); # don't carry on if there is nothing to look at

 my ($min, $max)=(100000000, 0);
 foreach my $type (keys %$hash) {
  foreach my $val (keys %{$hash->{$type}}) {
    next unless ($val);

    # this code is taken from the perl cookbook pg 44
    # it should detect for all nummbers
    if ($val !~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/) {
        undef $min;
        undef $max;
        last;
    }
    else {
      if ($val > $max) {$max=$val}
      if ($val < $min) {$min=$val}
    }
   }
 }
 # if $min and $max are defined then the value is a number
 # if not, then it is a string;
 if (defined $min && defined $max) {$self->{'value_format'}->{$tag} = ["float", $min, $max]} else {$self->{'value_format'}->{$tag}=["string"]}
 return $self->{'value_format'}->{$tag};
}

=head3 attribute_location

This is just an internal method to find the appropriate location of the attributes file depending on whether it is a peg, an rna, or a genome or whatever.

=cut

sub attribute_location {
    my ($self, $peg)=@_;
    return unless ($peg);
    my $type=$self->ftype($peg); # we need to know whether it is a peg, prophage, rna, etc

    my $location;
    if ($type)
    {
     my $genome = &genome_of($peg);
     $location="$FIG_Config::organisms/$genome/Features/$type/Attributes";
    }
    elsif ($peg =~ /^\d+\.\d+$/ && (-e "$FIG_Config::organisms/$peg"))
    {
     # $peg is a genome number and we know about it
     #$location="$FIG_Config::organisms/$peg/Attributes";
     # we want to put things in global again
     $location="$FIG_Config::global/Attributes/";
    }
    elsif (lc($peg) eq "subsystem")
    {
        $location="$FIG_Config::global/Attributes/";
    }
    else {
     print STDERR "Can't figure out what $peg is. It is neither a known feature or a genome id. Added to $FIG_Config::global/Attributes/\n";
     $location="$FIG_Config::global/Attributes/";
    }

    return $location;
}


=head3 add_cv_term

Add a controlled vocabulary term to a peg.  Pass in the peg, the vocab
name, the termId, and the term (see next paragraph).  returns error string
if problem, else returns nothing.

   my $status = $fig->add_cv_term( "master:EdF",
                                   "fig|9606.3.peg.26823", "MyVocab", "1234", "A thing of wonder.");
   if ($status) {print "error adding cv term: $status\n";}

Controlled vocabulary is read-only text associated with a peg.  Each
is a triple, namely (vocab name, termId, term text).  The termId is an
id that is used in the particulary vocabulary and the term text is the
actual term.  For example, the GO has the term "U12-type nuclear mRNA
branch site recognition" with termId GO:0000371.  Thus, the triplet is
(GO, GO:0000371, "U12-type nuclear mRNA branch site recognition").
Don't be confused by the GO: in GO:0000371.  We don't add the GO:.
That's just what GO decided to do.

termIds can not have ';' in them.

This routine encapsulates our present implementation via attributes.

=cut

sub add_cv_term {

    # $user not used yet but maybe should track who's doing the cv add?
    my ($self, $user, $peg, $vocab, $term_id, $term)=@_;

    $user =~ s/^\s*//g;
    $user =~ s/\s*$//g;
    $peg =~ s/^\s*//g;
    $peg =~ s/\s*$//g;
    $vocab =~ s/^\s*//g;
    $vocab =~ s/\s*$//g;
    $term_id =~ s/^\s*//g;
    $term_id =~ s/\s*$//g;
    $term =~ s/^\s*//g;
    $term =~ s/\s*$//g;


    if ( ! ($user && $peg && $vocab && $term_id && $term) ) {
        #print STDERR "add_cv_term: invalid arguments. All required, no empty strings\n";
        return "add_cv_term: invalid arguments. All required, no empty strings";
    }


    # make sure the key (the vocab name) is flagged as read only
    # and as CV in the global attributes meta data.  don't set this
    # if it's already there because otherwise a write of a new file
    # happens

    #key info returns empy hash rather than undef
    #orig my %key_info_hash = $self->key_info($vocab);

    my %key_info_hash;

    if ($self->key_info($vocab) ) {
        %key_info_hash = %{$self->key_info($vocab)};
    }

    if  ( ! keys %key_info_hash ) {
        print STDERR "add_cv_term: New CV name $vocab being added to key_info\n";

        $key_info_hash{"single"} = 0;
        $key_info_hash{"description"} = "Controled Vocabulary, $vocab";
        $key_info_hash{"is_cv"} = 1;
        $key_info_hash{"single"} = 0;
        $key_info_hash{"readonly"} = 1;

        $self->key_info($vocab,\%key_info_hash);

    } else {
        if (! $self->key_info($vocab)->{"is_cv"} ) {
            print STDERR "Error: attempt to use existing, non CV key as CV name\n";
            return "add_cv_term: Error: attempt to use existing, non CV key as CV name: $vocab";
        }
        print STDERR "add_cv_term: reusing existing CV name $vocab\n";
    }


    # shove in the attribuute
    #
    # we combine the term ID and term text into the attribute value
    # separated by a ";" which is forbidden from the termID (but not
    # from the term itself.

    my $x = $term_id . "; " . $term;

    $self->add_attribute($peg,$vocab, $x);
    return;
}

sub search_index_by_attribute {
    # please don't put a method between its description and the method. Honor the docs that we have.
    # Please add pod for these methods, too.

    # supports search_index method by finding attributes via the attribute table in
    # the database rather than via glimpse indexes.  This will go away with Bobs
    # migration to Lucene, but for now we've been asked to give immediate search
    # capability on attributes without rerunning index building.
    #
    # return array of (peg, org, aliasList, function) where we'll set aliasList to
    # the value of the alias and leave function blank.
    #
    # now case _in_sensitive
    #
    my($self,$searchTerm)=@_;
    return unless( $searchTerm);
    my @results;
    if ($self->{_ca}) {
        # Here we're using the new attribute system.
    } else {
        my $rdbH = $self->db_handle;

        # An error check to make sure that we are operating on the new version of attributes
        # If we are not, we will print an error and then return. Otherwise continue
        eval {$rdbH->SQL("SELECT genome,ftype,id,tag,val,url FROM attribute LIMIT 1")};
        if ($@) {return []}
        #if ($@) {print STDERR "Please rerun load_attributes to install the newest set of attributes\n"; return []}

        my $theTerm = uc( $searchTerm );
        my $relational_db_response=$rdbH->SQL("SELECT genome,ftype,id,tag,val from attribute WHERE UPPER(tag) LIKE '%$theTerm%' OR UPPER(val) LIKE '%$theTerm%'");

        my @results;
        foreach my $res (@$relational_db_response) {
            my ($genome,$ftype,$id, $tag, $value)=@$res;
            my $fid=$self->join_attribute_oid($genome,$ftype,$id);
            my $org = $self->genus_species( $self->genome_of($fid) );
            my @aliases = $self->feature_aliases($fid);
            my $a_string =join(" ",@aliases);
            #my $a_string = "test";
            push (@results, [$fid, $org, $a_string,"[attribute $tag] $value", $genome] );
            #the prior way - modified to accomodate consistent format of downloaded results
            #push (@results, [$fid, $org, "[attribute $tag] $value",""] );
        }
    }
    return @results;
}


sub find_by_attribute {
    # search by substrings in attribute values or attribute tags.
    # This might replace the present search-for-attributes that works by
    # glimpse.  The problem with the present approach is that you can't
    # search until you rebuild indices with make_attribute_index
    #


    my($self,$searchTerm)=@_;
    return unless( $searchTerm);
    my $rdbH = $self->db_handle;

    # An error check to make sure that we are operating on the new version of attributes
    # If we are not, we will print an error and then return. Otherwise continue
    eval {$rdbH->SQL("SELECT genome,ftype,id,tag,val,url FROM attribute LIMIT 1")};
    if ($@) {return []}
    #if ($@) {print STDERR "Please rerun load_attributes to install the newest set of attributes\n"; return []}

    my $relational_db_response=$rdbH->SQL("SELECT genome,ftype,id,tag,val from attribute WHERE tag LIKE '%$searchTerm%' OR val LIKE '%$searchTerm%'");
    my @results;

    foreach my $res (@$relational_db_response) {
        my ($genome,$ftype,$id, $tag, $value)=@$res;
        my $fid=$self->join_attribute_oid($genome,$ftype,$id);
        push (@results, [$fid, $tag, $value]);
    }
    return @results;
}



=head3 search_cv_file

Search a controlled vocabulary file for desired text. Pass the
name of the CV, e.g., "GO" or "HUGO" and get back a reference
to a list of results.  Each result is a line from the file,
and so is a tab-separated representation of the tripilet,
(CV_name, CV_id, CV_text)

Case insensitivee, substring.
=cut

sub search_cv_file
{
    my ($self, $cv,$search_term) =@_;
    my $file = $FIG_Config::global."/CV/cv_search_".$cv.".txt";
    if (! open(LOOKUP,"$file") ) {
        print STDERR "Search could not find vocabulary file, $file\n";
        return;
    }
    my @lines;
    while (<LOOKUP>) {
        chomp;
        push @lines, $_;
    }

    my @grep_results = grep(/$search_term/i,@lines);
    return [@grep_results];
}



################################# Indexing Features and Functional Roles  ####################################

=head3 search_index

    my ($pegs,$roles) = fig->search_index($pattern, $non_word_search, $user);

Find all pegs and roles that match a search pattern. The syntax of I<$pattern>
is deliberately left undefined so that we can change the underlying technology, but
a single word or phrase should work.

=over 4

=item pattern

A search pattern. In general, the pattern is a single word or phrase that is expected
to occur somewhere in a functional role, attribute key, or attribute value.

=item non_word_search (optional)

If specified, the pattern will be interpreted as a string instead of a series of
words.

=item user (optional)

If specified, the name of the current user. That user's annotation will be given precedence
when the functional role is determined.

=item RETURN

Returns a 2-tuple. The first element is a reference to a list of features. For each
feature, there is a tuple consisting of the (0) feature ID, (1) the organism name (genus
and species), (2) the aliases, (3) the functional role, and (4) the relevant annotator. The
second element in the returned tuple is a reference to a list of functional roles. All
the roles and features in the lists must match the pattern in some way.

=back

=cut

sub search_index {
    # Get the parameters.
    my ($self, $pattern, $non_word_search, $user) = @_;
    # Clean up the temporary directory to insure there's room for search results.
    &clean_tmp;
    # Convert the search pattern to Glimpse format. First, we convert spaces to semicolons.
    my $patternQ = $pattern;
    $patternQ =~ s/\s+/;/g;
    # Stop here to extract the search terms.
    my @words = split /;/, $pattern;
    Trace("Word list = (" . join(", ", @words) . ")") if T(Glimpse => 3);
    # Now escape the periods.
    $patternQ =~ s/\./\\./g;
    # Compute the glimpse directory. This facility is provided for testing purposes only.
    # If a "glimpse" member is specified in FIG_Config, then it will be presumed to contain
    # glimpse indexes. Thus, we can load a test index into a separate directory and twiddle
    # FIG_Config so we can run against the test index.
    my $dirName = (defined($FIG_Config::glimpse) ? $FIG_Config::glimpse : "$FIG_Config::data/Indexes");
    # Format the glimpse options. This is where the "non_word_search" parameter
    # is incorporated.
    my $glimpse_args = "-y -H \"$dirName\" -i";
    $glimpse_args .= " -w" unless $non_word_search;
    $glimpse_args .= " \'$patternQ\'";
    Trace("Search pattern = \"$pattern\", normalized to \"$patternQ\".") if T(Glimpse => 3);
    Trace("Glimpse parameters are: $glimpse_args") if T(Glimpse => 3);
    Trace("Glimpse directory is $FIG_Config::ext_bin") if T(Glimpse => 3);
    # Get the raw glimpse output. We also keep the error output for tracing purposes.
    my $errorFile = "$FIG_Config::temp/glimpseErrors$$.log";
    my @raw = `$FIG_Config::ext_bin/glimpse $glimpse_args 2>$errorFile`;
    # my @raw = `$FIG_Config::ext_bin/glimpse $glimpse_args`;
    my $rawCount = @raw;
    if ($rawCount == 0) {
        # No lines returned, so trace the error lines.
        my $errors = Tracer::GetFile($errorFile);
        Trace("Error lines from Glimpse:\n$errors") if T(Glimpse => 3);
    } else {
        Trace("$rawCount lines returned from glimpse.") if T(Glimpse => 3);
    }
    # Extract the feature lines from the raw data.
    my @pegs  =  map { $_ =~ /^\S+:\s+(\S.*\S)/; [split(/\t/,$1)] }
              grep { $_ =~ /^\S+peg.index/ } @raw;
    # Create a hash to hold the PEG data found so far.
    my %pegsFound = ();
    # Put the pegs found so far into the hash.
    for my $rawTuple (@pegs) {
        # Get this peg's data.
        my ($peg, $gs, $aliases, @funcs) = @{$rawTuple};
        # Only proceed if the peg exists.
        if (! $self->is_deleted_fid($peg)) {
            # Clean the glimpse markers out of the aliases. While we're at it, make
            # sure we have a string instead of an undef.
            if ($aliases) {
                $aliases =~ s/^aliases://;
            } else {
                $aliases = "";
            }
            # Process the functional assignments. Some of these will actually be
            # attribute key-value pairs. We'll create one list for stashing functional
            # assignments, and another for stashing attribute data. Note that we'll
            # only keep attributes that match one of the search words.
            my @functionList = ();
            my @attributeList = ();
            for my $func (@funcs) {
                Trace("$peg Function: $func") if T(Glimpse => 4);
                if ($func =~ /^function:\s*(.+)#(.+)$/) {
                    # Here we have a functional assignment. We push it onto the
                    # function list in the form (user, function).
                    push @functionList, [$2,$1];
                } elsif ($func =~ /^attribute:\s*(.+)$/) {
                    # Here we have an attribute. We only care if one of our
                    # search terms is in it.
                    Trace("Attribute entry $func.") if T(Glimpse => 4);
                    my $attributeAssignment = $1;
                    my $found = grep { $attributeAssignment =~ /$_/i } @words;
                    if ($found) {
                        push @attributeList, $attributeAssignment;
                    }
                }
            }
            # Find the desired functional role.
            my ($who, $function) = $self->choose_function($user, @functionList);
            # Store this peg in the hash.
            $pegsFound{$peg} = [$gs, $aliases, $function, $who, join("; ", @attributeList)];
        }
    }
    my $pegCount = keys %pegsFound;
    Trace("Raw glimpse results processed. $pegCount pegs found.") if T(Glimpse => 3);
    # Now form the list of PEGs from the hash.
    @pegs = map { [$_, @{$pegsFound{$_}}] } sort { &FIG::by_fig_id($a,$b) } keys %pegsFound;
    # PEGs are done, now do the roles.
    my @rolesT = grep { $_ =~ /^\S+role.index/ } @raw;
    my %roles  = map { $_ =~ /^\S+:\s+(\S.*\S)/; $1 => 1;} @rolesT;
    my @roles  = keys(%roles);
    # Return both lists.
    return ([@pegs],[@roles]);
}

=head3 choose_function

    my ($who, $function) = $fig->choose_function($user, @funcs);

Choose the best functional role from a list of role/user tuples. If a user is
specified, we look for one by that user. If that doesn't work, we look for one
by a master user. If THAT doesn't work, we take the first one.

=over 4

=item user

The name of the current user. If no user is active, specify either C<undef> or
a null string.

=item funcs

List of functional roles. Each role is represented by a 2-tuple consisting of the
user name followed by the role description.

=back

=cut

sub choose_function {
    # Get the parameters.
    my ($self, $user, @funcs) = @_;
    # We'll store the best role in here.
    my $function;
    # This will be used as an array index.
    my $i;
    # Get the number of functions.
    my $funCount = @funcs;
    # If a user was specified, choose his first assignment.
    if ($user) {
        # Find the first functional role for this user.
        for ($i = 0; ($i < $funCount) && ($funcs[$i]->[0] !~ /^$user/i); $i++) {}
        Trace("I = $i") if T(4);
        if ($i < $funCount) {
            $function = $funcs[$i];
        }
    }
    # If we didn't have a user or didn't find an assignment for this user, look
    # for a master user.
    if (! $function) {
        for ($i = 0; ($i < $funCount) && ($funcs[$i]->[0] !~ /^master/i); $i++) {}
        if ($i < $funCount) {
            $function = $funcs[$i];
        }
    }
    # If we still don't have a function, and a function exists, take the first one.
    if (! $function) {
        if ($funCount > 0) {
            $function = $funcs[0];
        } else {
            # No hope, return an empty list.
            $function = [];
        }
    }
    # Return the function found.
    return @{$function};
}
################################# Loading Databases  ####################################
=head3 load_all_list

    my @packages = FIG::load_all_list();

Return a list of the commands to be executed in order to load the SEED database.

=cut

sub load_all_list {
    my @packages = qw(load_peg_mapping
              index_contigs
              compute_genome_counts
              load_features
              index_sims
              index_translations
              add_assertions_of_function
              load_protein_families
              load_external_orgs
              index_annotations
              init_maps
              load_kegg
              load_distances
              make_indexes
              format_peg_dbs
              load_links
	      load_dlits
              index_subsystems
              load_bbhs
              load_literature
	      load_locks
              load_coupling
              load_chromosomal_clusters
	      load_go
	      load_id_correspondence
           );
#              load_ec_names
#              load_pch_pins
#              index_neighborhoods

    push(@packages, "pegs_in_conflict | peg_to_subsystems > $FIG_Config::global/conflicted.pegs");
    return @packages;
}

#=pod
#
#=head3 load_all
#
#usage: load_all
#
#This function is supposed to reload all entries into the database and do
#whatever is required to properly support indexing of pegs and roles.
#
#=cut

sub load_all {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($skip_to) = @_;
    my $start;
    my $just_list;

    #
    # If skip_to is numeric, start with that package.
    #
    # If it is the string "list", list the packages with their numbers.
    #

    if ($skip_to eq "list")
    {
        $just_list = 1;
    }
    elsif ($skip_to =~ /^\d+$/)
    {
        $start = $skip_to - 1;
    }
    else
    {
        $start = 0;
    }

    Trace("Loading SEED data.") if T(2);

    my @packages = load_all_list;

    my $pn = @packages;
    for my $i ($start..@packages - 1)
    {
        my $i1 = $i + 1;
        my $pkg = $packages[$i];

        my $date = `date`;
        chomp $date;
        print "$date:  Running $pkg ($i1 of $pn)\n";

        if (!$just_list)
        {
            &run($pkg);
        }
    }
    print "\n\nLoad complete.\n\n";
}

################################# Automated Assignments  ####################################

=head3 auto_assign

usage: $assignment = &FIG::auto_assign($peg,$seq)

This returns an automated assignment for $peg.  $seq is optional; if it is not
present, then it is assumed that similarities already exist for $peg.  $assignment is set
to either

    Function
or
    Function\tW

if it is felt that the assertion is pretty weak.

=cut

sub auto_assign {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($peg,$seq) = @_;

    my $cmd = $seq ? "echo \"$peg\" \"$seq\" | $FIG_Config::bin/auto_assign | $FIG_Config::bin/make_calls" : "echo \"$peg\" | $FIG_Config::bin/auto_assign | $FIG_Config::bin/make_calls";
#   print STDERR $cmd;
    my(@tmp) = `$cmd`;
    if ((@tmp == 1) && ($tmp[0] =~ /^\S+\t(\S.*\S)/))
    {
        return $1;
    }
    else
    {
        return "hypothetical protein";
    }
}

################################# Protein Families ####################################

=head2 Protein Families

In the protein families we have our own concept of an id that I have called an cid. This is entirely internal and does not map to any known database except our own, however it is used to store the correspondence between different protein families. Therefore, to find out what family any protein is in you need to convert that protein to an cid. You can start with a KEGG, COG, TIGR, SP, GI, or FIG id, and get an cid back. From there, you can find out what other proteins that cid maps to, and what families that protein is also in.

=head3 all_protein_families

usage: @all = $fig->all_protein_families

Returns a list of the ids of all of the protein families currently defined.

=cut

sub all_protein_families {
    my($self) = @_;

    my($relational_db_response);
    my $rdbH = $self->db_handle;

    if ($rdbH->table_exists('localfam_function') &&
        ($relational_db_response = $rdbH->SQL("SELECT DISTINCT family FROM localfam_function")) &&
        (@$relational_db_response >= 1))
    {
        return map { $_->[0] } @$relational_db_response;
    }
    return ();
}

=head3 families_for_protein

    my @families = $fig->families_for_protein($peg);

Return a list of all the families containing the specified protein.

=over 4

=item peg

ID of the PEG representing the protein in question.

=item RETURN

Returns a list of the IDs of the families containing the protein.

=back

=cut

sub families_for_protein {
 my($self, $peg)=@_;
 return () unless ($peg);
 my $cid=$self->prot_to_cid($peg);
 return unless ($cid);
 return $self->in_family($cid);
}

=head3 proteins_in_family

    my @proteins = $fig->proteins_in_family($family);

Return a list of every protein in a family.

=over 4

=item family

ID of the relevant protein family.

=item RETURN

Returns a list of all the proteins in the specified family.

=back

=cut

sub proteins_in_family {
 my($self, $family)=@_;
 return () unless ($family);
 my @prots;
 foreach my $cid ($self->ids_in_family($family)) {
  push @prots, $self->cid_to_prots($cid);
 }
 # note that some proteins may be duplicated, so we flatten this array and return only those things that are unique
 my %seen; # only return the first occurence of anyting.
 return grep {!$seen{$_}++} @prots;
}

=head3 family_function

    my $func = $fig->family_function($family);

Returns the putative function of all of the pegs in a protein family.  Remember, we
are defining "protein family" as a set of homologous proteins that have the
same function.

=over 4

=item family

ID of the relevant protein family.

=item RETURN

Returns the name of the function assigned to the members of the specified family.

=back

=cut

sub family_function {
    my($self,$family) = @_;
    return "" unless ($family);
    my($relational_db_response);

    my $rdbH = $self->db_handle;

    if ($rdbH->table_exists('localfam_function') &&
        ($relational_db_response = $rdbH->SQL("SELECT function from localfam_function WHERE family = '$family'")) &&
        (@$relational_db_response >= 1))
    {
        return $relational_db_response->[0]->[0];
    }
    return "";
}

=head3 sz_family

    my $n = $fig->sz_family($family);

Returns the number of proteins in a family.

=over 4

=item family

ID of the relevant protein family.

=item RETURN

Returns the number of proteins in the specified family.

=back

=cut

sub sz_family {
    my($self,$family) = @_;
    return 0 unless ($family);

    my($relational_db_response);
    my $rdbH = $self->db_handle;

    if ($rdbH->table_exists('localfam_function') &&
        ($relational_db_response = $rdbH->SQL("SELECT DISTINCT cid from localfam_cid WHERE family = '$family'")))
    {
        return scalar @$relational_db_response;
    }
    return 0;
}

=head3 ext_sz_family

usage: $n = $fig->ext_sz_family($family)

Returns the number of external IDs in $family.

=cut

sub ext_sz_family {
    my($self,$family) = @_;
    return 0 unless ($family);
    my @proteins=$self->ext_ids_in_family($family);
    return scalar(@proteins);
}

=head3 all_cids

usage: @all_cids=$fig->all_cids();

Returns a list of all the ids we know about.

=cut

sub all_cids {
    my($self) = @_;

    my($relational_db_response);
    my $rdbH = $self->db_handle;

    if ($rdbH->table_exists('localfam_cid') &&
        ($relational_db_response = $rdbH->SQL("SELECT DISTINCT cid FROM localfam_cid")) &&
        (@$relational_db_response >= 1))
    {
        return map { $_->[0] } @$relational_db_response;
    }
    return ();
}

=head3 ids_in_family

usage: @pegs = $fig->ids_in_family($family)

Returns a list of the cids in $family.

=cut

sub ids_in_family {
    my($self,$family) = @_;
    return () unless ($family);

    my($relational_db_response);
    my $rdbH = $self->db_handle;

    if ($rdbH->table_exists('localfam_function') &&
        ($relational_db_response = $rdbH->SQL("SELECT DISTINCT cid from localfam_cid WHERE family = '$family'")) &&
        (@$relational_db_response >= 1))
    {
        return map { $_->[0] } @$relational_db_response;
    }
    return ();
}

=head3 in_family

usage: @families = $fig->in_family($cid)

Returns an array containing the families containing an cid.

=cut

sub in_family {
    my($self,$cid) = @_;
    return () unless ($cid);

    my($relational_db_response);
    my $rdbH = $self->db_handle;

    if ($rdbH->table_exists('localfam_function') &&
        ($relational_db_response = $rdbH->SQL("SELECT DISTINCT family from localfam_cid WHERE cid = $cid")))
    {
     my %seen; # only return the first occurence of anyting.
     return grep {!$seen{$_}++} map { $_->[0] } @$relational_db_response;
    }
    return ();
}


=head3 ext_ids_in_family

usage: @exts = $fig->ext_ids_in_family($family)

Returns a list of the external ids in an external family name.

=cut

sub ext_ids_in_family {
    my($self,$family) = @_;
    return () unless ($family);

    my($relational_db_response);
    my $rdbH = $self->db_handle;

    if ($rdbH->table_exists('localid_map') &&
        ($relational_db_response = $rdbH->SQL("SELECT DISTINCT localid from localid_map WHERE family = '$family'")) &&
        (@$relational_db_response >= 1))
    {
        return map { $_->[0] } @$relational_db_response;
    }
    return ();
}

=head3 ext_in_family

usage: @ext_families = $fig->ext_in_family($id)

Returns an array containing the external families containing an id. The ID is the one from the original database (e.g. pfam|PB129746)

=cut

sub ext_in_family {
    my($self,$id) = @_;
    return () unless ($id);

    my($relational_db_response);
    my $rdbH = $self->db_handle;

    if ($rdbH->table_exists('localid_map') &&
        ($relational_db_response = $rdbH->SQL("SELECT DISTINCT family from localid_map WHERE localid = '$id'")))
    {
     my %seen; # only return the first occurence of anyting.
     return grep {!$seen{$_}++} map { $_->[0] } @$relational_db_response;
    }
    return ();
}

=head3 families_by_source

use: my @famlies = $fig->families_by_source('fig');

This use SQL to look up all the families that have a partial match to the argument supplied. It should be quicker than getting all families and parsing out the ones you want since it is done at the db level.

=cut

sub families_by_source {
 my ($self, $source)=@_;
 return () unless ($source);
 my($relational_db_response);
 my $rdbH = $self->db_handle;
 $source=lc($source);
 if (($relational_db_response= $rdbH->SQL("SELECT family from localfam_function WHERE family LIKE '$source\%'")) && $relational_db_response &&
    (@$relational_db_response >= 1))
    {
           return map { $_->[0] } @$relational_db_response;
    }
 else
    {
        return ();
    }
}

=head3 number_of_cids

use: my $number=$fig->number_of_cids

The number_of_ methods here all use SQL queries to count how many of each thing there are. This method just returns the number of cids

=cut

sub number_of_cids {
    my ($self)=@_;
    my($relational_db_response);
    my $rdbH = $self->db_handle;
    my $query="SELECT count(*) from (SELECT DISTINCT cid from localid_cid) as d";
    if (($relational_db_response= $rdbH->SQL($query)) && $relational_db_response) {return $relational_db_response->[0]->[0]}
    else {return undef}
}


=head3 number_of_families

use: my $number=$fig->number_of_families("fig");

This uses an SQL count method to count the number of families that match the given source. This should be a lot quicker than retrieving all families and then looping through them.

=cut

sub number_of_families {
    my ($self, $source)=@_;
    my($relational_db_response);
    my $rdbH = $self->db_handle;
    $source=lc($source);
    my $where="";
    $source && ($where .= " WHERE family LIKE '$source\%'");
    my $query="SELECT count(family) from (SELECT DISTINCT family from localfam_cid $where) as d";
    if (($relational_db_response= $rdbH->SQL($query)) && $relational_db_response) {return $relational_db_response->[0]->[0]}
    else {return undef}
}

=head3 number_of_proteins_in_families

use: my $number=$fig->number_of_proteins_in_families("fig", "distinct");

This uses and SQL count to count the number of proteins in families that match a given source. If distinct is true each protein will only be counted once, else the total number will be returned.

=cut

sub number_of_proteins_in_families {
    my ($self, $source, $distinct)=@_;
    my($relational_db_response);
    my $rdbH = $self->db_handle;
    $source=lc($source);
    my $query="SELECT count(localid) from ";
    my $where="";
    $source && ($where = "where localid like '$source\%'"); # only construct the where clause if we have a source, otherwise, we'll count everything
    $distinct ? ($query.="(SELECT DISTINCT localid from localid_map $where) as d") : ($query.="localid_map $where");
    if (($relational_db_response= $rdbH->SQL($query)) && $relational_db_response) {return $relational_db_response->[0]->[0]}
    else {return undef}
}


=head3 prot_to_cid

Convert a protein to a global ID
my $cid=$fig->prot_to_cid($proteinid)

$proteinid can be a FIG ID, a SP, tigr, or one of many other IDs

returns "" if not known

=cut


sub prot_to_cid {
    my($self,$prot) = @_;
    return "" unless ($prot);

    my($relational_db_response);
    my $rdbH = $self->db_handle;

    if ($rdbH->table_exists('localid_cid') &&
        ($relational_db_response = $rdbH->SQL("SELECT DISTINCT cid from localid_cid WHERE localid = '$prot'")) &&
        (@$relational_db_response == 1))
    {
        return $relational_db_response->[0]->[0];
    }
    return "";
}

=head3 cid_to_prots

Convert an internal ID to the proteins that map to that ID.
my @proteins=$fig->cid_to_prots($cid);

=cut

sub cid_to_prots {
    my($self,$cid) = @_;
    return () unless ($cid);

    my($relational_db_response);
    my $rdbH = $self->db_handle;

    if ($rdbH->table_exists('localid_cid') &&
        ($relational_db_response = $rdbH->SQL("SELECT DISTINCT localid from localid_cid WHERE cid = $cid")) &&
        (@$relational_db_response >= 1))
    {
        return map { $_->[0] } @$relational_db_response;
    }
    return ();
}


=head3 family_by_function

Get a list of families that have a partial match to a provided function.

E.g. my @families=$fig->family_by_function("histidine")

will return histidine kinase, histidine phosphatase, etc etc etc

=cut

sub family_by_function {
 my ($self, $func)=@_;
 return () unless ($func);
 my($relational_db_response);
 my $rdbH = $self->db_handle;
 $func=lc($func);

 if ($rdbH->table_exists('localfam_function') &&
     ($relational_db_response = $rdbH->SQL("SELECT DISTINCT family from localfam_function where lower(function) like '\%$func\%'")) &&
     (@$relational_db_response >= 1))
 {
     return map { $_->[0] } @$relational_db_response;
 }
 return ();
}

#
# PATtyfam lookup functions.
#

sub all_pattyfam_functions
{
    my($self) = @_;
    my $out;
    
    eval {
    my $res = $self->db_handle->SQL("SELECT DISTINCT(family_function) FROM family_membership");
    $out = $res;
    };
    return $out;
}

# "WHERE family = ?" for the WHERE, and instead of @$ids, just $famID.

sub pattyfam_members
{
    my($self, $fam_id) = @_;
    my @out;
    
    eval {
	my $res = $self->db_handle->SQL("SELECT fid, family, family_function FROM family_membership WHERE family = ?",
					undef, $fam_id);
	
	for my $ent (@$res)
	{
	    my($fid, $family, $fun) = @$ent;
	    
	    push @out, [$fid, $fun];
	}
    };
    return \@out;
}


sub pattyfams_for_proteins
{
    my($self, $ids) = @_;

    return if @$ids == 0;

    my %out;

    eval {
	my $qs = join(", ", map { "?" } @$ids);
	
	my $res = $self->db_handle->SQL("SELECT fid, family, family_function FROM family_membership WHERE fid IN ($qs)",
					undef, @$ids);
	
	for my $ent (@$res)
	{
	    my($fid, $family, $fun) = @$ent;
	    
	    push(@{$out{$fid}}, [$family, $fun]);
	}
    };
    return \%out;
}

sub pattyfams_for_function
{
    my($self, $fn) = @_;

    my $res;
    eval {

	my $op = $fn =~ /%/ ? 'LIKE' : '=';
	$res = $self->db_handle->SQL("SELECT fid, family, family_function FROM family_membership WHERE function $op ?",
					undef, $fn);
    };
    return $res;
}

################################# Abstract Set Routines  ####################################

=head2 Abstract Set Routines

=cut

sub all_sets {
    my($self,$relation,$set_name) = @_;
    my($relational_db_response);

    my $rdbH = $self->db_handle;

    if (($relational_db_response = $rdbH->SQL("SELECT DISTINCT $set_name FROM $relation")) &&
        (@$relational_db_response >= 1))
    {
        return map { $_->[0] } @$relational_db_response;
    }
    return ();
}

sub next_set {
    my($self,$relation,$set_name) = @_;
    my($relational_db_response);

    my $rdbH = $self->db_handle;

    if (($relational_db_response = $rdbH->SQL("SELECT MAX($set_name) FROM $relation")) &&
        (@$relational_db_response == 1))
    {
        return $relational_db_response->[0]->[0] + 1;
    }
}

sub ids_in_set {
    my($self,$which,$relation,$set_name) = @_;
    my($relational_db_response);

    my $rdbH = $self->db_handle;
    if (defined($which) && ($which =~ /^\d+$/))
    {
        if (($relational_db_response = $rdbH->SQL("SELECT id FROM $relation WHERE ( $set_name = $which)")) &&
            (@$relational_db_response >= 1))
        {
	    my @fids =  map { $_->[0] } @$relational_db_response;
	    my $del = $self->is_deleted_fid_bulk(@fids);
            return grep { ! $del->{$_} }
	    sort { by_fig_id($a,$b) }
	    @fids;
        }
    }
    return ();
}

sub in_sets {
    my($self,$id,$relation,$set_name) = @_;
    my($relational_db_response);

    if ($self->is_deleted_fid($id)) { return () }

    my $rdbH = $self->db_handle;

    if (($relational_db_response = $rdbH->SQL("SELECT $set_name FROM $relation WHERE ( id = \'$id\' )")) &&
        (@$relational_db_response >= 1))
    {
        return map { $_->[0] } @$relational_db_response;
    }
    return ();
}

sub sz_set {
    my($self,$which,$relation,$set_name) = @_;
    my($relational_db_response);

    my $rdbH = $self->db_handle;

    if (($relational_db_response = $rdbH->SQL("SELECT COUNT(*) FROM $relation WHERE ( $set_name = $which)")) &&
        (@$relational_db_response == 1))
    {
        return $relational_db_response->[0]->[0];
    }
    return 0;
}

sub delete_set {
    my($self,$set,$relation,$set_name) = @_;

#   print STDERR "deleting set $set\n";
    my $rdbH = $self->db_handle;

    return $rdbH->SQL("DELETE FROM $relation WHERE ( $set_name = $set )");
}

sub insert_set {
    my($self,$set,$ids,$relation,$set_name) = @_;
    my($id);

#   print STDERR "inserting set $set containing ",join(",",@$ids),"\n";
    my $rdbH = $self->db_handle;

    my @ids = grep { length($_) < 255 } @$ids;
    if (@ids < 2) { return 0 }

    my $rc = 1;
    foreach $id (@ids)
    {
        next if ($self->is_deleted_fid($id));
        if (! $rdbH->SQL("INSERT INTO $relation ( $set_name,id ) VALUES ( $set,\'$id\' )"))
        {
            $rc = 0;
        }
    }
#   print STDERR "    rc=$rc\n";
    return $rc;
}

sub in_set_with {
    my($self,$peg,$relation,$set_name) = @_;
    my($set,$id,%in);

    foreach $set ($self->in_sets($peg,$relation,$set_name))
    {
        foreach $id ($self->ids_in_set($set,$relation,$set_name))
        {
            $in{$id} = 1;
        }
    }
    return sort { &by_fig_id($a,$b) } keys(%in);
}


sub export_set {
    my($self,$relation,$set_name,$file) = @_;
    my($pair);

    my $rdbH = $self->db_handle;
    my $relational_db_response = $rdbH->SQL("SELECT $set_name, id FROM $relation");

    open(TMPSET,">$file")
        || die "could not open $file";
    flock(TMPSET,LOCK_EX) || confess "cannot lock $file";
    seek(TMPSET,0,2)      || confess "failed to seek to the end of the file";

    foreach $pair (sort { ($a->[0] <=> $b->[0]) or &by_fig_id($a->[1],$b->[1]) } @$relational_db_response)
    {
        if (! $self->is_deleted_fid($pair->[1]))
        {
            print TMPSET join("\t",@$pair),"\n";
        }
    }
    close(TMPSET);
    return 1;
}

################################# KEGG Stuff  ####################################

=head2 KEGG methods

=head3 all_compounds

    my @compounds = $fig->all_compounds();

Return a list containing all of the KEGG compounds.

=cut

sub all_compounds {
    my($self) = @_;

    my $rdbH = $self->db_handle;
    my $relational_db_response = $rdbH->SQL("SELECT DISTINCT cid FROM comp_name");
    if (@$relational_db_response > 0)
    {
        return sort map { $_->[0] } @$relational_db_response;
    }
    return ();
}

=head3 names_of_compound

    my @names = $fig->names_of_compound($cid);

Returns a list containing all of the names assigned to the specified KEGG compound. The list
will be ordered as given by KEGG.

=over 4

=item cid

ID of the desired compound.

=item RETURN

Returns a list of names for the specified compound.

=back

=cut

sub names_of_compound {
    my($self,$cid) = @_;

    my $rdbH = $self->db_handle;
    my $relational_db_response = $rdbH->SQL("SELECT pos,name FROM comp_name where cid = \'$cid\'");
    if (@$relational_db_response > 0)
    {
        return map { $_->[1] } sort { $a->[0] <=> $b->[0] } @$relational_db_response;
    }
    return ();
}

=head3 ids_of_compound

usage: @ids = $fig->ids_of_compound

Returns a list containing all of the ids assigned to the KEGG compounds.  The list
will be ordered as given by KEGG.

=cut

sub ids_of_compound {
    my($self,$name) = @_;

    my $rdbH = $self->db_handle;
    my $relational_db_response = $rdbH->SQL("SELECT pos,cid FROM comp_name where name = \'$name\'");
    if (@$relational_db_response > 0)
    {
        return map { $_->[1] } sort { $a->[0] <=> $b->[0] } @$relational_db_response;
    }
    return ();
}

=head3 ids_of_compound_like_name

usage: @ids = $fig->ids_of_compound_like_name($name)

Returns a list containing all of the ids assigned to the KEGG compounds that match $name.  The list
will be ordered as given by KEGG.

=cut

sub ids_of_compound_like_name {
    my($self,$name) = @_;

    # replace dashes with underscores, which will match any single character in the 'like' clause
    $name =~ s/-/_/g;
    my $rdbH = $self->db_handle;
    my $relational_db_response = $rdbH->SQL("SELECT pos,cid FROM comp_name where name ilike \'$name\'");
    if (@$relational_db_response > 0)
    {
        return map { $_->[1] } sort { $a->[0] <=> $b->[0] } @$relational_db_response;
    }
    return ();
}



=head3 comp2react

    my @rids = $fig->comp2react($cid);

Returns a list containing all of the reaction IDs for reactions that take $cid
as either a substrate or a product.

=cut

sub comp2react {
    my($self,$cid) = @_;

    my $rdbH = $self->db_handle;
    my $relational_db_response = $rdbH->SQL("SELECT rid FROM reaction_to_compound where cid = \'$cid\'");
    if (@$relational_db_response > 0)
    {
        return sort map { $_->[0] } @$relational_db_response;
    }
    return ();
}

=head3 valid_reaction_id

    my $flag = $fig->valid_reaction_id($rid);

Returns true iff the specified ID is a valid reaction ID.

This will become important as we include non-KEGG reactions

=over 4

=item rid

Reaction ID to test.

=item RETURN

Returns TRUE if the reaction ID is in the data store, else FALSE.

=back

=cut

sub valid_reaction_id
{
    my($self,$rid) = @_;

    my $rdbH = $self->db_handle;
    my $relational_db_response = $rdbH->SQL("SELECT rid FROM reaction_to_compound WHERE rid = '$rid'");
    return (@$relational_db_response > 0);
}

=head3 cas

    my $cas = $fig->cas($cid);

Return the Chemical Abstract Service (CAS) ID for the compound, if known.

=over 4

=item cid

ID of the compound whose CAS ID is desired.

=item RETURN

Returns the CAS ID of the specified compound, or an empty string if the CAS ID
is not known or does not exist.

=back

=cut

sub cas {
    my($self,$cid) = @_;

    my $rdbH = $self->db_handle;
    my $relational_db_response = $rdbH->SQL("SELECT cas FROM comp_cas where cid = \'$cid\'");
    if (@$relational_db_response == 1)
    {
        return $relational_db_response->[0]->[0];
    }
    return "";
}

=head3 cas_to_cid

    my $cid = $fig->cas_to_cid($cas);

Return the compound id (cid), given the Chemical Abstract Service (CAS) ID.

=over 4

=item cas

CAS ID of the desired compound.

=item RETURN

Returns the ID of the compound corresponding to the specified CAS ID, or an empty
string if the CAS ID is not in the data store.

=back

=cut

sub cas_to_cid {
    my($self,$cas) = @_;

    my $rdbH = $self->db_handle;
    my $relational_db_response = $rdbH->SQL("SELECT cid FROM comp_cas where cas = \'$cas\'");
    if (@$relational_db_response == 1)
    {
        return $relational_db_response->[0]->[0];
    }
    return "";
}

=head3 all_reactions

    my @rids = $fig->all_reactions();

Return a list containing all of the KEGG reaction IDs.

=cut

sub all_reactions {
    my($self) = @_;

    my $rdbH = $self->db_handle;
    my $relational_db_response = $rdbH->SQL("SELECT DISTINCT rid FROM reaction_to_compound");
    if (@$relational_db_response > 0)
    {
        return sort map { $_->[0] } @$relational_db_response;
    }
    return ();
}

=head3 reversible

    my $flag = $fig->reversible($rid);

Return TRUE if the specified reaction is reversible. A reversible reaction has no main
direction. The connector is symbolized by C<< <=> >> instead of C<< => >>.

=over 4

=item rid

ID of the ralevant reaction.

=item RETURN

Returns TRUE if the specified reaction is reversible, else FALSE. If the reaction
does not exist, returns TRUE.

=back

=cut

sub reversible {
    my ($self, $rid) = @_;

    my $rdbH = $self->db_handle;
    my $relational_db_response = $rdbH->SQL("SELECT reversible FROM reversible where rid = \'$rid\'");
    if (@$relational_db_response == 1)
    {
        return $relational_db_response->[0]->[0];
    }
    return 1;
}

=head3 reaction_direction

    my $rev = $fig->reaction_direction($rid);

Returns an array of triplets mapping from reactions in the context of maps to reversibility.

=over 4

=item rid

ID of the relevant reaction.

=item RETURN

Return C<< B >> if the reaction proceeds in both directions, C<< L >> if it proceeds from right
to left, or C<< R >> if it proceeds from left to right (by convention the "substrates"
are on the left and the "products" are on the right).

=back

=cut

sub reaction_direction {
    my ($self, $rid) = @_;
    my @results = ();

    my $rdbH = $self->db_handle;
    my $relational_db_response = $rdbH->SQL("SELECT rid, mapid, direction FROM reaction_direction where rid = \'$rid\'");

    if (@$relational_db_response > 0)
    {
	foreach my $res (@$relational_db_response) {
	    my ($rid, $mapid, $rev)=@$res;
	    push (@results, [$rid, $mapid, $rev]);
	}
    }

    return @results;
}

=head3 reaction2comp

    my @tuples = $fig->reaction2comp($rid, $which, $paths);

Return the substrates or products for a reaction.  In any event (i.e.,
whether you ask for substrates or products), you get back a list of
3-tuples.  Each 3-tuple will contain

    [$cid,$stoich,$main]

Stoichiometry indicates how many copies of the compound participate in
the reaction. It is normally numeric, but can be things like "n" or "(n+1)".
$main is 1 iff the compound is considered "main" or "connectable".

=over 4

=item rid

ID of the reaction whose compounds are desired.

=item which

TRUE if the products (right side) should be returned, FALSE if the substrates
(left side) should be returned.

=item paths

Optional list of paths to check whether compound is "main"

=item RETURN

Returns a list of 3-tuples. Each tuple contains the ID of a compound, its
stoichiometry, and a flag that is TRUE if the compound is one of the main
participants in the reaction.  If paths are specified, the flag indicates
whether the compound is main in any of the specified paths.

=back

=cut

sub reaction2comp {
    my($self,$rid,$which,$paths) = @_;

    my $rdbH = $self->db_handle;
    my $relational_db_response_not_main = $rdbH->SQL("SELECT cid,stoich,main FROM reaction_to_compound where rid = \'$rid\' and setn = \'$which\' and main = \'0\'");
    my $relational_db_response_main = $rdbH->SQL("SELECT distinct cid,stoich,main FROM reaction_to_compound where rid = \'$rid\' and setn = \'$which\' and main = \'1\'");

    if (@$relational_db_response_not_main > 0 || @$relational_db_response_main > 0)
    {
	my @tuples_to_return = @$relational_db_response_not_main;

	if (! $paths || scalar @$paths == 0)
	{
	    push @tuples_to_return, @$relational_db_response_main;
	}
	else
	{
	    my $inner_paths_string = join "','", @$paths;

	    foreach my $tuple (@$relational_db_response_main)
	    {
		my $relational_db_response_main_path = $rdbH->SQL("SELECT cid,stoich,main FROM reaction_to_compound where rid = \'$rid\' and setn = \'$which\' and main = \'1\' and cid = \'$tuple->[0]\' and path in \(\'$inner_paths_string\'\)");

		push @tuples_to_return, [$tuple->[0], $tuple->[1], @$relational_db_response_main_path > 0 ? "1" : "0"];
	    }
	}

	return sort { $a->[0] cmp $b->[0] } map { $_->[1] =~ s/\s+//g; $_ }  @tuples_to_return;
    }

    return ();
}

=head3 catalyzed_by

    my @ecs = $fig->catalyzed_by($rid);

Return the ECs (roles) that are reputed to catalyze the reaction.  Note that we are currently
just returning the ECs that KEGG gives.  We need to handle the incompletely specified forms
(e.g., 1.1.1.-), but we do not do it yet.

=over 4

=item rid

ID of the reaction whose catalyzing roles are desired.

=item RETURN

Returns the IDs of the roles that catalyze the reaction.

=back

=cut

sub catalyzed_by {
    my($self,$rid) = @_;

    my $rdbH = $self->db_handle;
    my $relational_db_response = $rdbH->SQL("SELECT role FROM reaction_to_enzyme where rid = \'$rid\'");
    if (@$relational_db_response > 0)
    {
        return sort map { $_->[0] } @$relational_db_response;
    }
    return ();
}

=head3 catalyzes

    my @ecs = $fig->catalyzes($role);

Returns the reaction IDs of the reactions catalyzed by the specified role (normally an EC).

=over 4

=item role

ID of the role whose reactions are desired.

=item RETURN

Returns a list containing the IDs of the reactions catalyzed by the role.

=back

=cut

sub catalyzes {
    my ($self, $role) = @_;

    my $rdbH = $self->db_handle;
    $role = quotemeta $role;
    my $relational_db_response = $rdbH->SQL("SELECT rid FROM reaction_to_enzyme where role = \'$role\'");
    if (@$relational_db_response > 0)
    {
        return sort map { $_->[0] } @$relational_db_response;
    }
    return ();
}


=head3 displayable_reaction

    my $displayString = $fig->displayable_reaction($rid)

Returns a string giving the displayable version of a reaction.

=cut

sub displayable_reaction {
    my($self,$rid) = @_;

    my @tmp = `grep $rid $FIG_Config::data/KEGG/reaction_name.lst`;
    if (@tmp > 0)
    {
        chomp $tmp[0];
        return $tmp[0];
    }
    return $rid;
}

=head3 all_maps

    my @maps = $fig->all_maps();

Return all of the KEGG maps in the data store.

=cut

sub all_maps {
    my($self) = @_;

    my $cached = $self->cached("_ec_map_cache");
    my $pairs = $cached->{pairs};
    if (!ref($pairs))
    {
	my $rdbH = $self->db_handle;
	$pairs = $rdbH->SQL("SELECT ec, map FROM ec_map ");
	$cached->{pairs} = $pairs;
    }

    my %maps;
    $maps{$_->[1]} = 1 for @$pairs;
    return sort keys %maps;
}

sub roles_for_prot {
    my($self, $prot) = @_;

    $prot = quotemeta $prot;

    my $rdbH = $self->db_handle;
    my $relational_db_response = $rdbH->SQL("SELECT role FROM roles WHERE prot='$prot' ");
    if (@$relational_db_response > 0)
    {
        return map { $_->[0] =~ s/\s+$//; $_->[0] } @$relational_db_response;
    }
    return ();
}

sub prots_for_role {
    my($self, $role) = @_;

    $role = quotemeta $role;

    my $rdbH = $self->db_handle;
    my $relational_db_response = $rdbH->SQL("SELECT DISTINCT prot FROM roles WHERE role='$role' AND prot LIKE 'fig|%' AND NOT prot LIKE 'fig|9999999%' ");
    if (@$relational_db_response > 0)
    {
        return map { $_->[0] } @$relational_db_response;
    }
    return ();
}

=head3 ec_to_maps

    my @maps = $fig->ec_to_maps($ec);

Return the set of maps that contain a specific functional role. The role can be
specified by an EC number or a full-blown role ID.

=over 4

=item ec

The EC number or role ID of the role whose maps are desired.

=item RETURN

Returns a list of the IDs for the maps that contain the specified role.

=back

=cut

sub ec_to_maps {
    my($self,$ec) = @_;

    my $rdbH = $self->db_handle;
    my $relational_db_response = $rdbH->SQL("SELECT map FROM ec_map WHERE ( ec = \'$ec\' )");
    if (@$relational_db_response > 0)
    {
        return map { $_->[0] } @$relational_db_response;
    }
    return ();
}

=head3 role_to_maps

This is an alternate name for L</ec_to_maps>.

=cut

sub role_to_maps {
    my ($self, $role) = @_;
    return $self->ec_to_maps($role);
}

=head3 map_to_ecs

    my @ecs = $fig->map_to_ecs($map);

Return the set of functional roles (usually ECs) that are contained in the functionality
depicted by a map.

=over 4

=item map

ID of the KEGG map whose roles are desired.

=item RETURN

Returns a list of EC numbers for the roles in the specified map.

=back

=cut

sub map_to_ecs {
    my($self,$map) = @_;

    my $rdbH = $self->db_handle;
    my $relational_db_response = $rdbH->SQL("SELECT ec FROM ec_map WHERE ( map = \'$map\' )");
    if (@$relational_db_response > 0)
    {
        return map { $_->[0] } @$relational_db_response;
    }
    return ();
}

=head3 map_name

    my $name = $fig->map_name($map);

Return the descriptive name covering the functionality depicted by the specified map.

=over 4

=item map

ID of the map whose description is desired.

=item RETURN

Returns the descriptive name of the map, or an empty string if no description is available.

=back

=cut

sub map_name {
    my($self,$map) = @_;

    my $cache = $self->cached("_map_name");
    if (!%$cache)
    {
	my $res = $self->db_handle->dbh->selectall_hashref(qq(SELECT map, mapname FROM map_name), 'map');
	%$cache = %$res;
    }
    return $cache->{$map}->{mapname};
}

################################# Functional Roles  ####################################
=head2 Functional Roles

=head3 neighborhood_of_role

usage: @roles = $fig->neighborhood_of_role($role)

Returns a list of functional roles that we consider to be "the neighborhood" of $role.

=cut

sub neighborhood_of_role {
    my($self,$role) = @_;
    my($readC);

    my $file = "$FIG_Config::global/role.neighborhoods";
    my $rdbH = $self->db_handle;
    my $roleQ = quotemeta $role;
    my $relational_db_response = $rdbH->SQL("SELECT seek, len FROM neigh_seeks WHERE role = \'$roleQ\' ");
    if (@$relational_db_response == 1)
    {
        my($seek,$ln) = @{$relational_db_response->[0]};
        my $fh   = $self->openF($file);
        seek($fh,$seek,0);
        my $readN = read($fh,$readC,$ln-1);
        ($readN == ($ln-1))
            || confess "could not read the block of sims at $seek for $ln - 1 characters; $readN actually read from $file\n$readC";
        return grep { $_ && ($_ !~ /^\/\//) } split(/\n/,$readC);
    }
    return ();
}

=head3 roles_of_function

    my @roles = $fig->roles_of_function($func);

Returns a list of the functional roles implemented by the specified function. This method
parses the role data out of the function name, and does not require access to the database.

=over 4

=item func

Name of the function whose roles are to be parsed out.

=item RETURN

Returns a list of the roles performed by the specified function.

=back

=cut

sub roles_of_function {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my $func = (@_ == 1) ? $_[0] : $_[1];

    $func =~ s/\s*[\!\#].*$//;
    my %roles = map { $_ => 1 } (split(/\s*;\s+|\s+[\@\/]\s+/,$func),($func =~ /\d+\.\d+\.\d+\.\d+/g),$func);
    return sort keys(%roles);
}

sub function_to_subsystems {
    my($self,$func) = @_;

    my %subs;
    my @roles = $self->roles_of_function($func);
    if (@roles > 0)
    {
        foreach my $role (@roles)
        {
            foreach my $sub ($self->role_to_subsystems($role))
            {
                $subs{$sub} = 1;
            }
        }
    }
    return sort keys(%subs);
}

=head3 protein_subsystem_to_roles

    my $roles = $fig->protein_subsystem_to_roles($peg, $subsystem);

Return the roles played by a particular PEG in a particular subsytem. If the protein is not part of the
subsystem, an empty list will be returned.

=over 4

=item peg

ID of the protein whose role is desired.

=item subsystem

Name of the relevant subsystem.

=item RETURN

Returns a reference to a list of the roles performed by the specified PEG in the specified subsystem.

=back

=cut

sub protein_subsystem_to_roles {
    my($self,$prot, $subsystem) = @_;

    my($relational_db_response);
    my $rdbH = $self->db_handle;
    my $query = "SELECT role FROM subsystem_index WHERE protein=? AND subsystem=?";
    $relational_db_response = $rdbH->SQL($query, undef, $prot, $subsystem);
    my $result = [];
    if (scalar(@$relational_db_response)) {
	@$result = map { $_->[0] } @$relational_db_response;
    }
    return $result;
}

sub role_to_subsystems {
    my($self,$role) = @_;

    my($relational_db_response);
    my $rdbH = $self->db_handle;
    my $query;
    if ($FIG_Config::exclude_experimental_subsystems)
    {
	$query = qq(SELECT distinct subsystem
		    FROM subsystem_index i JOIN subsystem_metadata m ON i.subsystem = m.subsystem
		    WHERE  role = ? AND m.class_1 NOT LIKE 'Experimental%');
    }
    else
    {
	$query = "SELECT distinct subsystem FROM subsystem_index  WHERE  role = ?";
    }
    return (($relational_db_response = $rdbH->SQL($query, undef, $role)) && (@$relational_db_response >= 1)) ?
        map { $_->[0] } @$relational_db_response : ();
}

=head3 is_BRC_genome

$fig->is_BRC_genome($genome)
returns true if $genome is an BRC genome

=cut

sub is_BRC_genome {
    my($self,$org) = @_;

    return (-e "$FIG_Config::organisms/$org/BRC") ? 1 : 0;
}

=head3 is_NMPDR_genome

$fig->is_NMPDR_genome($genome)
returns true if $genome is an NMPDR genome

=cut

sub is_NMPDR_genome {
    my($self,$org) = @_;

    return (-e "$FIG_Config::organisms/$org/NMPDR") ? 1 : 0;
}

=head3 seqs_with_role

    my @pegs = $fig->seqs_with_role($role,$who);

Return a list of the pegs that implement $role.  If $who is not given, it
defaults to "master".  The system returns all pegs with an assignment made by
either "master" or $who (if it is different than the master) that implement $role.
Note that this includes pegs for which the "master" annotation disagrees with that
of $who, the master's implements $role, and $who's does not.

=cut

sub seqs_with_role {
    my($self,$role,$who,$genome) = @_;
    my($relational_db_response,$query);

    my $roleQ = quotemeta $role;

    $who = $who ? $who : "master";
    my $rdbH = $self->db_handle;

    my $who_cond;
    if ($who eq "master")
    {
        $who_cond = "( made_by = \'master\' OR made_by = \'unknown\' )";
    }
    else
    {
        $who_cond = "( made_by = \'master\' OR made_by = \'$who\' OR made_by = \'unknown\')";
    }

    if (! $genome)
    {
        $query = "SELECT distinct prot FROM roles  WHERE (( role = \'$roleQ\' ) AND $who_cond )";
    }
    else
    {
        $query = "SELECT distinct prot FROM roles  WHERE (( role = \'$roleQ\' ) AND $who_cond AND (org = \'$genome\'))";
    }
    return (($relational_db_response = $rdbH->SQL($query)) && (@$relational_db_response >= 1)) ?
        grep { ! $self->is_deleted_fid($_) } map { $_->[0] } @$relational_db_response : ();
}

=head3 seqs_with_roles_in_genomes

usage: $result = $fig->seqs_with_roles_in_genomes($genomes,$roles,$made_by)

This routine takes a pointer to a list of genomes ($genomes) and a pointer to a list of
roles ($roles) and looks up all of the sequences that connect to those roles according
to either the master assignments or those made by $made_by.  Again, you will get assignments
for which the "master" assignment connects, but the $made_by does not.

A hash is returned.  The keys to the hash are genome IDs for which at least one sequence
was found.  $result->{$genome} will itself be a hash, assuming that at least one sequence
was found for $genome.  $result->{$genome}->{$role} will be set to a pointer to a list of
2-tuples.  Each 2-tuple will contain [$peg,$function], where $function is the one for
$made_by (which may not be the one that connected).

=cut

sub seqs_with_roles_in_genomes {
    my($self,$genomes,$roles,$made_by) = @_;
    my($genome,$role,$roleQ,$role_cond,$made_by_cond,$query,$relational_db_response,$peg,$genome_cond,$hit);
    my $rdbH = $self->db_handle;
    my $result = {}; # foreach $genome ($self->genomes) { $result->{$genome} = {} }
    if (! $made_by) { $made_by = 'master' }
    if ((@$genomes > 0) && (@$roles > 0))
    {
        $genome_cond = "(" . join(" OR ",map { "( org = '$_' )" } @$genomes) . ")";
        $role_cond   = "(" . join(" OR ",map { $roleQ = quotemeta $_; "( role = '$roleQ' )" } @$roles) . ")";
        $made_by_cond = ($made_by eq 'master') ? "(made_by = 'master')" : "(made_by = 'master' OR made_by = '$made_by')";
        $query = "SELECT distinct prot, role FROM roles  WHERE ( $genome_cond AND $role_cond )";
        if (($relational_db_response = $rdbH->SQL($query)) && (@$relational_db_response >= 1))
        {
            foreach $hit (@$relational_db_response)
            {
                ($peg,$role) = @$hit;
                if (! $self->is_deleted_fid($peg))
                {
                    $genome = $self->genome_of($peg);
                    push(@{ $result->{$genome}->{$role} },[$peg,scalar $self->function_of($peg,$made_by)]);
                }
            }
        }
    }
    return $result;
}

=head3 largest_clusters

usage: @clusters = $fig->largest_clusters($roles,$user)

This routine can be used to find the largest clusters containing some of the
designated set of roles.  A list of clusters is returned.  Each cluster is a pointer to
a list of pegs.

=cut

sub largest_clusters {
    my($self,$roles,$user,$sort_by_unique_functions) = @_;
    my($genome,$x,$role,$y,$peg,$loc,$contig,$beg,$end,%pegs,@pegs,$i,$j);

    my $ss = $self->seqs_with_roles_in_genomes([$self->genomes],$roles,$user);

    my @clusters = ();

    foreach $genome (keys(%$ss))
    {
        my %pegs;
        $x = $ss->{$genome};
        foreach $role (keys(%$x))
        {
            $y = $x->{$role};
            foreach $peg (map { $_->[0] } @$y)
            {
                if ($loc = $self->feature_location($peg))
                {
                    ($contig,$beg,$end) = $self->boundaries_of($loc);
                    $pegs{$peg} = [$peg,$contig,int(($beg + $end) / 2)];
                }
            }
        }

        @pegs = sort { ($pegs{$a}->[1] cmp $pegs{$b}->[1]) or ($pegs{$a}->[2] <=> $pegs{$b}->[2]) } keys(%pegs);
        $i = 0;
        while ($i < $#pegs)
        {
            for ($j=$i+1; ($j < @pegs) && &close_enough_locs($pegs{$pegs[$j-1]},$pegs{$pegs[$j]}); $j++) {}
            if ($j > ($i+1))
            {
                push(@clusters,[@pegs[$i..$j-1]]);
            }
            $i = $j;
        }
    }
    if ($sort_by_unique_functions)
    {
        @clusters = sort { $self->unique_functions($b,$user) <=> $self->unique_functions($a,$user) } @clusters;
    }
    else
    {
        @clusters = sort { @$b <=> @$a } @clusters;
    }
    return @clusters;
}

sub unique_functions {
    my($self,$pegs,$user) = @_;
    my($peg,$func,%seen);

    foreach $peg (@$pegs)
    {
        if ($func = $self->function_of($peg,$user))
        {
            $seen{$func} = 1;
        }
    }
    return scalar keys(%seen);
}

sub close_enough_locs {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($x,$y) = @_;

    return (($x->[1] eq $y->[1]) && (abs($x->[2] - $y->[2]) < 5000));
}

sub candidates_for_role {
    my($self,$role,$genome,$cutoff,$user) = @_;
    my($peg);

    $user = $user ? $user : "master";

    my @cand = map { $_->[0] }
               sort { $a->[1] <=> $b->[1] }
               map { $peg = $_; [$peg,$self->crude_estimate_of_distance($genome,&FIG::genome_of($peg))] }
               $self->seqs_with_role($role,$user);

    return $self->candidates_for_role_from_known($genome,$cutoff,\@cand);
}

sub candidates_for_role_from_known {
    my($self,$genome,$cutoff,$known) = @_;
    my($peg);

    my @cand = @$known;
    my $hits = {};
    my $seen = {};
    my $how_many = (@cand > 10) ? 9 : $#cand;
    &try_to_locate($self,$genome,$hits,[@cand[0..$how_many]],$seen,$cutoff);
    if (keys(%$hits) == 0)
    {
        splice(@cand,0,$how_many+1);
        &try_to_locate($self,$genome,$hits,\@cand,$seen,$cutoff);
    }
    return sort {$hits->{$a}->[0] <=> $hits->{$b}->[0]} keys(%$hits);
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

=head2 Bidirectional Best Hits

=head3 best_bbh_candidates

usage: @candidates = $fig->best_bbh_candidates($genome,$cutoff,$requested,$known)

This routine returns a list of up to $requested candidates from $genome.  A candidate is a BBH
against one of the PEGs in genomes from the list given by@$known.
Each entry in the list is a 3-tuple:

    [CandidatePEG,KnownBBH,Pscore]

=cut

sub best_bbh_candidates {
    my($self,$genome,$cutoff,$requested,$known,$frac_match) = @_;
    my($i,$j,$k,$sim,@sims,$peg,$id2,$genome2,$sim_back);
    my($bbh,%seen,%computed_sims,$genome1);

    $frac_match = defined($frac_match) ? $frac_match : 0.7;
    my @got = ();
    my @cand = $self->candidates_for_role_from_known($genome,$cutoff,$known);
    if (@cand > 0)
    {
        my %genomes = map { $genome1 = &FIG::genome_of($_); $genome1 => 1 } @$known;
        my %pegs    = map { $_ => 1 } @$known;
        for ($i=0; (@got < $requested) && ($i < @cand); $i++)
        {
            $peg = $cand[$i];
            undef %seen;
            @sims = grep { $genomes{&FIG::genome_of($_->id2)} } $self->sims($peg,1000,$cutoff,"fig");
            $bbh = 0;
            for ($j=0; (! $bbh) && ($j < @sims); $j++)
            {
                $sim = $sims[$j];
                $id2 = $sim->id2;
                $genome2 = &FIG::genome_of($id2);
                if (! $seen{$genome2})
                {
                    if ($pegs{$id2})
                    {
                        if (! defined($sim_back = $computed_sims{$id2}))
                        {
                            my @sims_back = $self->sims($id2,1000,$cutoff,"fig");
                            for ($k=0; ($k < @sims_back) && (&FIG::genome_of($sims_back[$k]->id2) ne $genome); $k++) {}
                            if ($k < @sims_back)
                            {
                                $sim_back = $computed_sims{$id2} = $sims_back[$k];
                            }
                            else
                            {
                                $sim_back = $computed_sims{$id2} = 0;
                            }
                        }
                        if ($sim_back)
                        {
                            if (($sim_back->id2 eq $peg) && $self->ok_match($sim_back,$frac_match))
                            {
                                $bbh = [$id2,$sim_back->psc];
                            }
                        }
                    }
                    $seen{$genome2} = 1;
                }
            }

            if ($bbh)
            {
                push(@got,[$peg,@$bbh]);
            }
        }
    }
    return @got;
}




=pod

=head3 best_bbh_candidates_additional

usage: @candidates = $fig->best_bbh_candidates_additional($genome,$cutoff,$requested,$known)

This routine returns a list of up to $requested candidates from $genome.  A candidate is a BBH
against one of the PEGs in genomes from the list given by@$known.
The method collects additional information from the similarities and is used in
the subsystem extension.
Each entry in the list is a 10-tuple:

    [CandidatePEG,KnownBBH,Pscore,fraction, b1, e1, b2, e2, ln1, ln2]

=cut

sub best_bbh_candidates_additional {
    my($self,$genome,$cutoff,$requested,$known,$frac_match) = @_;
    my($i,$j,$k,$sim,@sims,$peg,$id2,$genome2,$sim_back);
    my($bbh,%seen,%computed_sims,$genome1);

    $frac_match = defined($frac_match) ? $frac_match : 0.7;
    my @got = ();
    my @cand = $self->candidates_for_role_from_known($genome,$cutoff,$known);
    if (@cand > 0)
    {
        my %genomes = map { $genome1 = &FIG::genome_of($_); $genome1 => 1 } @$known;
        my %pegs    = map { $_ => 1 } @$known;
        for ($i=0; (@got < $requested) && ($i < @cand); $i++)
        {
            $peg = $cand[$i];
            undef %seen;
            @sims = grep { $genomes{&FIG::genome_of($_->id2)} } $self->sims($peg,1000,$cutoff,"fig");
            $bbh = 0;
            for ($j=0; (! $bbh) && ($j < @sims); $j++)
            {
                $sim = $sims[$j];
                $id2 = $sim->id2;
                $genome2 = &FIG::genome_of($id2);
                if (! $seen{$genome2})
                {
                    if ($pegs{$id2})
                    {
                        if (! defined($sim_back = $computed_sims{$id2}))
                        {
                            my @sims_back = $self->sims($id2,1000,$cutoff,"fig");
                            for ($k=0; ($k < @sims_back) && (&FIG::genome_of($sims_back[$k]->id2) ne $genome); $k++) {}
                            if ($k < @sims_back)
                            {
                                $sim_back = $computed_sims{$id2} = $sims_back[$k];
                            }
                            else
                            {
                                $sim_back = $computed_sims{$id2} = 0;
                            }
                        }
                        if ($sim_back)
                        {
                            if (($sim_back->id2 eq $peg) && $self->ok_match($sim_back,$frac_match))
                            {
                                my $frac = $self->min(($sim_back->e1+1 - $sim_back->b1) / $sim_back->ln1, ($sim_back->e2+1 - $sim_back->b2) / $sim_back->ln2);

                                $bbh = [$id2,$sim_back->psc,$frac,$sim_back->b1, $sim_back->e1, $sim_back->b2, $sim_back->e2, $sim_back->ln1, $sim_back->ln2 ];
                            }
                        }
                    }
                    $seen{$genome2} = 1;
                }
            }

            if ($bbh)
            {
                push(@got,[$peg,@$bbh]);
            }
        }
    }
    return @got;
}



sub ok_match {
    my($self,$sim,$frac_match) = @_;

    my $ln1 = $sim->ln1;
    my $ln2 = $sim->ln2;
    my $b1  = $sim->b1;
    my $e1  = $sim->e1;
    my $b2  = $sim->b2;
    my $e2  = $sim->e2;

    return (((($e1 - $b1) / $ln1) >= $frac_match) &&
            ((($e2 - $b2) / $ln2) >= $frac_match))
}

sub external_calls {
    my($self,$pegs) = @_;
    my($peg,$func);

    open(TMP,">/tmp/pegs.$$") || die "could not open /tmp/pegs.$$";
    foreach $peg (@$pegs)
    {
        print TMP "$peg\n";
    }
    close(TMP);
    open(TMP,">/tmp/parms.$$") || die "could not open /tmp/parms.$$";
    print TMP "no_fig\t1\n";
    close(TMP);

    my %call = map { chop; ($peg,$func) = split(/\t/,$_) }
                `$FIG_Config::bin/auto_assign /tmp/parms.$$ < /tmp/pegs.$$ 2> /dev/null | $FIG_Config::bin/make_calls`;
    unlink("/tmp/pegs.$$","/tmp/parms.$$");
    return map { $call{$_} ? [$_,$call{$_}] : [$_,"hypothetical protein"] } @$pegs;
}

use SameFunc;

sub same_func_why {
    my($self,$f1,$f2) = @_;

    return &SameFunc::same_func_why($f1,$f2);
}
sub same_func {
    my($self,$f1,$f2) = @_;

    return &SameFunc::same_func($f1,$f2);
}

#################################   DNA sequence Stuff ####################################

=head2 DNA Sequences

=head3 extract_seq

usage: $seq = &FIG::extract_seq($contigs,$loc)

This is just a little utility routine that I have found convenient.  It assumes that
$contigs is a hash that contains IDs as keys and sequences as values.  $loc must be of the
form

       Contig_Beg_End

where Contig is the ID of one of the sequences; Beg and End give the coordinates of the sought
subsequence.  If Beg > End, it is assumed that you want the reverse complement of the subsequence.
This routine plucks out the subsequence for you.

=cut

sub extract_seq {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($contigs,$loc) = @_;
    my($contig,$beg,$end,$contig_seq);
    my($plus,$minus);

    $plus = $minus = 0;
    my $strand = "";
    my @loc = split(/,/,$loc);
    my @seq = ();
    foreach $loc (@loc)
    {
        if ($loc =~ /^\S+_(\d+)_(\d+)$/)
        {
            if ($1 < $2)
            {
                $plus++;
            }
            elsif ($2 < $1)
            {
                $minus++;
            }
        }
    }
    if ($plus > $minus)
    {
        $strand = "+";
    }
    elsif ($plus < $minus)
    {
        $strand = "-";
    }

    foreach $loc (@loc)
    {
        if ($loc =~ /^(\S+)_(\d+)_(\d+)$/)
        {
            ($contig,$beg,$end) = ($1,$2,$3);

            my $len = length($contigs->{$contig});
            if (!$len)
            {
                carp "Undefined or zero-length contig $contig";
                return "";
            }

            if (($beg > $len) || ($end > $len))
            {
                carp "Region $loc out of bounds (contig len=$len)";
            }
            else
            {
                if (($beg < $end) || (($beg == $end) && ($strand eq "+")))
                {
                    push(@seq,substr($contigs->{$contig},$beg-1,($end+1-$beg)));
                }
                else
                {
                    $strand = "-";
                    push(@seq,&reverse_comp(substr($contigs->{$contig},$end-1,($beg+1-$end))));
                }
            }
        }
    }
    return join("",@seq);
}




=head3 contigs_of

    my @contig_ids = $fig->contigs_of($genome);

Returns a list of all of the contigs occurring in the designated genome.

=over 4

=item genome

ID of the genome whose contigs are desired.

=item RETURN

Returns a list of the IDs for the contigs occurring in the specified genome.

=back

=cut

sub contigs_of {
    my($self,$genome) = @_;
    return $self->all_contigs($genome);
}

=head3 number_of_contigs

usage: $n=$fig->number_of_contigs($genome)

This uses the SQL count function to count the numbmer of contigs. It should be a lot faster than pulling all the contigs and counting them.

In fact, it causes about a 10-fold increase in speed! Compare fig n_contigs and fig number_of_contigs

=cut

sub number_of_contigs {
    my ($self, $genome)=@_;
    my($rdbH,$relational_db_response);

    my $cached = $self->cached("_contig_counts");
    my $n = $cached->{$genome};
    return $n if defined($n);

    $rdbH = $self->db_handle;
    if (defined($genome))
    {
	my $counts = $rdbH->SQL(qq(SELECT genome, COUNT(contig) FROM contig_lengths GROUP BY genome));

	$cached->{$_->[0]} = $_->[1] for @$counts;
	return $cached->{$genome};
    }
    return undef;
}



=head3 all_contigs

usage: @contig_ids = $fig->all_contigs($genome)

Returns a list of all of the contigs occurring in the designated genome.

=cut
#: Return Type @;
sub all_contigs {
    my($self,$genome) = @_;
    my($rdbH,$relational_db_response);

    $rdbH = $self->db_handle;
    if (defined($genome))
    {
        if ($relational_db_response = $rdbH->SQL("SELECT DISTINCT contig FROM contig_lengths WHERE ( genome = \'$genome\' )"))
        {
            return map { $_->[0] } @$relational_db_response;
        }
    }
    return undef;
}

=head3 contig_ln

usage: $n = $fig->contig_ln($genome,$contig)

Returns the length of $contig from $genome.

=cut

sub contig_ln {
    my($self,$genome,$contig) = @_;
    my($rdbH,$relational_db_response);

    my $cached = $self->cached("_contig_lengths");

    my $ln = $cached->{$genome, $contig};

    return $ln if defined($ln);

    $rdbH = $self->db_handle;
    if (defined($genome) && defined($contig))
    {
        if (($relational_db_response = $rdbH->SQL("SELECT len FROM contig_lengths WHERE ( genome = \'$genome\' ) and ( contig = \'$contig\' )")) &&

            (@$relational_db_response == 1))
        {
            $ln = $relational_db_response->[0]->[0];
	    $cached->{$genome, $contig} = $ln;
	    return $ln;
        }
    }
    return undef;
}

sub contig_ln_bulk {
    my($self,$gc_list) = @_;
    my($rdbH,$relational_db_response);

    $rdbH = $self->db_handle;
    my $ret = {};
    if (ref($gc_list) eq 'ARRAY' && @$gc_list > 0)
    {
	my $where = join(" OR ",
			 map { "(genome = ? AND contig = ?)" } @$gc_list);

	my @vals = map { @$_ } @$gc_list;

        $relational_db_response = $rdbH->SQL(qq(SELECT genome, contig, len
						FROM contig_lengths
						WHERE $where), undef, @vals);
	for my $ent (@$relational_db_response)
	{
	    my($genome, $contig, $len) = @$ent;
	    $ret->{$genome}->{$contig} = $len;
        }
    }
    return $ret;
}


=head3 get_dna_seq

    my $seq = $fig->get_dna_seq($fid);

Returns the DNA sequence for an FID

=over 4

=item fid

FIG identifier of the feature whose sequence is desired

=item RETURN

DNA sequence

=back

=cut

sub get_dna_seq {
    my ($self, $fid) = @_;

    my $genome    = $self->genome_of( $fid );
    my @locations = $self->feature_location( $fid );

    my $seq = $self->dna_seq($genome, @locations);

    return $seq;
}


=head3 dna_seq

usage: $seq = $fig->dna_seq($genome,@locations)

Returns the concatenated subsequences described by the list of locations.  Each location
must be of the form

    Contig_Beg_End

where Contig must be the ID of a contig for genome $genome.  If Beg > End the location
describes a stretch of the complementary strand.

=cut

#: Return Type $;
sub dna_seq {
    my($self,$genome,@locations) = @_;
    my(@pieces,$loc,$contig,$beg,$end,$ln,$rdbH);

    @locations = map { split(/,/,$_) } @locations;
    @pieces = ();
    foreach $loc (@locations)
    {
        if ($loc =~ /^(\S+)_(\d+)_(\d+)$/)
        {
            ($contig,$beg,$end) = ($1,$2,$3);
            $ln = $self->contig_ln($genome,$contig);

            if (! $ln) {
                carp "dna_seq($genome, $loc): contig length undefined";
                return "";
            }

            if (&between(1,$beg,$ln) && &between(1,$end,$ln))
            {
                if ($beg < $end)
                {
                    push(@pieces, $self->get_dna($genome,$contig,$beg,$end));
                }
                else
                {
                    push(@pieces, &reverse_comp($self->get_dna($genome,$contig,$end,$beg)));
                }
            }
        }
    }
    return lc(join("",@pieces));
}

sub get_dna {
    my($self,$genome,$contig,$beg,$end) = @_;
    my $relational_db_response;

    my $rdbH = $self->db_handle;
    my $indexpt = int(($beg-1)/10000) * 10000;
    if (($relational_db_response = $rdbH->SQL("SELECT startN,fileno,seek FROM contig_seeks WHERE ( genome = \'$genome\' ) AND ( contig = \'$contig\' ) AND ( indexpt = $indexpt )")) &&
        (@$relational_db_response == 1))
    {
        my($startN,$fileN,$seek) = @{$relational_db_response->[0]};
        my $fh = $self->openF($self->N2file($fileN));
        if (seek($fh,$seek,0))
        {
            my $chunk = "";
            read($fh,$chunk,int(($end + 1 - $startN) * 1.03));
#           print STDERR "genome=$genome contig=$contig beg=$beg end=$end startN=$startN chunk=$chunk\n";
            $chunk =~ s/\s//g;
            my $ln = ($end - $beg) + 1;
            if (length($chunk) >= $ln)
            {
                return lc(substr($chunk,(($beg-1)-$startN),$ln));
            }
        }
    }
    return undef;
}

#################################   Taxonomy  ####################################

=head2 Taxonomy

=head3 taxonomy_of

usage: $taxonomy = $fig->taxonomy_of($genome_id)

Returns the taxonomy of the specified genome.  Gives the taxonomy down to
genus and species.

=cut

sub taxonomy_of :Scalar {
    my($self,$genome) = @_;
    my($ans);
    my $taxonomy = $self->cached('_taxonomy');

    $ans = $taxonomy->{$genome};

    if (!defined($ans)) {
	if (keys(%$taxonomy) == 0) {
	    my $rdbH = $self->db_handle;
	    my $relational_db_response = $rdbH->SQL("SELECT genome,taxonomy  FROM genome");
	    my $pair;
	    foreach $pair (@$relational_db_response) {
		my ($db_genome, $db_taxonomy) = @$pair;

		$db_taxonomy =~ s/^\s*//o;
		$db_taxonomy =~ s/Candidatus\s+//og;
		$db_taxonomy =~ s/\s+/ /og;
		$db_taxonomy =~ s/\s*$//o;

		$taxonomy->{$db_genome} = $db_taxonomy;
	    }
	    $ans = $taxonomy->{$genome};
	}
    }

    if (!$ans) {
	#warn "No taxonomy found for $genome\n";
    }
    else
    {
	$ans =~ s/^\s*//o;
	$ans =~ s/Candidatus\s*//og;
	$ans =~ s/\s+/ /og;
	$ans =~ s/\s*$//o;
    }

    return $ans;
}


=head3 get_taxonomy_id_of

usage: $taxonomyID = $fig->get_taxonomy_id_of($genome_id)

Returns the taxonomy ID of the specified genome. If no taxonomy ID is found the genome id without ".\d+" suffix will be returned.

=cut

sub get_taxonomy_id_of{
    my($self,$genome) = @_;

    my $tax_id = undef;
    if (-f "$FIG_Config::organisms/$genome/TAXONOMY_ID"){
	open (TAX , "$FIG_Config::organisms/$genome/TAXONOMY_ID") or die "Can't open $FIG_Config::organisms/$genome/TAXONOMY_ID\n";
	$tax_id = <TAX>;
	chomp $tax_id;
	close (TAX);
    }
    else{
	($tax_id) = $genome =~ /(\d+)\.\d+/;
    }

    return $tax_id;
}

=head3 set_taxonomy_id_for

usage: $taxonomyID = $fig->set_taxonomy_id_for($genome_id)

Sets the taxonomy id for genome.

=cut

sub set_taxonomy_id_for{
    my($self,$genome) = @_;

    my $tax_id = undef;
    if (-d "$FIG_Config::organisms/$genome/"){
	open (TAX , ">$FIG_Config::organisms/$genome/TAXONOMY_ID") or die "Can't open $FIG_Config::organisms/$genome/TAXONOMY_ID\n";
	print TAX "$tax_id\n" ;
	close (TAX);
    }
    else{
	print STDERR "No directory $FIG_Config::organisms/$genome/\n";
    }

    return $tax_id;
}



=head3 taxonomy_list

usage: $taxonomy = $fig->taxonomy_list()

Returns the taxonomy list of all organisms in a hash ref.  Gives the taxonomy down to
genus and species.

=cut

sub taxonomy_list {
    my($self) = @_;
    my $taxonomy = $self->cached('_taxonomy');

    if (keys(%$taxonomy) == 0)
    {
	my $rdbH = $self->db_handle;
	my $relational_db_response = $rdbH->SQL("SELECT genome,taxonomy  FROM genome");
	my $pair;
	foreach $pair (@$relational_db_response)
	{
	    $taxonomy->{$pair->[0]} = $pair->[1];
	}
    }
    return $taxonomy;
}



=head3 is_bacterial

usage: $fig->is_bacterial($genome)

Returns true iff the genome is bacterial.

=cut

sub is_bacterial :Scalar {
    my($self,$genome) = @_;

    return ($self->taxonomy_of($genome) =~ /^\s*Bacteria/o) ? 1 : 0;
}


=head3 is_archaeal

usage: $fig->is_archaeal($genome)

Returns true iff the genome is archaeal.

=cut

sub is_archaeal :Scalar {
    my($self,$genome) = @_;

    return ($self->taxonomy_of($genome) =~ /^\s*Archaea/o) ? 1 : 0;
}


=head3 is_prokaryotic

usage: $fig->is_prokaryotic($genome)

Returns true iff the genome is prokaryotic

=cut

sub is_prokaryotic :Scalar {
    my($self,$genome) = @_;

    return ($self->taxonomy_of($genome) =~ /^\s*(Archaea|Bacteria)/o) ? 1 : 0;
}


=head3 is_eukaryotic

usage: $fig->is_eukaryotic($genome)

Returns true iff the genome is eukaryotic

=cut

sub is_eukaryotic :Scalar {
    my($self,$genome) = @_;

    my $tax = $self->taxonomy_of($genome);
    if (! $tax) { return 0 }
    return ($tax =~ /^\s*Eukaryota/o) ? 1 : 0;
}


=head3 is_viral

usage: $fig->is_viral($genome)

Returns true iff the genome is viral

=cut

sub is_viral :Scalar {
    my($self,$genome) = @_;

    return ($self->taxonomy_of($genome) =~ /^\s*Vir/o) ? 1 : 0;
}


=head3 is_plasmid

usage: $fig->is_plasmid($genome)

Returns true iff the genome is marked as being a plasmid

=cut

sub is_plasmid :Scalar {
    my ($self, $genome) = @_;

    my $taxonomy = $self->taxonomy_of($genome);

    return ($taxonomy =~ m/[Pp]lasmid/o);
}


=head3 is_environmental

usage: $fig->is_environmental($genome)

Returns true if the genome is from an environmental sample

=cut

sub is_environmental :Scalar {
    my($self,$genome) = @_;
    return ($self->taxonomy_of($genome) =~ /environmental samples/io) ? 1 : 0;
}


=head3 sort_genomes_by_taxonomy

usage: @genomes = $fig->sort_genomes_by_taxonomy(@list_of_genomes)

This routine is used to sort a list of genome IDs to put them
into taxonomic order.

=cut

sub sort_genomes_by_taxonomy {
    my($self,@genomes) = @_;

    return map     { $_->[0] }
           sort    { $a->[1] cmp $b->[1] }
           map     { [$_,$self->taxonomy_of($_)] }
           @genomes;
}

=head3 crude_estimate_of_distance

usage: $dist = $fig->crude_estimate_of_distance($genome1,$genome2)

There are a number of places where we need estimates of the distance between
two genomes.  This routine will return a value between 0 and 1, where a value of 0
means "the genomes are essentially identical" and a value of 1 means
"the genomes are in different major groupings" (the groupings are archaea, bacteria,
euks, and viruses).  The measure is extremely crude.

=cut

sub crude_estimate_of_distance :Scalar {
    my($self,$genome1,$genome2) = @_;
    my($i,$v,$d,$dist);

    if ($genome1 > $genome2) { ($genome1,$genome2) = ($genome2,$genome1) }

    my $relational_db_response;
    my $rdbH = $self->db_handle;

    if (($relational_db_response = $rdbH->SQL("SELECT dist FROM distances WHERE ( genome1 = \'$genome1\' ) AND ( genome2 = \'$genome2\' ) ")) &&
        (@$relational_db_response == 1))
    {
        return $relational_db_response->[0]->[0];
    }
    return $self->crude_estimate_of_distance1($genome1,$genome2);
 }

sub crude_estimate_of_distance1 :Scalar {
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

=head3 sort_fids_by_taxonomy

usage: @sorted_by_taxonomy = $fig->sort_fids_by_taxonomy(@list_of_fids)

Sorts a list of feature IDs based on the taxonomies of the genomes that contain the features.

=cut

sub sort_fids_by_taxonomy {
    my($self,@fids) = @_;

    return map     { $_->[0] }
           sort    { $a->[1] cmp $b->[1] }
           map     { [$_,$self->taxonomy_of(&genome_of($_))] }
           @fids;
}


# RAE. Sometimes we want to do the building tree for all genomes, not just complete ones.
# Therefore, I broke this into two sections, one that should retain all the function of
# build_tree_of_complete and the other that does the calculation


sub build_tree_of_complete {
    my($self,$min_for_label) = @_;
    return $self->build_tree_of_all($min_for_label, "complete");
}

sub build_tree_of_all {
    my($self, $min_for_label, $complete)=@_;

    #
    # Find a cached version of the tree if it exists already. We will leak
    # memory if we don't do this, because trees do not deallocate due to circular data structures.
    #

    my $cache = $self->cached('_precomputed_trees');
    my $res = $cache->{$min_for_label, $complete};
    if (!defined($res))
    {
	$res = [$self->build_tree_of_all_real($min_for_label, $complete)];
	$cache->{$min_for_label, $complete} = $res;
	weaken($cache->{$min_for_label, $complete});
    }
    return @$res;
}

sub build_tree_of_all_real {
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
    &limit_labels($tree,$min_for_label);
    unlink("/tmp/tree$$");
    return ($tree,&tips_of_tree($tree));
}

sub get_taxonomy_tree {
    my($self) = @_;

    my $relational_db_response;
    my $rdbH = $self->db_handle;

    if (($relational_db_response = $rdbH->SQL("SELECT genome, taxonomy FROM genome ")) && (@$relational_db_response > 0)) {

      my $tree = {};
      foreach my $element (@$relational_db_response) {
	if ($element->[0] !~ /^99999/) {
	  my @tax_list = map { '{"' . $_ . '"}' } split("; ", $element->[1]);
	  for (my $i=0; $i<scalar(@tax_list); $i++) {
	    my @x = @tax_list;
	    splice(@x, $i + 1);
	    my $a = '$tree->' . join('->', @x);
	    eval 'unless (exists(' . $a . ')) { ' . $a . '= {}; }';
	  }
	}
      }

      return $tree;
    } else {
      return undef;
    }
}

sub limit_labels {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($tree,$min_for_label) = @_;

    my($children) = &tree_utilities::node_pointers($tree);
    if (@$children == 1)
    {
        return 1;
    }
    else
    {
        my $n = 0;
        my $i;
        for ($i=1; ($i < @$children); $i++)
        {
            $n += &limit_labels($children->[$i],$min_for_label);
        }
        if ($n < $min_for_label)
        {
            $tree->[0] = "";
        }
        return $n;
    }
}

sub taxonomic_groups_of_all {
    my($self,$min_for_labels) = @_;

    my($tree,undef) = $self->build_tree_of_all($min_for_labels);
    return &taxonomic_groups($tree);
}

sub taxonomic_groups_of_complete {
    my($self,$min_for_labels) = @_;

    my($tree,undef) = $self->build_tree_of_complete($min_for_labels);
    return &taxonomic_groups($tree);
}

sub taxonomic_groups {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($tree) = @_;

    my($groups,undef) = &taxonomic_groups_and_children($tree);
    return $groups;
}

sub taxonomic_groups_and_children {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($tree) = @_;
    my($ids1,$i,$groupsC,$idsC);

    my $ptrs   = &tree_utilities::node_pointers($tree);
    my $ids    = [];
    my $groups = [];

    if (@$ptrs > 1)
    {
        $ids1 = [];
        for ($i=1; ($i < @$ptrs); $i++)
        {
            ($groupsC,$idsC) = &taxonomic_groups_and_children($ptrs->[$i]);
            if (@$groupsC > 0)
            {
                push(@$groups,@$groupsC);
            }
            push(@$ids1,@$idsC);
        }

        if ($tree->[0])
        {
            push(@$groups,[$tree->[0],$ids1]);
        }
        push(@$ids,@$ids1);
    }
    elsif ($tree->[0])
    {
        push(@$ids,$tree->[0]);
    }

    return ($groups,$ids);
}

################################# Literature Stuff  ####################################

=head2 Literature Methods

=cut

sub get_titles_by_gi {
    my($self,$gi) = @_;

    &verify_existence_of_literature;

    $gi =~ s/^gi\|//;
    my $relational_db_response;
    my $rdbH = $self->db_handle;

    if (($relational_db_response = $rdbH->SQL("SELECT pmid,title FROM literature_titles WHERE ( gi = '$gi' ) ")) &&
        (@$relational_db_response > 0))
    {
        return sort { $a->[1] cmp $b->[1] } @$relational_db_response;
    }
    else
    {
        return ();
    }
}

sub get_titles_by_peg {
    my($self,$peg) = @_;
    my $gi;

    &verify_existence_of_literature;

    my @gis = grep { $_ =~ /^gi\|/ } $self->feature_aliases($peg);
    if (@gis > 0)
    {
        my $relational_db_response;
        my $rdbH = $self->db_handle;
        my $constraint = join(" OR ", map { $gi = ($_ =~ /gi\|(\S+)/) ? $1 : $_;  "( gi = '$gi' )" } @gis);
        if (($relational_db_response = $rdbH->SQL("SELECT pmid,title FROM literature_titles WHERE ( $constraint ) ")) &&
            (@$relational_db_response > 0))
        {
            return sort { $a->[1] cmp $b->[1] } @$relational_db_response;
        }
        else
        {
            return ();
        }
    }
    return ();
}

sub get_title_by_pmid {
    my($self,$pmid) = @_;

    &verify_existence_of_literature;

    $pmid =~ s/^.*\|//;
    my $relational_db_response;
    my $rdbH = $self->db_handle;

    if (($relational_db_response = $rdbH->SQL("SELECT DISTINCT title FROM literature_titles WHERE ( pmid = '$pmid' ) ")) &&
        (@$relational_db_response == 1))
    {
        return $relational_db_response->[0]->[0];
    }
    else
    {
        return "";
    }
}

sub verify_existence_of_literature {

    if (! -d "$FIG_Config::global/Literature")
    {
        mkdir("$FIG_Config::global/Literature",0777);
        system "touch $FIG_Config::global/Literature/gi_pmid_title";
        system "$FIG_Config::bin/load_literature";
    }
}

################################# Subsystems  ####################################

=head3 active_subsystems

    my $ssHash = $fig->active_subsystems($genome, $allFlag);

Get all the subsystems in which a genome is present. The return value is a hash
which maps each subsystem name to the code for the variant used by the specified
genome.

=over 4

=item genome

ID of the genome whose subsystems are desired.

=item allFlag (optional)

If TRUE, all subsystems are returned, with unknown variants marked by a variant
code of C<-1> and iffy variants marked by a code of C<0>. If FALSE or omitted,
only subsystems in which the variant is definitively known are returned. The
default is FALSE.

=back

=cut

sub active_subsystems
{
    my($self, $genome, $all) = @_;

    my $out = {};
    if ($FIG_Config::use_subsystem_estimates)
    {
	$out = $self->active_subsystems_estimate($genome, $all);
    }

    my $allout = $self->active_subsystems_actual($genome, $all);
    $out->{$_} = $allout->{$_} foreach keys %$allout;
    return $out;
}


sub active_subsystems_actual {
    my($self,$genome,$all) = @_;
    my($active,$file,$variant);

    $active = {};

    my $dh;
    if (!opendir($dh, "$FIG_Config::data/Subsystems"))
    {
	warn "Cannot open subsystem dir $FIG_Config::data/Subsystems: $!";
	return $active;
    }

    while (my $ssname = readdir($dh))
    {
	if (open(my $fh, "<", "$FIG_Config::data/Subsystems/$ssname/spreadsheet"))
	{
	    while (<$fh>)
	    {
		if (/^$genome\t(\S+)/)
		{
		    my $variant = $1;
		    last if (!$all && ($variant eq '0' || $variant eq '-1' || $variant eq '*0' || $variant eq '*-1'));
		    $active->{$ssname} = $variant;
		}
	    }
	    close($fh);
	}
    }
    closedir($dh);

#     foreach $_ (`grep \"^$genome\" $FIG_Config::data/Subsystems/*/spreadsheet`)
#     {
#         if (($_ =~ /^(.*?)\/spreadsheet:$genome\t(\S+)/))
#         {
#             next if (!($all) && (($2 eq '0') || ($2 eq '-1') || ($2 eq '*0') || ($2 eq '*-1')));
#             $file = $1;
#             $variant = $2;
#             if ($file =~ /^.*?([^\/]+)$/)
#             {
#                 $active->{$1} = $variant;
#             }
#         }
#     }
    return $active;
}

sub active_subsystems_estimate
{
    my($self, $genome, $all) = @_;

    my $dir = $self->organism_directory($genome);
    my $fh;
    my $sfh;
    my $ret = {};

    if (!open($sfh, "<", "$dir/Subsystems/subsystems"))
    {
	warn "active_subsystems_estimate: No subsystems file found\n";
	return $ret;
    }
    #
    # Read variant codes.
    #
    while (<$sfh>)
    {
	chomp;
	my($ss, $var) = split(/\t/);
	next if (!$all && ($var eq '0' || $var eq '-1' || $var eq '*0' || $var eq '*-1'));
	$ss =~ s/\s+/_/g;
	$ret->{$ss} = $var;
    }
    close($sfh);
    return $ret;
}

=head2 Subsystem Methods

=cut

sub exportable_subsystem {
    my($self,$ssa) = @_;
    my(%seqs,@genomes);

    my $spreadsheet = [];
    my $notes = [];

    $ssa =~ s/[ \/]/_/g;
    my $subsys_dir = "$FIG_Config::data/Subsystems/$ssa";
    if (open(SSA,"<$subsys_dir/spreadsheet"))
    {
        #
        # Push the subsystem metadata.
        #
        my $version = $self->subsystem_version($ssa);
        my $exchangable = $self->is_exchangable_subsystem($ssa);
        push(@$spreadsheet,"$ssa\n$version\n$exchangable\n");

        my $curation = "0000000000\tmaster:unknown\tstarted\n";
	if  (open(my $log, "<","$FIG_Config::data/Subsystems/$ssa/curation.log"))
	{
	    while (defined($_ = <$log>))
	    {
		if (/started/)
		{
		    $curation = $_;
		}
	    }
	    close($log);
        }

        push(@$spreadsheet,$curation,"//\n");

        #
        # Roles
        #

        while (defined($_ = <SSA>) && ($_ !~ /^\/\//))
        {
            push(@$spreadsheet,$_);
        }
        push(@$spreadsheet,"//\n");

        #
        # Subsets
        #

        while (defined($_ = <SSA>) && ($_ !~ /^\/\//))
        {
            push(@$spreadsheet,$_);
        }
        push(@$spreadsheet,"//\n");

        #
        # The spreadsheet itself.
        # Collect the pegs referenced into %seqs.
        #
        while (defined($_ = <SSA>))
        {
            push(@$spreadsheet,$_);
            chomp;
            my @flds = split(/\t/,$_);
            my $genome = $flds[0];
            push(@genomes,$genome);
            my($i,$id);
            for ($i=2; ($i < @flds); $i++)
            {
                if ($flds[$i])
                {
                    my @entries = split(/,/,$flds[$i]);
                    foreach $id (@entries)
                    {
                        my $type = ($id =~ /^(\S+)\.(\d+)$/) ? $1 : "peg";
                        my $n    = ($id =~ /(\d+)$/) ? $1 : "";
                        if ($type && $n)
                        {
                            $seqs{"fig\|$genome.$type.$n"} = 1;
                        }
                    }
                }
            }
        }
        push(@$spreadsheet,"//\n");

        #
        # Assignments and aliases.
        #

        my($fid);
        foreach $fid (sort { &FIG::by_fig_id($a,$b) } keys(%seqs))
        {
            my @aliases = grep { $_ =~ /^(sp\||gi\||pirnr\||kegg\||N[PGZ]_)/ } $self->feature_aliases($fid);

            my $alias_txt = join(",",@aliases);
            my $genome = $self->genome_of($fid);
            my $gs_txt = $self->genus_species($genome);
            my $func_txt = scalar $self->function_of($fid);
            my $location = $self->feature_location($fid);
            my %seen;
            my @checksums = map { [ $_, $self->contig_md5sum( $genome, $_ ) ] }
                            grep { $_ && ( ! $seen{ $_ }++ ) }
                            map  { m/^(\S+)_\d+_\d+$/ }
                            split(/,/,$location);
                            my @loc = split( /,/, $location );
            my $checksum = join(";",map { join(",",@$_) } @checksums);

            push(@$spreadsheet, join("\t", ($fid,
                                            $alias_txt,
                                            $gs_txt,
                                            $func_txt),
                                            $location,
                                            $checksum) . "\n");
        }
        push(@$spreadsheet,"//\n");

        #
        # sequence data
        #

        foreach $fid (sort { &FIG::by_fig_id($a,$b) } keys(%seqs))
        {
            my $aliases = $self->feature_aliases($fid);
            my $seq = (&ftype($fid) eq "peg") ? $self->get_translation($fid) :
                                                $self->dna_seq(&genome_of($fid),
                                                               scalar $self->feature_location($fid));
            push(@$spreadsheet,">$fid $aliases\n");
            my($i,$ln);
            $ln = length($seq);
            for ($i=0; ($i < $ln); $i += 60)
            {
                if (($ln - $i) > 60)
                {
                    push(@$spreadsheet,substr($seq,$i,60) . "\n");
                }
                else
                {
                    push(@$spreadsheet,substr($seq,$i) . "\n");
                }
            }
        }
        close(SSA);

        push(@$spreadsheet,"//\n");

        #
        # Notes file
        #

        if (open(NOTES,"<$FIG_Config::data/Subsystems/$ssa/notes"))
        {
            while (defined($_ = <NOTES>))
            {
                push(@$notes,$_);
            }
            close(NOTES);
        }

        if ($notes->[$#{$notes}] ne "\n") { push(@$notes,"\n") }
        push(@$notes,"//\n");

        #
        # And tag the reactions onto the end. This is fudging the format a little bit, but
        # it should let older parsers handle the subsystems with extra sections.
        #

        if (open(REACTIONS, "<$FIG_Config::data/Subsystems/$ssa/reactions"))
        {
            while (<REACTIONS>)
            {
                push(@$notes, $_);
            }
        }

        #
        # And here we break compatibility. If we have diagrams,
        # save the diagram images.
        #

        if (opendir(D, "$subsys_dir/diagrams"))
        {
            my @ids = grep { not m/^\./ and -d "$subsys_dir/diagrams/$_" } readdir(D);
            closedir(D);

            for my $id (@ids)
            {
                my $ddir = "$subsys_dir/diagrams/$id";

                my $name = &FIG::file_head("$ddir/NAME", 1);
                chomp($name);

                if ($name)
                {
                    push(@$notes, "//diagram:$id:name\t$name\n");
                    push(@$notes, "//end\n");
                }

                #
                # Find the diagram image.
                #

                my @images = <$ddir/diagram.{png,gif,jpg,html}>;

                for my $img_file (@images)
                {
                    if (open(DIAGRAM, "<$img_file"))
                    {
                        my $size = -s DIAGRAM;
                        my $base = basename($img_file);
                        push(@$notes, "//diagram:$id:diagram=$base\t$size\n");
                        my $buf;
                        while (read(DIAGRAM, $buf, 60*57))
                        {
                            my $enc = encode_base64($buf);
                            #
                            # Feh, escape the start of the lines.
                            #
                            $enc =~ s/^/B:/mg;
                            push(@$notes, $enc);
                        }
                        close(DIAGRAM);
                        push(@$notes, "//end\n");
                    }
                }

            }
        }
    }
    return ($spreadsheet,$notes);
}

sub usable_subsystem {
    my($self,$sub,$no_cluster_based) = @_;

    if ($self->is_private_subsystem($sub))
    {
	return 0;
    }

    my $cat = $self->subsystem_classification($sub);
    #print STDERR "$sub is $cat->[0]\n";
    return $self->usable_subsystem_classification($cat, $no_cluster_based);
}

sub usable_subsystem_classification
{
    my($self, $cat, $no_cluster_based) = @_;
    return (defined($cat->[0]) && $cat->[0] ne '' &&
	    ($cat->[0] !~ /experimental/i) &&
	    ((! $no_cluster_based) || ($cat->[0] !~ /Clustering-based/)) &&
	    ($cat->[0] !~ /delete/i)
	   );
}

=head3 is_experimental_subsystem

This states if a subsystem is experimental, what would be the
opposite of usable.

=cut

sub is_experimental_subsystem {
    my ( $self, $sub ) = @_;

    my $is_usable = $self->usable_subsystem( $sub );
    return 0 if ( $is_usable );
    return 1;
}

sub is_cluster_based_subsystem {
    my ( $self, $sub ) = @_;

    my $cat = $self->subsystem_classification($sub);
    return (defined($cat->[0]) && ($cat->[0] ne '') &&
	    ($cat->[0] =~ /Clustering-based/));
}

sub is_exchangable_subsystem :Scalar {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my $ssa = (@_ == 1) ? $_[0] : $_[1];
    $ssa =~ s/[ \/]/_/g;
    if (open(TMP,"<$FIG_Config::data/Subsystems/$ssa/EXCHANGABLE"))
    {
        my $line;
        $line = <TMP>;
        if ($line && ($line =~ /^(\S+)/) && $1)
        {
            return 1;
        }
        close(TMP);
    }
    else {
        return 1;
    }
    return 0;
}

=head3 is_private_subsystem

This states if a subsystem is private, meaning that it cannot be be exported.
This is just the opposite of exchangable.

=cut

sub is_private_subsystem :Scalar {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my $ssa = (@_ == 1) ? $_[0] : $_[1];

    if ( is_exchangable_subsystem( $ssa ) ) {
	return 0;
    }
    return 1;
}

sub all_exchangable_subsystems {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);

    my @exchangable = ();
    if (opendir(SUB,"$FIG_Config::data/Subsystems"))
    {
        push(@exchangable,grep { ($_ !~ /^\./) && &is_exchangable_subsystem($_) } readdir(SUB));
        closedir(SUB);
    }
    return @exchangable;
}


=head3 nmpdr_subsystem

Gets and sets whether the subsystem should be published with the NMPDR. Specifically writes a file called NMPDR in the subsystem directory.

Use:

$fig->nmpdr_subsystem($ssa, 1); # to set it as an nmpdr subsystem
$fig->nmpdr_subsystem($ssa, -1); # to set it as NOT an nmpdr subsystem
$fig->nmpdr_subsystem($ssa); # to test whether it is an nmpdr subsystem


=cut

sub nmpdr_subsystem {
    my ($self, $ssa, $nmpdr)=@_;
    if (defined $nmpdr && $nmpdr > 0)
    {
        open(OUT, ">".$FIG_Config::data."/Subsystems/$ssa/NMPDR") || die "Can't write to ". $FIG_Config::data."/Subsystems/$ssa/NMPDR";
        print OUT $ssa;
        close OUT;
        return 1;
    }
    elsif (defined $nmpdr && $nmpdr < 0)
    {
        unlink($FIG_Config::data."/Subsystems/$ssa/NMPDR");
        return 0;
    }

    if (-e $FIG_Config::data."/Subsystems/$ssa/NMPDR") {return 1}
    else {return 0}
}

=head3 distributable_subsystem

Gets and sets whether the subsystem is freely distributable and should be included in new releases.

Use:

$fig->distributable_subsystem($ssa, 1); # to set it as a distributable subsystem
$fig->distributable_subsystem($ssa, -1); # to set it as NOT a distributable subsystem
$fig->distributable_subsystem($ssa); # to test whether it is a distributable subsystem


=cut

sub distributable_subsystem {
    my ($self, $ssa, $distributable)=@_;
    if (defined $distributable && $distributable > 0)
    {
        open(OUT, ">".$FIG_Config::data."/Subsystems/$ssa/DISTRIBUTE") || die "Can't write to ". $FIG_Config::data."/Subsystems/$ssa/DISTRIBUTE";
        print OUT $ssa;
        close OUT;
        return 1;
    }
    elsif (defined $distributable && $distributable < 0)
    {
        unlink($FIG_Config::data."/Subsystems/$ssa/DISTRIBUTE");
        return 0;
    }

    if (-e $FIG_Config::data."/Subsystems/$ssa/DISTRIBUTE") {return 1}
    else {return 0}
}




=head3 all_subsystems

    my @names = $fig->all_subsystems();

Return a list of all of the subsystems in the data store.

=cut

sub all_subsystems {

    my($self) = @_;

    if (ref($self) =~ /^FIG/ && $self->{sdata})
    {
	return $self->{sdata}->all_subsystems();
    }

    my @subsystems = ();
    if (opendir(SUB,"$FIG_Config::data/Subsystems"))
    {
        push(@subsystems,grep { ($_ !~ /^\./) } readdir(SUB));
        closedir(SUB);
    }
    return @subsystems;
}

sub all_subsystems_detailed {

    my($self) = @_;

    return $self->{sdata}->all_subsystems_detailed();
}

sub subsystem_metadata_update
{
    my($self, $name, $field, $value) = @_;
    my %fields = map { $_ => 1 } qw( classification class_1 class_2 curator creation_date
				    last_update version exchangable);
    $fields{$field} or confess "subsystem_metadata_update: invalid field $field";
    my $dbh = $self->db_handle->{_dbh};
    eval {
	local($dbh->{RaiseError})=1;
	my $n = $dbh->do(qq(UPDATE subsystem_metadata
		    SET $field = ?
		    WHERE subsystem = ?),
		 undef, $value, $name);
    };
    $self->flush_subsystem_cache();
}

sub flush_subsystem_cache
{
    my($self) = @_;
    return $self->{sdata}->flush_cache();
}

sub variant_code {
    my($self,$subsys,$genome) = @_;
    my $subsystem = new Subsystem( $subsys, $self, 0 );
    if ($subsystem && (my $idx = $subsystem->get_genome_index($genome)))
    {
	return $subsystem->get_variant_code($idx);
    }
    else
    {
	print "$subsys is probably not a subsystem or $genome is not in it\n";
    }
}

=head3 all_usable_subsystems

    my @names = $fig->all_usable_subsystems();

Return a list of all of the subsystems in the data store that are "usable", that is,
not experimental or deleted.

Use the subsystem information cache if valid.

=cut

sub all_usable_subsystems
{
    my($self) = @_;

    my $cache = $self->get_valid_cache_file("subsys/usable_subsystems");

    my @subsystems = ();

    if ($cache)
    {
	#warn "Reading from cache\n";
	while (<$cache>)
	{
	    chomp;
	    push(@subsystems, $_);
	}
	$cache->close();
    }
    else
    {
	#warn "reading from dir\n";
	if (opendir(SUB,"$FIG_Config::data/Subsystems"))
	{
	    push(@subsystems, grep { ($_ !~ /^\./) and $self->usable_subsystem($_) } readdir(SUB));
	    closedir(SUB);
	}
    }
    return @subsystems;
}




=head3 index_subsystems

Run indexing on one or more subsystems. If no subsystems are defined we will reindex the whole thing. Otherwise we will only index the defined subsystem. Note that this method just launches index_subsystems as a background job. Returns the job of the child process.

$pid=$fig->index_subsystems("Alkanesulfonates Utilization"); # do only Alkanesulfonates Utilization
$pid=$fig->index_subsystems(@ss); # do subsystems in @ss
$pid=$fig->index_subsystems(); # do all subsystems

=cut

sub index_subsystems {
 my ($self, @ss)=@_;
 print STDERR "Trying $FIG_Config::bin/index_subsystems @ss\n";
 return $self->run_in_background(
  sub {
   my $cmd="$FIG_Config::bin/index_subsystems @ss";
   print "Will run '$cmd'\n";
   &run($cmd);
   print "finished.\n";
   }
 );
}

sub delete_subsystem {
    my($self,$sub) = @_;

    if (my $subF = $self->get_subsystem($sub))
    {
	$subF->delete_indices();
	&FIG::run("rm -r \"$FIG_Config::data/Subsystems/$sub\"");
    }
}

sub rename_subsystem {
    my($self,$from,$to) = @_;

    if (my $subF = $self->get_subsystem($from))
    {
	$subF->delete_indices();
	if (rename("$FIG_Config::data/Subsystems/$from","$FIG_Config::data/Subsystems/$to"))
	{
	    &FIG::run("$FIG_Config::bin/index_subsystems \"$to\"");
	}
    }
}


=head3 perform_subsystem_salvage

    my $glist = [['273035.1', '273035.4']];
    my $pmap = { 'fig|273035.1.peg.1' => 'fig|273035.4.peg.4', ... };
    $fig->perform_subsystem_salvage($glist, $pmap);

For each subsystem in this SEED, perform a subsystem salvage operation for each old-genome / new-genome pair in $glist.
This operation will determine if the old genome exists in the subsystem. If it does, the new genome is
added to the subsystem, and we attempt to map the pegs from the cells in the old subsystem's row to the new
subsystem. If all pegs map, we copy the variant code for the genome. If all cells did not map, we prepend
a * to the variant code before copying.

=cut

sub perform_subsystem_salvage
{
    my($fig, $genome_pairs, $map) = @_;

    for my $ssname ($fig->all_subsystems())
    {
	my $ss = new Subsystem($ssname, $fig);

	if (!defined($ss))
	{
	    warn "Subsystem $ssname not found during perform_subsystem_salvage()\n";
	    next;
	}

	for my $gpair (@$genome_pairs)
	{
	    my($g, $ng) = @$gpair;

	    my $idx = $ss->get_genome_index($g);
	    if (defined($idx))
	    {
		print "Salvaging $ssname for old_genome=$g new_genome=$ng\n";
		my $row = $ss->get_row($idx);

		my $new_idx = $ss->get_genome_index($ng);

		if (defined($new_idx))
		{
		    warn "$ng was already present as $new_idx\n";
		    next;
		}

		$new_idx = $ss->add_genome($ng);

		if (!defined($new_idx))
		{
		    die "Subsystem $ss add_genome($ng) failed";
		}

		my @nrow;
		my $pegcount = 0;
		my $mapcount = 0;
		for (my $ci = 0; $ci < @$row; $ci++)
		{
		    my $c = $row->[$ci];

		    my @nc;

		    for my $p (@$c)
		    {
			$pegcount++;
			my $new = $map->{$p};
			if ($new)
			{
			    $mapcount++;
			    push(@nc, $new);
			}
		    }

		    $ss->set_pegs_in_cell($new_idx, $ci, \@nc);
		}

		my $vc = $ss->get_variant_code($idx);
		#
		# Always put the star in.
		#
		#if ($mapcount < $pegcount)
		{
		    $vc = "*$vc";
		}
		$ss->set_variant_code($new_idx, $vc);
		$ss->write_subsystem(1);
	    }
	}
    }
}



=head3 all_constructs

Hmmm...

=cut

sub all_constructs {
    my($self) = @_;

    my @subsystems = ();
    if (opendir(SUB,"$FIG_Config::data/Subsystems"))
    {
        push(@subsystems,grep { ($_ !~ /^\./) } readdir(SUB));
        closedir(SUB);
    }

    my @c;
    for my $subname (@subsystems)
    {
        $subname =~ s/[ \/]/_/g;
        my $cfile = "$FIG_Config::data/Subsystems/$subname/constructs";
        if (-f $cfile)
        {
            my $sub = $self->get_subsystem($subname);
            my @a = Construct::parse_constructs_file($cfile, $sub);
            my $l = [];

            for my $con (@a)
            {
                my($cname, $list) = @$con;
                my $nreqs = [];

                for my $req (@$list)
                {
                    if ($req->[0] eq 'R')
                    {
                        push(@$nreqs, ['R', $req->[2]]);
                    }
                    else
                    {
                        push(@$nreqs, $req);
                    }
                }
                push(@$l, [$cname, $nreqs]);
            }
            push(@c, [$subname, $l]);
        }
    }
    return @c;
}

=head3 subsystem_version

 my $version=subsystem_version($subsystem_name)

Returns the current version of the subsystem.

=cut

sub subsystem_version :Scalar {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my $ssa = (@_ == 1) ? $_[0] : $_[1];
    $ssa =~ s/[ \/]/_/g;

    if (open(VER,"<$FIG_Config::data/Subsystems/$ssa/VERSION"))
    {
        my $ver = <VER>;
        close(VER);
        if ($ver =~ /^(\S+)/)
        {
            return $1;
        }
    }
    return 0;
}

=head3 subsystem_classification

 Get or set the classification of the subsystem. Added by RAE in response to the changes made on seed wiki
 If a reference to an array is supplied it is saved as the new classification of the subsystem.
 Regardless, the current classification is returned as a reference to an array. There is no control over what the things are.
 Returns a reference to an empty array if a valid subsystem is not supplied, or if no classification is known

 The classification is stored as a \t separated list of things in $subsys/CLASSIFICATION. There is no control over what the things are.


=cut

sub subsystem_classification {
 my ($self, $ssa, $classification)=@_;
 $ssa =~ s/[ \/]/_/g;

 my $return=['', ''];

 if ($ssa && $classification->[0]) {
    return $return unless (-e "$FIG_Config::data/Subsystems/$ssa/");
    if (open(SSA,">$FIG_Config::data/Subsystems/$ssa/CLASSIFICATION")) {
     print SSA join("\t", @$classification), "\n";
    }
    close SSA;
    return $classification;
 }

 return $self->{sdata}->subsystem_classification($ssa);
 #
 # See if we have populated the subsystem classification cache, and use that if so.
 #
 my $ss_cache = $self->cached('_ss_classification');
 if (my $res = $ss_cache->{$ssa})
 {
     return [split(/\t/, $res)];
 }

 # using get_subsystem is really slow, and so we are going to cat the file and return that

 #return $subsys->get_classification;
 if (open(SSA,"<$FIG_Config::data/Subsystems/$ssa/CLASSIFICATION")) {
    my @line;
    while (my $x = <SSA>) {
     chomp $x;
     my @thisline=split(/\t/,$x);
     if ($thisline[0] || $thisline[1]) {@line=@thisline}
    }
    close(SSA);
    $line[0]='' unless (defined $line[0]);
    $line[1]='' unless (defined $line[1]);
    return [$line[0], $line[1]];
 }
 else
 {
    return ['', ''];
 }
}

=head3 all_subsystem_classifications

    my @classifications = $fig->all_subsystem_classifications();

Return a list of all the subsystem classifications. Each element in the
list will contain a main subsystem class and a basic subsystem class.
The resulting list enables us to determine easily what the three-level
subsystem tree would look like.

=cut

sub all_subsystem_classifications {
    my $self=shift;

    #
    # Ensure we have cached subsystem metadata.
    #

    my $ss_meta = $self->load_subsystem_classification();

    my %out = map { $_ => 1 } values %$ss_meta;
    return map { [split(/\t/)] } keys %out;
}

sub load_subsystem_classification
{
    my($self) = @_;

    #
    # Is it cached?
    #

    my $meta;
    if (! $FIG_Config::disable_subsystem_classification_cache)
    {
	$meta = $self->cached("_ss_classification");
	return $meta if %$meta;
    }

    #
    # Try the subsystem_metadata table first.
    #

    my $res;
    eval {
	my $rdbH = $self->db_handle;
	my $dbh = $rdbH->{_dbh};

	local($dbh->{RaiseError})=1;

	my $expt = "";
	if ($FIG_Config::exclude_experimental_subsystems)
	{
	    $expt = " WHERE class_1 NOT LIKE 'Experimental%'";
	}

	$res = $dbh->selectall_arrayref(qq(SELECT subsystem, classification
					      FROM subsystem_metadata $expt));
    };
    if ($res)
    {
	map { $meta->{$_->[0]} = $_->[1] } @$res;
	return $meta;
    }

    for my $ss ($self->all_subsystems)
    {
	my $cl=join "\t", @{$self->subsystem_classification($ss)};
	$meta->{$ss} = $cl;
    }
    return $meta;
}

=head3 subsystem_curator

usage: $curator = $fig->subsystem_curator($subsystem_name)

Return the curator of a subsystem.

=cut

sub subsystem_curator :Scalar {
    my($self, $ssa) = @_;
    my($who) = "";

    $ssa =~ s/[ \/]/_/g;

    return $self->{sdata}->curator($ssa);

    if (open(DATA,"<$FIG_Config::data/Subsystems/$ssa/curation.log"))
    {
        while (defined($_  = <DATA>))
        {
            if ($_ =~ /^\d+\t(\S+)\s+started/)
            {
                $who = $1;
            }
        }
        close(DATA);
    }
    $who =~ s/^master://i;
    return $who;
}

sub reset_subsystem_curator :Scalar {
    my($self, $ssa, $who) = @_;

    $ssa =~ s/[ \/]/_/g;

    if ($who && open(LOG,">>$FIG_Config::data/Subsystems/$ssa/curation.log"))
    {
        $who = ($who =~ /^master:/) ? $who : "master:$who";
        my $time = time;
        print LOG "$time\t$who\tstarted\n";
        close(LOG);
	#
	# Be sure to load ss AFTER update to ensure we catch the new value.
	#
	my $ss = Subsystem->new($ssa, $self);
	$ss->db_sync();
	$self->flush_subsystem_cache();
        return 1;
    }
    return 0;
}

=head3 subsystem_info

Returns the number of diagrams of the passed subsystem.

=cut

sub subsystem_num_new_diagrams {
  my($self,$ssa) = @_;

  my $diag_dir = "$FIG_Config::data/Subsystems/$ssa/diagrams";
  if (opendir(DIR, $diag_dir))
    {
      my @diagrams = grep { /^d/ && -d "$diag_dir/$_" } readdir(DIR);
      my $counter = 0;
      closedir DIR;
      foreach my $d ( @diagrams ) {
         my $image_map = "$diag_dir/$d/diagram.html";
         if ($image_map) {
            if ( open(IN, "$image_map") ) {
               my $header = <IN>;
               close(IN);

               if ($header =~ /\<map name=\"GraffleExport\"\>/) {
                  $counter++;
               }
            }
         }
      }
      return $counter;
    }
  else
    {
      return 0;
    }
}


sub subsystem_num_diagrams {
  my($self,$ssa) = @_;

  my $diag_dir = "$FIG_Config::data/Subsystems/$ssa/diagrams";
  if (opendir(DIR, $diag_dir))
    {
      my @diagrams = grep { /^d/ && -d "$diag_dir/$_" } readdir(DIR);
      closedir DIR;
      return scalar(@diagrams);
    }
  else
    {
      return 0;
    }
}
=head3 subsystem_info

usage: ($version, $curator, $pedigree, $roles) = $fig->subsystem_info($subsystem_name)

Return information about the given subsystem.

$roles is a list of tuples (abbrev, name).


=cut

sub subsystem_info {
    my($self,$ssa) = @_;
    my($version, $curator, $pedigree, $roles);;

    $ssa =~ s/[ \/]/_/g;

    $roles = [];

    $version = $self->subsystem_version($ssa);
    $curator = $self->subsystem_curator($ssa);

    if (open(CUR, "<$FIG_Config::data/Subsystems/$ssa/curation.log"))
    {
        local($/);
        $pedigree = <CUR>;
    }

    if (open(SSA,"<$FIG_Config::data/Subsystems/$ssa/spreadsheet"))
    {
        #
        # The spreadsheet appears to be of the form
        #
        #   role-abbr role name
        #   ...
        #   //
        #   Something about subsets
        #   //
        #   genome-id spreadsheet-info
        #

        local $/ = "//";

        my($chunk);

        if (defined($chunk = <SSA>))
        {
            for $_ (split(/\n/, $chunk))
            {
                chomp;
                if (/^(\S[^\t]*)\s+(\S.*\S)\s*$/)
                {
                    push(@$roles, [$1, $2]);
                }
            }
        }
        close(SSA);
    }

    return ($version, $curator, $pedigree, $roles);
}

sub curation_history {
    my ( $self, $ssa ) = @_;
    $ssa =~ s/[ \/]/_/g;
    my $curhash;

    if (open(CUR, "<$FIG_Config::data/Subsystems/$ssa/curation.log")) {
        while ( <CUR> ) {
            my ( $ts, $cur, $what ) = split( "\t", $_ );
            $curhash->{ $ts }->{ 'curator' } = $cur;
            $curhash->{ $ts }->{ 'what' } = $what;
        }
    }
    return $curhash;
}

=head3 subsystems_for_genome

usage: @subsystems = $fig->subsystems_for_genome($genome, $all)

Return the list of subsystems in which the genome has been entered.

@subsystems is a list of subsystem names.

It will only return those genomes with a variant code other than 0 or -1,
unless the $all argument is "true" (in which case all subsystems are returned).

If $all is 2 then it will return all subsystems with a variant code other than -1.

=cut
#: Return Type $@@;


sub subsystems_for_genome {
    my($self,$genome, $all) = @_;

    if (! $self->is_genome($genome)) { return () }
    my $rdbH = $self->db_handle;


    # There are some legacy seed instances lacking the variant field in subsystem_index, so
    # trap that error and return an empty list.

    my $subsystem_data;

    {
        my $dbh = $rdbH->{_dbh};
        local $dbh->{RaiseError} = 1;
        local $dbh->{PrintError} = 0;

        my $sql="SELECT DISTINCT subsystem from subsystem_index WHERE (protein like 'fig\|$genome.peg.%'";
        if (defined($all) && ($all == 2)) {
	  $sql .= " AND (variant != '-1')";
	  $sql .= " AND (variant != '*-1')";
	} elsif (!$all) {$sql .= " AND (variant != '-1' AND variant != '0' AND variant != '*-1' AND variant != '*0')"}
        $sql .= ")";

        eval {
            $subsystem_data = $rdbH->SQL($sql);
        };
    }

    if ($@ =~ /variant/) {
        return [];
    }

    return  map { $_->[0] } @$subsystem_data;
}


=head3 subsystem_genomes

usage: $genomes = $fig->subsystem_genomes($subsystem_name, $all)

Return the list of genomes in the subsystem.

$genomes is a list of tuples (genome_id, name)

unless ($all) is set to true it will only return those genomes with a variant code other thaN
0 OR -1.

=cut
#: Return Type $@@;

sub subsystem_genomes :Scalar {
    my($self,$ssa,$all) = @_;
    my $fileName = "$FIG_Config::data/Subsystems/$ssa/spreadsheet";
    my $genomes = $self->readSpreadsheetForGenomes($fileName, $all);
    return $genomes;
}

=head3 readSpreadsheetForGenomes

    my $genomeList = $fig->readSpreadsheetForGenomes($fileName, $all);

Read the genomes from a specific subsystem file. This allows the client to get
the genome data for a backup subsystem.

=over 4

=item fileName

Name of the subsystem spreadsheet file.

=item all

If TRUE, all genomes will be read. Otherwise, only those genomes with a specific variant code
(i.e. not 0 or -1) will be returned.

=item RETURN

Returns a reference to a list of 2-tuples, each consisting of a genome ID and the genome's name.

=back

=cut

sub readSpreadsheetForGenomes {
    my ($self, $fileName, $all) = @_;
    my($genomes);

    $genomes = [];

    if (open(SSA,"<$fileName"))
    {
        #
        # The spreadsheet appears to be of the form
        #
        #   role-abbr role name
        #   ...
        #   //
        #   Something about subsets
        #   //
        #   genome-id spreadsheet-info
        #

        local $/ = "//";

        my($chunk);

        if (defined($chunk = <SSA>))
        {
        }
        if (defined($chunk = <SSA>))
        {
        }
        local $/ = "\n";

        while (<SSA>)
        {
            chomp;
            s/^\s*//;
            s/\s*$//;
            next if $_ eq "";
            if (($_ =~ /^(\d+\.\d+)\s+(\S+)/) && ($all || ($2 && ($2 ne "-1"))))
            {
                my $genome = $1;
                if ($self->is_genome($genome))
                {
                    my $name = $self->genus_species($genome);
                    push(@$genomes, [$genome, $name]);
                }
            }
        }
        close(SSA);
    }

    return $genomes;
}

#
#    @pegs              = $fig->pegs_in_subsystem_cell($subsystem, $genome,$role)
#    @roles              = $fig->subsystem_to_roles($subsystem)
#    @maps             = $fig->role_to_maps($role)
#    @subsystems = $fig->peg_to_subsystems($peg);

=head3 get_subsystem

    my $subsysObject = $fig->get_subsystem($name, $force_load);

Return a subsystem object for manipulation of the named subsystem. If the
subsystem does not exist, an undefined value will be returned.

=over 4

=item name

Name of the desired subsystem.

=item force_load

TRUE to reload the subsystem from the data store even if it is already cached in
memory, else FALSE.

=item RETURN

Returns a blessed object that allows access to subsystem data, or an undefined
value if the subsystem does not exist.

=back

=cut

sub get_subsystem :Scalar
{
    my($self, $subsystem, $force_load) = @_;
    my $sub;

    $subsystem =~ s/[ \/]/_/g;
    my $cache = $self->cached('_Subsystems');
    if ($force_load || !($sub = $cache->{$subsystem}))
    {
        $sub = new Subsystem($subsystem, $self);
	if ($sub)
	{
	    $cache->{$subsystem} = $sub;
# 	    weaken $cache->{$subsystem};
	}
    }
    return $sub;
}

=head3 clear_subsystem_cache

    $fig->clear_subsystem_cache();

Delete all subsystems from the subsystem cache. This is not normally
needed, because the cache is kept fairly small. However, in cases where
all of the subsystems are needed, the cache grows by more than a
gigabyte, and because the subsystems point back to the FIG object, the
memory is not cleaned up properly. Calling this mehtod before you release
the FIG object removes that problem.

=cut

sub clear_subsystem_cache {
    # Get the parameters.
    my ($self) = @_;
    # Get the subsystem cache.
    my $cache = $self->cached('_Subsystems');
    # Loop through the cached subsystems, deleting them.
    # Each one we delete subtracts one from the reference counter
    # on FIG, thus making our memory safer.
    for my $sub (keys %$cache) {
        delete $cache->{$sub};
    }
}

sub clear_caches
{
    my($self) = @_;
    delete $self->{gdata};
    delete $self->{sdata};
    $self->clear_subsystem_cache();
}

=head3 subsystem_to_roles

    my @roles = $fig->subsystem_to_roles($subsysID);

Return a list of the roles for the specified subsystem.

=over 4

=item subsysID

Name (ID) of the subsystem whose roles are to be listed.

=item RETURN

Returns a list of roles.

=back

=cut

sub subsystem_to_roles
{
    my($self, $subsystem) = @_;
    $subsystem =~ s/[ \/]/_/g;

    my $sub = $self->get_subsystem($subsystem);

    return () unless $sub;

    return $sub->get_roles();
}

sub is_aux_role_in_subsystem {
    my($self,$subsystem,$role) = @_;

    my $rdbH = $self->db_handle();
    if (! $rdbH->table_exists('aux_roles')) { return 0 }
    $subsystem =~ s/\s/_/g;
    my $nameQ = quotemeta $subsystem;
    my $roleQ = quotemeta $role;
    my $q = "SELECT subsystem FROM aux_roles WHERE subsystem = '$nameQ' AND role = '$roleQ'";

    my $relational_db_response;
    return (($relational_db_response = $rdbH->SQL($q)) && (@$relational_db_response > 0));

    # my $subO = $self->get_subsystem($subsystem);

    # return $subO ? $subO->is_aux_role($role) : 0;
}

sub pegs_in_subsystem_cell
{
    my($self, $subsystem, $genome, $role) = @_;
    $subsystem =~ s/[ \/]/_/g;

    my $sub = $self->get_subsystem($subsystem);

    return undef unless $sub;
    return grep { ! $self->is_deleted_fid($_) } $sub->get_pegs_from_cell($genome, $role);
}

sub get_clearinghouse :Scalar
{
    my($self, $url) = @_;

    if (defined($self->{_clearinghouse}))
    {
        return $self->{_clearinghouse};
    }

    if (!$ClearinghouseOK)
    {
        warn "Error: Clearinghouse code not available.\n";
        return undef;
    }

    if ($url eq "")
    {
	# $url = "http://www.mcs.anl.gov/~olson/SEED/api.cgi";
	$url = "http://pubseed.theseed.org/legacy_clearinghouse/api.cgi";
    }

    my $ch = new Clearinghouse($url);
    $self->{_clearinghouse} = $ch;

    return $ch;
}

sub publish_subsystem_to_clearinghouse
{
    my ($self, $ssa, $url, $newsys) = @_;
    my ($id, $token);

    $ssa =~ s/[ \/]/_/g;

    my $ch = $self->get_clearinghouse($url);

    if (!defined($ch))
    {
        warn "Cannot publish: clearinghouse not available\n";
        return undef;
    }

    my($version, $curator, $pedigree, $roles) = $self->subsystem_info($ssa);

    my $genomes = $self->subsystem_genomes($ssa);

    my @genome_names = ();
    for my $g (@$genomes)
    {
        push(@genome_names, $g->[1]);
    }

    my $seed_id = $self->get_seed_id();

    my $time = int(time());

    my $returnlog = '';

    if ( !$newsys ) {
        print "publishing: ss=$ssa version=$version time=$time curator=$curator seed_id=$seed_id\n";
    }
    else {
        $returnlog .= "publishing: ss=$ssa version=$version time=$time curator=$curator seed_id=$seed_id<BR>";
    }
    my $ret = $ch->publish_subsystem($ssa, $version, $time, $curator, $pedigree, $seed_id,
                                     $roles, \@genome_names);

    ($id, $token, $url) = @$ret;
    if ( !$newsys ) {
        print "Got id  $id token $token url $url\n";
    }
    else {
        $returnlog .= "Got id  $id token $token url $url<BR>";
    }


    #
    # Retrieve the package
    #

    if ( !$newsys ) {
        print "Packaging...\n";
    }
    else {
        $returnlog .= "Packaging...<BR>";
    }

    my($spreadsheet, $notes) = $self->exportable_subsystem($ssa);
    my $package = join("", @$spreadsheet, @$notes);
    if ( !$newsys ) {
        print "Sending...\n";
    }
    else {
        $returnlog .= "Sending...<BR>";
    }
    $ch->upload_subsystem_package($url, $package);

    if ( !$newsys ) {
        return 1;
    }
    else {
        return $returnlog;
    }
}

#
# Feh - for credentials handling it's easier to set up subclass of LWP::UserAgent.
#

{
    package FigUserAgent;
    use base 'LWP::UserAgent';

    sub new
    {
	my($class, $user, $pass, @rest) = @_;

	my $self = LWP::UserAgent->new(@rest);
	$self->{_fig_saved_creds} = [$user, $pass];
	return bless $self, $class;
    }

    sub get_basic_credentials
    {
	my($self, $realm, $uri, $isproxy) = @_;
	return @{$self->{_fig_saved_creds}};
    }
}

=head3 install_subsystem_directory_on_server

Install the given local subsystem directory on the SEED at the URL
provided. If authentication is required, the given username and
password will be used.

Uses an HTTP POST of the tarfile of the contents of the local directory to the
install_subsystem_dir.cgi CGI script.

=cut

sub install_subsystem_directory_on_server
{
    my($self, $dir, $server_url, $username, $password) = @_;

    my $url = "$server_url/install_subsystem_dir.cgi";

    if (! -d $dir)
    {
	die "Subsystem directory $dir does not exist";
    }
    if (! -f "$dir/spreadsheet")
    {
	die "Subsystem directory $dir does not appear to contain a subsystem";
    }

    my $ssa = basename($dir);

    #
    # Create compressed tarfile.
    #
    my $tarfile = "$FIG_Config::temp/subsys.$$.tgz";
    &run("tar -c -z -f $tarfile -C $dir .");

    my $form = [ssa => $ssa,
		tarfile => [$tarfile]];


    my $ua = new FigUserAgent($username, $password);
    my $res = $ua->post($url, $form, 'Content-type' => 'form-data');

    unlink($tarfile);

    if ($res->is_success)
    {
	warn "Successful post: " . $res->content . "\n";
	return;
    }

    die "Failure posting request: " . $res->status_line . "\n" . $res->content;
}


sub all_subsystems_with_roles
{
    my($self) = @_;

    my $rdbH = $self->db_handle;

    my $q;

    if ($FIG_Config::exclude_experimental_subsystems)
    {
	$q = qq(SELECT DISTINCT i.subsystem, i.role
		FROM subsystem_index i JOIN subsystem_metadata m ON i.subsystem = m.subsystem
		WHERE m.class_1 NOT LIKE 'Experimental%');
    }
    else
    {
	$q = "SELECT DISTINCT i.subsystem, i.role  FROM subsystem_index i";
    }
    
    my %seen;

    if (my $relational_db_response = $rdbH->SQL($q))
    {
        my $pair;
        foreach $pair (@$relational_db_response)
        {
            my $key = $pair->[0];
	    push (@{$seen{$key}}, $pair->[1]);
	}
    }
    return \%seen;
}


#
# Return the list of subsystems this peg appears in.
# Each entry is a pair [subsystem, role].
#

=head3 subsystems_for_peg

 Return the list of subsystems and roles that this peg appears in.
 Returns an array. Each item in the array is
 a reference to a tuple of subsystem and role.  If the second argument ($noaux)
 is "true", only roles playing non-auxiliary roles will be returned.
 If the third argument ($active_only), only roles where the peg's row
 has a non-0 and non-"-1" variant are returned.

=cut

sub subsystems_for_peg
{
    my($self, $peg,$noaux, $active_only) = @_;

    if ($self->is_deleted_fid($peg)) { return () }

    ($peg =~ /^fig\|\d+\.\d+\.\w+.\d+$/) or return;

    my $rdbH = $self->db_handle;

    my $q = "SELECT subsystem, role, variant FROM subsystem_index WHERE protein = ?";
    if ($active_only)
    {
	$q .= " AND (variant != '-1' and variant != '0' and variant != '*-1' and variant != '*0' ) ";
    }

    if (my $relational_db_response = $rdbH->SQL($q, undef, $peg))
    {
        my %seen;
        my @in;
        my $pair;
        foreach $pair (@$relational_db_response)
        {
            $pair->[0] =~ s/ /_/g;
            my $key = join("\t",@$pair);
            if (! $seen{$key})
            {
		if ((!$active_only) || ($pair->[2] !~ /^\*?(0|-1)$/))
		{
		    push(@in,$pair);
		    $seen{$key} = 1;
		}
            }
        }

	if ($noaux)
	{
	    #
	    # Note that we can handle this case more directly by joining to the aux_roles
	    # table
	    # SELECT i.subsystem, i.role
	    # FROM subsystem_index i LEFT JOIN aux_roles a ON i.role = a.role
	    # WHERE protein = ? AND a.subsystem IS NULL
	    #

	    my @nonaux = ();
	    foreach my $x (@in)
	    {
		if (! $self->is_aux_role_in_subsystem($x->[0],$x->[1]))
		{
		    push(@nonaux,$x);
		}
	    }
	    return @nonaux;
	}
	else
	{
	    return @in;
	}
    }
    else
    {
        return ();
    }
}

=head3 subsystems_for_peg_complete

 Return the list of subsystems that this peg appears in.
 Returns an array. Each item in the array is
 a reference to a tuple of subsystem, role, variant and is_auxiliary.

=cut

sub subsystems_for_peg_complete
{
    my($self, $peg) = @_;

    if ($self->is_deleted_fid($peg)) { return () }

    ($peg =~ /^fig\|\d+\.\d+\.\w+\.\d+$/) or return;

    my $rdbH = $self->db_handle;

    my $q = "SELECT subsystem, role, variant FROM subsystem_index WHERE protein = '$peg'";

    if (my $relational_db_response = $rdbH->SQL($q))
    {
        my %seen;
        my @in;
        my $pair;
        foreach $pair (@$relational_db_response)
        {
            $pair->[3] = $self->is_aux_role_in_subsystem($pair->[0],$pair->[1]);
            $pair->[0] =~ s/ /_/g;
            my $key = join("\t",@$pair);
            if (! $seen{$key})
            {
                push(@in,$pair);
                $seen{$key} = 1;
            }
        }
        return @in;
    }
    else
    {
        return ();
    }
}

=head3 subsystems_for_pegs_complete

 Return the list of subsystems, roles and variants that the pegs appear in.
 Returns a hash keyed by peg. Each item in the hash is a reference to a tuple
 of subsystem, role and variant. If the last argument ($include_aux)
 is "true", also roles playing auxiliary roles will be returned.

=cut

sub subsystems_for_pegs_complete {
    my ($self, $pegs, $include_aux) = @_;

    my $rdbH = $self->db_handle;
    my %results;

    my $q = "SELECT subsystem, role, variant, protein FROM (SELECT subsystem, role, variant, protein FROM subsystem_index WHERE protein IN ('".join("', '", @$pegs)."')) AS t1 LEFT JOIN deleted_fids ON t1.protein=deleted_fids.fid WHERE deleted_fids.fid IS NULL";
    unless ($include_aux) {
	$q = "SELECT t2.subsystem, t2.role, t2.variant, t2.protein FROM (".$q.") AS t2 LEFT JOIN aux_roles ON (t2.subsystem=aux_roles.subsystem AND t2.role=aux_roles.role) WHERE aux_roles.role IS NULL";
    }
    my $ret = $rdbH->SQL($q);

    foreach my $row (@$ret) {
	if (exists($results{$row->[3]})) {
	    push(@{$results{$row->[3]}}, [ $row->[0], $row->[1], $row->[2] ]);
	} else {
	    $results{$row->[3]} = [ [ $row->[0], $row->[1], $row->[2] ] ];
	}
    }

    return %results;
}

=head3 subsystems_for_peg_complete_estimate

 Return the list of subsystems and roles that this peg appears in.
 Returns an array. Each item in the array is
 a reference to a tuple of subsystem and role.
 Data is taken from the subsystems estimate data in orgdir/Subsystems.

=cut

sub subsystems_for_pegs_complete_estimate
{
    my($self, $pegs) = @_;

    #
    # Partition pegs on organisms.
    #
    my %peg_genomes;
    $peg_genomes{&FIG::genome_of($_)}->{$_} = 1 for @$pegs;

    # we want the behavior of capturing everything that is really in a subsystem or
    # estimated to be in a subsystem (RAO)
    my %results = $self->subsystems_for_pegs_complete($pegs);

    for my $genome (keys %peg_genomes)
    {
	my $pegs = $peg_genomes{$genome};

	my $dir = $self->organism_directory($genome);
	my $fh;
	if (!open($fh, "<", "$dir/Subsystems/bindings"))
	{
	    warn "subsystems_for_pegs_estimate: No bindings file found in $dir\n";
	    next;
	}

	my $sfh;
	if (!open($sfh, "<", "$dir/Subsystems/subsystems"))
	{
	    warn "subsystems_for_pegs_estimate: No subsystems file found\n";
	    return ();
	}
	#
	# Read variant codes.
	#
	my %variant;
	while (<$sfh>)
	{
	    chomp;
	    my($ss, $var) = split(/\t/);
	    $ss =~ s/\s+/_/g;
	    $variant{$ss} = $var;
	}
	close($sfh);

	while (<$fh>)
	{
	    chomp;
	    my($ss, $role, $peg) = split(/\t/);
	    next unless $pegs->{$peg};
	    $ss =~ s/\s+/_/g;
	    if (! &already_in($results{$peg},$ss))
	    {
		push(@{$results{$peg}}, [$ss, $role, $variant{$ss}]);
	    }
	}
	close($fh);
    }
    return (%results);
}

sub already_in {
    my($xL,$ss) = @_;

    my $i;
    if (! $xL) { return 0 }
    for ($i=0; ($i < @$xL) && ($xL->[$i] ne $ss); $i++) {}
    return ($i < @$xL);
}

=head3 subsystems_for_peg

 Return the list of subsystems and roles that this peg appears in.
 Returns an array. Each item in the array is
 a reference to a tuple of subsystem and role.  If the last argument ($noaux)
 is "true", only roles playing non-auxiliary roles will be returned.

=cut

sub subsystems_for_pegs
{
    my($self, $pegs,$noaux) = @_;
    my $rdbH = $self->db_handle;
    my %results;

    foreach my $peg (@$pegs){
	if ($self->is_deleted_fid($peg)) { next; }

	($peg =~ /^fig\|\d+\.\d+\.\w+\.\d+$/) or next;

	my $q = "SELECT subsystem, role FROM subsystem_index WHERE protein = '$peg'";

	if (my $relational_db_response = $rdbH->SQL($q))
	{
	    my %seen;
	    my @in;
	    my $pair;
	    foreach $pair (@$relational_db_response)
	    {
		$pair->[0] =~ s/ /_/g;
		my $key = join("\t",@$pair);
		if (! $seen{$key})
		{
		    push(@in,$pair);
		    $seen{$key} = 1;
		}
	    }
	    if ($noaux)
	    {
		my @nonaux = ();
		foreach my $x (@in)
		{
		    if (! $self->is_aux_role_in_subsystem($x->[0],$x->[1]))
		    {
			push(@nonaux,$x);
		    }
		}
		push (@{$results{$peg}}, @nonaux);
	    }
	    else
	    {
		push (@{$results{$peg}}, @in);
	    }
	}
	else
	{
	    push (@{$results{$peg}}, ());
	}
    }
    return (%results);
}

=head3 subsystems_roles

Return the list of subsystems and roles for every peg in subsystems
Returns an array. Each item in the array is
a reference to a three-ple of subsystem, role, and peg.

=cut

sub subsystems_roles
{
    my($self) = @_;

    my $rdbH = $self->db_handle;

    my $q = "SELECT subsystem, role, protein FROM subsystem_index";

    if (my $relational_db_response = $rdbH->SQL($q))
    {
        my %seen;
        my @in;
        my $pair;
        foreach $pair (@$relational_db_response)
        {
            my $key = join("\t",@$pair);
            if (! $seen{$key})
            {
                push(@in,$pair);
                $seen{$key} = 1;
            }
        }
        return @in;
    }
    else
    {
        return ();
    }
}


=head3 subsystems_for_role

Return a list of subsystems, roles, and proteins containing a given role

Returns an array. Each item in the array is a reference to a three-ple of subsystem, role, and peg.

=cut

sub subsystems_for_role
{
    my($self, $role) = @_;

    my $rdbH = $self->db_handle;

    my $roleQ = quotemeta $role;
    my $q = "SELECT subsystem, role, protein FROM subsystem_index WHERE role = \'$roleQ\'";

    if (my $relational_db_response = $rdbH->SQL($q))
    {
        my %seen;
        my @in;
        my $pair;
        foreach $pair (@$relational_db_response)
        {
            my $key = join("\t",@$pair);
            if (! $seen{$key})
            {
                push(@in,$pair);
                $seen{$key} = 1;
            }
        }
        return @in;
    }
    else
    {
        return ();
    }
}

=head3 subsystems_for_ec

Return a list of subsystems, roles, and proteins containing an EC number.

Returns an arrray. Each item in the array is a reference to a three-ple of subsystem, role, and peg.

=cut

sub subsystems_for_ec
{
    my($self, $ec) = @_;

    my $rdbH = $self->db_handle;

    my $q = "SELECT DISTINCT subsystem, role, protein FROM subsystem_index WHERE role like \'\%$ec\%\'";

    my $relational_db_response;
    if (($relational_db_response = $rdbH->SQL($q)) &&
            (@$relational_db_response > 0))
    {
        return @$relational_db_response;
    }
    else
    {
        return ();
    }
}



=head3 assigned_pegs_in_subsystems

Return list of [peg, function, ss, role in ss].

=cut

sub assigned_pegs_in_subsystems
{
    my($self, $genome) = @_;

    my @result = ();
    for my $peg ($self->pegs_of($genome))
    {
        my $fn = $self->function_of($peg);
        next if $fn eq "";
        next if $self->hypo($fn);

        my $rdbH = $self->db_handle;

        my $q = "SELECT subsystem, role FROM subsystem_index WHERE protein = '$peg'";

        if (my $relational_db_response = $rdbH->SQL($q))
        {
            my $pair;

            foreach $pair (@$relational_db_response)
            {
                  my ($ss, $role) = @$pair;

                  push(@result, [$peg, $fn, $ss, $role]);
             }

         }
    }
    return @result;
}


sub role_to_pegs {
    my($self,$role) = @_;

    my $rdbH = $self->db_handle;
    # $role =~ s/\'/\\\'/g;
    # my $q    = "SELECT protein FROM subsystem_index WHERE role = '$role'";
    # if (my $relational_db_response = $rdbH->SQL($q))
    # {
    #     return map { $_->[0] } @$relational_db_response;
    # }
    # return ();

    map { $_->[0] } @{ $rdbH->SQL( "SELECT prot FROM roles WHERE role = ?",
                                   undef, $role
                                 ) || []
                     };
}


sub peg_to_roles_in_subsystems {
    my($self,$peg) = @_;

    my $rdbH = $self->db_handle;
    my $q    = "SELECT subsystem, role FROM subsystem_index WHERE protein = '$peg'";
    if (my $relational_db_response = $rdbH->SQL($q))
    {
        return @$relational_db_response;
    }
    return ();
}

=head3 assigned_pegs_not_in_ss

Return all pegs with non-hypothetical assignments that are not in ss.

=cut

sub assigned_pegs_not_in_ss
{
    my($self, $genome) = @_;

    my @result = ();
    for my $peg ($self->pegs_of($genome))
    {
        my $fn = $self->function_of($peg);
        next if $fn eq "";
        next if $self->hypo($fn);

        my @subs = $self->subsystems_for_peg($peg);
        if (@subs < 1)
        {
            push(@result, [$peg, $fn, "No Subsytem", "No Role"]);
        }
    }
    return @result;
}

=head3 assigned_pegs

Return list of [peg, function, ss, role in ss] for every non-hypo protein regardless of being in ss

=cut

sub assigned_pegs
{
    my($self, $genome) = @_;

    my @result = ();
    for my $peg ($self->pegs_of($genome))
    {
        my $fn = $self->function_of($peg);
        next if $fn eq "";
        next if $self->hypo($fn);

        my $rdbH = $self->db_handle;

        my $q = "SELECT subsystem, role FROM subsystem_index WHERE protein = '$peg'";

        if (my $relational_db_response = $rdbH->SQL($q))
        {
            my $pair;

            if(@$relational_db_response > 0)
            {

               foreach $pair (@$relational_db_response)
               {
                  my ($ss, $role) = @$pair;

                  push(@result, [$peg, $fn, $ss, $role]);
               }
            }

            else
            {
             push(@result, [$peg, $fn, "No Subsystem", "No Role"]);
            }

        }
    }
    return @result;
}

sub ok_to_auto_update_subsys {
    my($self,$subsystem, $alter) = @_;

    # if alter > 0 we create the file. If alter < 0 we delete the file
    if (defined $alter && $alter > 0)
    {
    	open(OUT, ">$FIG_Config::data/Subsystems/$subsystem/ok.to.auto.update")
			|| die "We can't open the file $FIG_Config::data/Subsystems/$subsystem/ok.to.auto.update";
	print OUT "$subsystem\n";
	close OUT;
    }
    elsif (defined $alter && $alter < 0)
    {
    	unlink "$FIG_Config::data/Subsystems/$subsystem/ok.to.auto.update";
    }

    return -e "$FIG_Config::data/Subsystems/$subsystem/ok.to.auto.update";
}

=head3 subsystem_roles

Return a list of all roles present in locally-installed subsystems.
The return is a hash keyed on role name with each value a list
of subsystem names.

=cut

sub subsystem_roles
{
    my($self) = @_;

    my $rdbH = $self->db_handle;

    my $q = "SELECT distinct subsystem, role FROM subsystem_index";

    my $ret = {};

    if (my $relational_db_response = $rdbH->SQL($q))
    {
        foreach my $pair (@$relational_db_response)
        {
            my($subname, $role) = @$pair;
            push(@{$ret->{$role}}, $subname);
        }
    }

    return $ret;
}

#
# Return just the list of subsystems the peg appears in.
#
#      @subs = $fig->peg_to_subsystems($peg,"no-aux") will give only subsystems
#
# in which the PEG connect to a role that is not marked as "AUX"

sub peg_to_subsystems
{
    my($self, $peg, $noaux, $active_only) = @_;

    if ($self->is_deleted_fid($peg)) { return () }
    my @subs;
    my %in = map { $_->[0] =~ s/ /_/g; $_->[0] => 1 } $self->subsystems_for_peg($peg,$noaux, $active_only);
    return sort keys(%in);
}

sub pegs_to_usable_subsystems
{
    my($fig, $peg_list) = @_;

    my $dbh = $fig->db_handle->{_dbh};
    my $peg_str = join(", ", map { $dbh->quote($_) } @$peg_list);

    my $res = $dbh->selectall_arrayref(qq(SELECT DISTINCT si.protein, si.subsystem
			       FROM subsystem_index si
			           JOIN subsystem_metadata m ON si.subsystem = m.subsystem
			           LEFT JOIN deleted_fids df ON si.protein = df.fid
			       WHERE si.protein IN ($peg_str) AND
			             df.fid IS NULL AND
			             m.class_1 <> '' AND
			             m.class_1 NOT LIKE 'experimental%' COLLATE latin1_swedish_ci AND
			             m.class_1 NOT LIKE '%delete%' COLLATE latin1_swedish_ci));

    my $out = {};
    map { push @{$out->{$_->[0]}}, $_->[1] } @$res;
    return $out;
}

sub write_subsystem_spreadsheet {
    my($self,$ssa,$roles,$genomes,$pegs_in_cells) = @_;
    my(@genomes,$genome,$role,@pegs,$pair,$gs);

    $ssa =~ s/[ \/]/_/g;
    &verify_dir("$FIG_Config::data/Subsystems/$ssa");
    open(SSA,">$FIG_Config::data/Subsystems/$ssa/spreadsheet") || die "Cannot open $FIG_Config::data/Subsystems/$ssa/spreadsheet";
    foreach $pair (@$roles)
    {
        print SSA join("\t",@$pair),"\n";
    }
    print SSA "//\n";
    print SSA "All\n\nAll\n//\n";
    @genomes = map { $_->[1] }
               sort { ($a->[0] cmp $b->[0]) or ($a->[1] <=> $b->[1]) }
               map {$genome = $_; $gs = $self->genus_species($genome); [$gs,$genome] }
               @$genomes;
    foreach $genome (@genomes)
    {
        print SSA "$genome\t0";
        foreach $role (@$roles)
        {
            $_ = $pegs_in_cells->{"$genome\t$role->[1]"};
            @pegs = $_ ? sort { &by_fig_id($a,$b) } @{$_} : ();
            print SSA "\t",join(",",map { $_ =~ /^fig\|\d+\.\d+\.peg\.(\d+)/; $1 } @pegs);
        }
        print SSA "\n";
    }
    close(SSA);
    chmod(0777,"$FIG_Config::data/Subsystems/$ssa");
}

=head3 get_genome_subsystem_count

    my $num_subsytems = $fig->get_genome_subsystem_count($genomeID);

Return the number of subsystems of the genome identified by $genomeID.

=over 4

=item genomeID

ID of the genome whose number of subsystems is to be returned.

=item RETURN

Returns the number of subsystems.

=back

=cut

sub get_genome_subsystem_count {
  # Get the parameters.
  my ($self, $genomeID) = @_;
  # Declare the return variable.
  my $retVal;
  # Get the database handle.
  my $rdbH = $self->db_handle;

  my $dbh = $rdbH->{_dbh};

  $retVal = $rdbH->SQL(qq(SELECT COUNT(DISTINCT subsystem)
                                            FROM subsystem_index
                                            WHERE (protein like 'fig\|$genomeID.peg.%' AND
                                                   variant != '-1' AND variant ne '*-1')
                                           ));
  return $retVal->[0]->[0];
}

=head3 get_all_subsystem_pegs

    my @pegData = $fig->get_all_subsystem_pegs($genomeID);

Return the subsystems, roles, and variant codes for all features in the
specified genome. Unlike L</get_genome_subsystem_data>, this method
returns all pegs, regardless of the variant code.

=over 4

=item genomeID

ID of the relevant genome.

=item RETURN

Returns a hash that maps each subsystem ID to a list of 3-tuples, each
consisting of a role ID, a peg ID, and a variant code.

=back

=cut

sub get_all_subsystem_pegs {
    # Get the parameters.
    my ($self, $genomeID) = @_;
    # Get the database handle.
    my $rdbH = $self->db_handle;
    # Declare the return variable.
    my $qry;

    if ($FIG_Config::exclude_experimental_subsystems)
    {
	$qry = qq(SELECT i.subsystem, i.role, i.protein, i.variant
		  FROM subsystem_index i JOIN subsystem_metadata m ON i.subsystem = m.subsystem
		  WHERE protein LIKE ? ORDER BY subsystem, role
		  AND m.class_1 NOT LIKE 'Experimental%');
    }
    else
    {
	$qry = "SELECT subsystem, role, protein, variant " .
                             "FROM subsystem_index WHERE protein LIKE ? ORDER BY subsystem, role";
    }
    
    my $subList = $rdbH->SQL($qry, undef, "fig|$genomeID.%");
    # Form it into a hash of subsystems.
    my %retVal;
    for my $tuple (@$subList) {
        # Yank out the subsystem ID.
        my $subsysID = shift @$tuple;
        # Push the remainder of the tuple into the list for the specified subsystem.
        push @{$retVal{$subsysID}}, $tuple;
    }
    # Return the result.
    return %retVal;
}


=head3 get_genome_subsystem_data

    my $roleList = $fig->get_genome_subsystem_data($genomeID);

Return the roles and pegs for a genome's participation in subsystems. The
subsystem name, role ID, and feature ID will be returned for each of
the genome's subsystem-related PEGs.

=over 4

=item genomeID

ID of the genome whose PEG breakdown is desired.

=item RETURN

Returns a pointer to a list of 4-tuples. Each tuple consists of a subsystem name, a role ID,
a feature ID, and the variant code.

=back

=cut

sub get_genome_subsystem_data {
    # Get the parameters.
    my ($self, $genomeID,$all) = @_;

    my %out;
    if ($FIG_Config::use_subsystem_estimates)
    {
	my $active = $self->active_subsystems_estimate($genomeID, $all);
	if (open(my $bfh, "<", $self->organism_directory($genomeID) . "/Subsystems/bindings"))
	{
	    while (<$bfh>)
	    {
		chomp;
		my($ss, $role, $fid) = split(/\t/);
		if ($active->{$ss} && !$self->is_deleted_fid($fid))
		{
		    # push(@$out, [$ss, $role, $fid, $active->{$ss}]);
		    $out{$ss, $role, $fid} = $active->{$ss};
		}
	    }
	    close($bfh);
	}
    }

    # Declare the return variable.
    my $retVal;
    # Get the database handle.
    my $rdbH = $self->db_handle;

    #
    # For now need to try with variant first, then back off to not using variant
    # if we hit a database error.
    #

    {
        my $dbh = $rdbH->{_dbh};
        local $dbh->{RaiseError} = 1;
        local $dbh->{PrintError} = 0;
	my $constraint = "WHERE (protein like 'fig\|$genomeID.peg.%')";
	if (! $all)
	{
	    $constraint = $constraint . " AND (variant != '-1') AND (variant != '0')";
	    $constraint = $constraint . " AND (variant != '*-1') AND (variant != '*0')";
	}

        eval {
            $retVal = $rdbH->SQL(qq(SELECT DISTINCT subsystem,role,protein,variant
                                            FROM subsystem_index
                                            $constraint
                                           ));

        };
    }

    if ($@ =~ /variant/)
    {
        $retVal = $rdbH->SQL(qq(SELECT DISTINCT subsystem,role,protein,variant
                                        FROM subsystem_index
                                        WHERE (protein like 'fig\|$genomeID.peg.%')
                                       ));
    }

    for my $ent (@$retVal)
    {
	my($ss, $role, $fid, $var) = @$ent;
	$out{$ss, $role, $fid} = $var;
    }


    #
    # We've merged the estimates and real in %out; now reformat into the real output list.
    #

    @$retVal = ();

    for my $k (sort keys %out)
    {
	my($ss, $role, $fid) = split(/$;/, $k);
	push(@$retVal, [$ss, $role, $fid, $out{$k}]);
    }
     
    # Return the result.
    return $retVal;
}

=head3 mark_subsystems_modified

    $fig->mark_subsystems_modified()

Update the timestamp on the subsystem-modified flag stored in FIG/var. This will
trigger a subsequent subsystem information call that uses the FIG.pm subsystem
data cache to reload the cache from the database.

=cut

sub mark_subsystems_modified
{
    my($self) = @_;
    my $f = $FIG_Config::var . "/ss_modified";
    unlink($f);
    open(F, ">", $f);
    print F time . "\n";
    close(F);
}

=head3 get_genome_stats

    my ($gname,$szdna,$pegs,$rnas,$taxonomy) = $fig->get_genome_stats($genomeID);

Return basic statistics about a genome.

=over 4

=item genomeID

ID of the relevant genome.

=item RETURN

Returns a 5-tuple containing the genome name, number of base pairs, number of PEG
features, number of RNA features, and the taxonomy string.

=back

=cut

sub get_genome_stats {
    # Get the parameters.
    my ($self, $genomeID) = @_;

    if (!$self->is_genome($genomeID))
    {
	return ();
    }

    my $rdbH = $self->db_handle;
    my $relational_db_response = $rdbH->SQL("SELECT gname,szdna,pegs,rnas,taxonomy FROM genome WHERE genome = '$genomeID'");
    my($db_gname, $db_szdna, undef, undef, $tax) = @{$relational_db_response->[0]};
    #
    # Need to patch with the actual counts of the PEGs and RNAs.
    #

    my $rna_count = $self->all_features($genomeID, 'rna');
    my $peg_count = $self->all_features($genomeID, 'peg');

    return ($db_gname, $db_szdna, $peg_count, $rna_count, $tax);
}

=head3 get_genome_assignment_data

    my $roleList = $fig->get_genome_subsystem_data($genomeID);

Return the functional assignments and pegs for a genome. The feature ID and assigned
function will be returned for each of the genome's PEGs.

=over 4

=item genomeID

ID of the genome whose PEG breakdown is desired.

=item RETURN

Returns a list of 2-tuples. Each tuple consists of a peg ID and its master
functional assignment.

=back

=cut

sub get_genome_assignment_data {
    # Get the parameters.
    my ($self, $genomeID) = @_;

    my @p = $self->pegs_of($genomeID);
    my $f = $self->function_of_bulk(\@p);

    return [ grep { not $self->is_deleted_fid($_->[0]) } map { [ $_, $f->{$_} ] } keys %$f];


    # Get the database handle.
    my $rdbH = $self->db_handle;
    # Get the data.
    my $retVal = $rdbH->SQL("SELECT prot,assigned_function FROM assigned_functions WHERE ( prot like 'fig\|$genomeID.peg.%' AND made_by = 'master' )");
    # Return it.

    my @filtered = grep  { not $self->is_deleted_fid($_->[0]) } @$retVal;

    return [@filtered];
}

sub get_basic_statistics
{
    my($self, $genome) = @_;

    #
    # Check cache.
    #

    my $fh = $self->get_valid_cache_file("$genome/basic_statistics");
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


    my $fh_cache = $self->write_cache_file("$genome/basic_statistics");
    if ($fh_cache)

    {
	while (my($k, $v) = each %$statistics)
	{
	    print $fh_cache join("\t", $k, $v), "\n";
	}
	close($fh_cache);
    }

    return $statistics;
}


sub get_peg_statistics {
    my ($self, $genome) = @_;

    #
    # Check cache.
    #

    my $fh = $self->get_valid_cache_file("$genome/peg_statistics");
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
    my $assignment_data = $self->get_genome_assignment_data($genome);

    my $hypo_sub = 0;
    my $hypo_nosub = 0;
    my $nothypo_sub = 0;
    my $nothypo_nosub = 0;
    my %in = map { $_->[2] => 1 } @$subsystem_data;
    my $in = keys(%in);

    my %sscount = map { $_->[0] => 1 } @$subsystem_data;

    foreach $_ (@$assignment_data)
    {
	my($peg,$func) = @$_;
	my $is_hypo = &FIG::hypo($func);

	if    ($is_hypo && $in{$peg})           { $hypo_sub++ }
	elsif ($is_hypo && ! $in{$peg})         { $hypo_nosub++ }
	elsif ((! $is_hypo) && (! $in{$peg}))   { $nothypo_nosub++ }
	elsif ((! $is_hypo) && $in{$peg})       { $nothypo_sub++ }
    }
    my $tot = $hypo_sub + $nothypo_sub + $hypo_nosub + $nothypo_nosub;

    my ($fracHS, $fracNHS, $fracHNS, $fracNHNS);

    if ($tot == 0) {
	$fracHS = sprintf "%.2f", 0.0;
	$fracNHS = sprintf "%.2f", 0.0;
	$fracHNS = sprintf "%.2f", 0.0;
	$fracNHNS = sprintf "%.2f", 0.0;
    } else {
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


    my $fh_cache = $self->write_cache_file("/$genome/peg_statistics");
    if ($fh_cache)

    {
	while (my($k, $v) = each %$statistics)
	{
	    print $fh_cache join("\t", $k, $v), "\n";
	}
	close($fh_cache);
    }

    return $statistics;
}

################################ Caching #################################
#
# Code for supporting caching of  commonly-used data.
#
# Caches are updated by the update-caches script, and should be done
# on mirrored machines each time mirroring is done, and on non-mirrored
# machines at some interval appropriate to the load on the machine.
#

=head3 get_valid_cache_file

If the given cache file (name is relative to the FIG cache directory) exists and
is less than a day old (Parameterize this sometime!) open and return a filehandle.

=cut

sub get_valid_cache_file
{
    my($self, $file) = @_;

    my $dir = $self->get_cache_directory();
    my $path = "$dir/$file";

    my $fh = new FileHandle($path);
    if ($fh)
    {
	my @s = stat($fh);
	my $age = time - $s[9];
	if ($age > 86400)
	{
	    $fh->close();
	    return undef;
	}
	else
	{
	    return $fh;
	}
    }
    else
    {
	return undef;
    }
}

sub write_cache_file
{
    my($self, $file) = @_;

    my $dir = $self->get_cache_directory();
    my $path = "$dir/$file";
    my $sdir = dirname($path);
    &FIG::verify_dir($sdir);

    my $fh = new FileHandle(">$path");
    return $fh;
}

sub get_cache_directory {
    my($self) = @_;
    my $cache_dir;

    my $cinfo = $self->cached('_cache_info');
    if (defined($cache_dir = $cinfo->{directory})) {
	return $cache_dir;
    }
    else {
	$cache_dir = $FIG_Config::cache_dir;
    }

    if (!$cache_dir)
    {
	$cache_dir = "$FIG_Config::var/seed_cache";
    }

    &FIG::verify_dir($cache_dir);

    $cinfo->{directory} = $cache_dir;
    return $cache_dir;
}





################################# Dlit Stuff       ###################################

=head3 add_dlit

    $rc = $fig->add_dlit(
		          -status   => 'D',       # required
		          -peg      => $peg,      # or -md5 => $md5,  # one is required
		          -pubmed   => $pubmed,   # required
		          -curator  => 'RossO',   # required
		          -go       => '',        # default = ''
		          -override => 1);        # default = 0

This adds a dlit tuple.  The currently supported arguments are

    -status =>          ' '  for not curated
                        'D'  for dlit (direct literature on role)
                        'G'  for genome data (propagates to all ' ' entries for this article)
                        'N'  for not relevant
                        'R'  for relevant, but not dlit

    -md5    =>          supply an md5 hash code for the peg, not the id.

    -peg    =>          the peg being connected to literature.  This peg will
                        be treated as a representative of the set that have the
                        same protein sequence.

    -pubmed =>          pubmed ID (all numeric, but stored as string)

    -curator =>         curator making the assertion (30 char max)

    -go     =>          an optional list of 3-character codes separated by commas

    -override =>        0 -> if there is an existing tuple, ignore this request
                        1 -> if there is an existing tuple, replace it

The returned value will be

                        0 -> the tuple was not inserted
                        1 -> the tuple was inserted
=cut

sub add_dlit {
    return Dlits::add_dlit(@_);
}


=head3 all_dlits

    $dlits = $fig->all_dlits();

Returns a reference to an array of all current dlit data.

The returned value is

    [ [ status, md5_hash, pubmed, curator, go_code ], ... ]
=cut

sub all_dlits {
    return Dlits::all_dlits(@_);
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
    return Dlits::add_title(@_);
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
    return Dlits::update_title(@_);
}


=head3 get_title

    $title = $fig->get_title( $pubmed_id )

Get a title for a literature id

Returned value:

    $title   upon success
    undef    upon failure

=cut

sub get_title {
    return Dlits::get_title(@_);}


=head3 all_titles

    [ [ id, title ], ... ] = $fig->all_titles()

Get all pubmed_id, title pairs

Returned value:

    [ [ id, title ], ... ]   upon success
    []                       upon failure

=cut

sub all_titles {
    return Dlits::all_titles(@_);
}



=head3 get_dlits_for_peg

    $dlits = $fig->get_dlits_for_peg();

Returns a reference to an array of current dlit data for a peg.

The returned value is

    [ [ status, md5_hash, pubmed, curator, go_code ], ... ]
=cut

sub get_dlits_for_peg {
    return Dlits::get_dlits_for_peg(@_);
}

=head3 get_dlits_for_pegs

    $dlits = $fig->get_dlits_for_pegs(); 

Returns a reference to an array of current dlit data for a list of pegs.

The returned value is

    [ [ status, md5_hash, pubmed, curator ], ... ]
=cut

sub get_dlits_for_pegs {
    return Dlits::get_dlits_for_pegs(@_);
}

################################# PEG Translation  ####################################

=head2 PEG Translations

=cut

sub translate_pegs {
    my($self,$pegs,$seq_of, $cb) = @_;
    my($seq,$aliases,$pegT,%to,%sought,@keys,$peg,$alias);

    $cb = sub {} unless ref($cb) eq "CODE";

    my $tran_peg = {};

    my $n = scalar keys (%$pegs);
    my $idx = 0;

    foreach $peg (keys(%$pegs))
    {
        $idx++;
        &$cb("$idx of $n") if $idx % 100 == 0;

	if ($self->ftype($peg) eq 'peg')
	{
	    #
	    # First, see if the peg of the same name locally has the same
	    # last 10 chars.
	    #
	    if (($seq = $self->get_translation($peg)) &&
		(length($seq) > 10) && (length($seq_of->{$peg}) > 10) &&
		(uc substr($seq,-10) eq substr($seq_of->{$peg},-10)))
	    {
		$tran_peg->{$peg} = $peg;
	    }
	    else
	    {
		#
		# Otherwise, search for a local peg that has the same alias
		# as this peg. (Canonicalize based on the original source)
		#
		($aliases,undef,undef) = @{$pegs->{$peg}};
		undef %to;
		foreach $alias (split(/,/,$aliases))
		{
		    if ($pegT = $self->by_alias($alias))
		    {
			$to{$pegT}++;
		    }
		}

		#
		# If we have a unique answer, we are done.
		# Otherwise mark this one as needing more search.
		#
		if ((@keys = keys(%to)) == 1)
		{
		    $tran_peg->{$peg} = $keys[0];
		}
		else
		{
		    $sought{$peg} = 1;
		}
	    }
	}
	else
	{
	    #
	    # This is some other sort of feature.
	    #
	    # Just check to see that the local feature of the same name has the same DNA
	    # sequence. If not, we don't match, and we probably can't.
	    #

	    my $local_dna = $self->dna_seq($peg);
	    if ($local_dna eq $seq_of->{$peg})
	    {
		$tran_peg->{$peg} = $peg;
	    }
	    else
	    {
		warn "no local match for $peg $seq_of->{$peg}\n";
	    }
        }
    }

    if ((scalar keys(%sought)) > 0)
    {
        &tough_search($self,$pegs,$seq_of,$tran_peg,\%sought);
    }
    return $tran_peg;
}

=head3 tough_search($pegs, $seq_of, $tran_peg, $sought)

$pegs - not used
$seq_of - hash from peg to peg sequence
$tran_peg - hash into which translated pegs are placed
$sought - hash keyed on the list of pegs we're looking for.

=cut

sub tough_search {
    my($self,$pegs,$seq_of,$tran_peg,$sought) = @_;
    my($peg,$seq,%needH,%needT,%poss,$id,$sub,$to,$x,$genome);

    #
    # Construct %needT, key is 50-bases from tail of sequence, values are pegs from
    # the list of pegs we're seeking.
    #
    # %needH is the same, but keyed on 50 bases from the head of the sequence.
    #
    foreach $peg (keys(%$sought))
    {
        if (($seq = $seq_of->{$peg}) && (length($seq) > 50))
        {
            my $sub = substr($seq,-50);
            push(@{$needT{$sub}},$peg);
            $sub = substr($seq,0,50);
            push(@{$needH{$sub}},$peg);
        }
    }
#   print STDERR &Dumper(\%needT,\%needH);
    open(NR,"<$FIG_Config::global/nr")
        || die "could not open $FIG_Config::global/nr";
    $/ = "\n>";
    while (defined($_ = <NR>))
    {
        chomp;
        if ($_ =~ /^>?(\S+)[^\n]*\n(.*)/s)
        {
            $id  =  $1;
            $seq =  $2;
            if (length($seq) >= 50)
            {
                $sub = uc substr($seq,-50);
                if ((($x = $needT{uc substr($seq,-50)}) && (@$x == 1)) ||
                    (($x = $needH{uc substr($seq,0,50)}) && (@$x == 1)))
                {
                    $peg = $x->[0];
                    my @same = grep { $_ =~ /^fig/ } map { $_->[0] } $self->mapped_prot_ids($id);
                    if (@same > 0)
                    {
                        push(@{$poss{$peg}},@same);
                    }
                }
            }
        }
    }
    close(NR);
    $/ = "\n";

#   print STDERR &Dumper(\%poss);
    foreach $peg (keys(%poss))
    {
#       print STDERR "processing $peg\n";
        $x = $poss{$peg};
        if (@$x == 1)
        {
            $tran_peg->{$peg} = $x->[0];
            delete $sought->{$peg};
        }
        elsif ($genome = $self->probable_genome($self->genome_of($peg),$tran_peg))
        {
#           print STDERR "    mapped to genome $genome\n";
            my $genomeQ = quotemeta $genome;
            my @tmp = grep { $_ =~ /^fig\|$genomeQ\./ } @$x;
            if (@tmp == 1)
            {
                $tran_peg->{$peg} = $tmp[0];
                delete $sought->{$peg};
            }
            else
            {
#               print STDERR &Dumper(["failed",$peg,$genome,\@tmp]);
            }
        }
        else
        {
#           print STDERR "could not place genome for $peg\n";
        }
    }

    foreach $peg (keys(%$sought))
    {
        print STDERR "failed to map $peg\n";
    }
}

sub probable_genome {
    my($self,$genome,$tran_peg) = @_;
    my($peg,%poss,@poss);

    my $genomeQ = quotemeta $genome;
    foreach $peg (grep { $_ =~ /^fig\|$genomeQ\./ } keys(%$tran_peg))
    {
        $poss{$self->genome_of($tran_peg->{$peg})}++;
    }
    @poss = sort { $poss{$b} <=> $poss{$a} } keys(%poss);
    if ((@poss == 1) || ((@poss > 1) && ($poss{$poss[0]} > $poss{$poss[1]})))
    {
        return $poss[0];
    }
    else
    {
#       print STDERR &Dumper(["could not resolve",\%poss,$genome]);
        return undef;
    }
}

=head3 find_genome_by_content

Find a genome given the number of contigs, number of nucleotides, and
checksum. We pass in a potential name for the genome as a quick
starting check.

=cut

sub find_genome_by_content
{
    my($self, $genome, $n_contigs, $n_nucs, $cksum) = @_;

    my(@genomes);

    my $gbase = (split(/\./, $genome))[0];

    #
    # Construct the list of genomes so that we first check ones with the same
    # base-part as the one passed in.
    #
    for ($self->genomes())
    {
        if (/^$gbase\./)
        {
            unshift(@genomes, $_);
        }
        else
        {
            push(@genomes, $_);
        }
    }


    for my $genome (@genomes)
    {
        if (open(my $cfh, "<$FIG_Config::organisms/$genome/COUNTS"))
        {
            if (defined($_ = <$cfh>))
            {
                my($cgenome, $cn_contigs, $cn_nucs, $ccksum) = split(/\t/);

                if ($cgenome eq $genome and
                    $cn_contigs == $n_contigs and
                    $cn_nucs == $n_nucs and
                    $ccksum == $cksum)
                {
                    return $genome;
                }
            }
            close($cfh);
        }
    }

    return undef;
}

################################# Support for PEG Links  ####################################

=head2 Links

=cut

sub peg_links {
    my($self,$fid) = @_;

    return $self->fid_links($fid);
}

=head3 fid_links

    my @links = $fig->fid_links($fid);

Return a list of hyperlinks to web resources about a specified feature.

=over 4

=item fid

ID of the feature whose hyperlinks are desired.

=item RETURN

Returns a list of raw HTML strings representing hyperlinks to web pages relating to
the specified feature.

=back

=cut
#: Return Type @;
sub fid_links {
    my($self,$fid,$aliases) = @_;
    my($i,$got,$genome,$fidN);

    if ($self->is_deleted_fid($fid)) { return () }
    my @links = ();

    my @aliases = ref($aliases) ? @$aliases : $self->feature_aliases($fid);

    my $link_file;
    for $link_file (("$FIG_Config::global/fid.links","$FIG_Config::global/peg.links"))
    {
        if (open(GLOBAL,"<$link_file"))
        {
            while (defined($_ = <GLOBAL>))
            {
                chop;
                my($pat,$link) = split(/\t/,$_);
                for ($i=0,$got=0; (! $got) && ($i < @aliases); $i++)
                {
                    if ($aliases[$i] =~ /$pat/)
                    {
                        push(@links,eval "\"$link\"");
                        $got = 1;
                    }
                }
            }
            close(GLOBAL);
        }
    }

    my $relational_db_response;
    my $rdbH = $self->db_handle;

    if (($relational_db_response = $rdbH->SQL("SELECT link FROM fid_links WHERE ( fid = \'$fid\' )")) &&
        (@$relational_db_response > 0))
    {
        push(@links, map { $_->[0] } @$relational_db_response);
    }
    return sort { $a =~ /\>([^\<]+)\<\/a\>/; my $l1 = $1;
                  $b =~ /\>([^\<]+)\<\/a\>/; my $l2 = $1;
                  $a cmp $b } @links;
}

=head3 fids_with_link_to

    my @links = $fig->fids_with_link_to("text");

Return a list of tples of [fid, link] where text is a free-text string that will match to the URL. You can use this to get all the links that point to PIR, for example to identify all proteins that are members of PIR superfamilies.

=over 4

=item text

A free-text match to the URL. The match is made using the SQL "like" command, so try to be as specific as possible.

=item RETURN

Returns a list tuples of [fid, link]

=back

=cut

sub fids_with_link_to {
    my($self,$text) = @_;

    my $relational_db_response;
    my $rdbH = $self->db_handle;

    if (($relational_db_response = $rdbH->SQL("SELECT fid,link FROM fid_links WHERE ( link like \'\%$text\%\' )")) &&
        (@$relational_db_response > 0))
    {
        return $relational_db_response;
    }
    return ();
}


sub add_peg_links {
    my($self,@links) = @_;

    return $self->add_fid_links(@links);
}

sub add_fid_links {
    my($self,@links) = @_;
    my($fid,$link,$pair,$i);

    my $relational_db_response;
    my $rdbH = $self->db_handle;

    foreach $pair (@links)
    {
        ($fid,$link) = @$pair;

        if (($relational_db_response = $rdbH->SQL("SELECT link FROM fid_links WHERE ( fid = \'$fid\' )")) &&
            (@$relational_db_response > 0))
        {
            for ($i=0; ($i < @$relational_db_response) && ($relational_db_response->[$i]->[0] ne $link); $i++) {}
            if ($i == @$relational_db_response)
            {
                &add_fid_link($self,$fid,$link);
            }
        }
        else
        {
            &add_fid_link($self,$fid,$link);
        }
    }
}

sub add_fid_link {
    my($self,$fid,$link) = @_;

    if ($self->is_deleted_fid($fid)) { return }
    my $rdbH = $self->db_handle;

    ($fid =~ /^fig\|(\d+\.\d+)\.([^.]+)\.\d+$/) || confess "bad fid $fid";
    my $genome = $1;
    my $type   = $2;

    $rdbH->SQL("INSERT INTO fid_links ( fid,link ) VALUES ( \'$fid\',\'$link\' )");
    if ( $genome && $type && open( TMP, ">>$FIG_Config::organisms/$genome/Features/$type/links" ) )
    {
        print TMP "$fid\t$link\n";
        close(TMP);
        chmod 0777,"$FIG_Config::organisms/$1/Features/$type/links";
    }
}

sub delete_peg_link {
    my($self,$peg,$link) = @_;

    $self->delete_fid_link($peg,$link);
}

sub delete_fid_link {
    my($self,$fid,$link) = @_;
    my($i);

    if ($self->is_deleted_fid($fid)) { return }
    my $genome = $self->genome_of($fid);

    ($fid =~ /^fig\|\d+\.\d+\.([^.]+)\.\d+$/) || confess "bad fid $fid";
    my $type = $1;

    my $rdbH = $self->db_handle;
    $rdbH->SQL("DELETE FROM fid_links WHERE ( fid = \'$fid\' AND link = \'$link\' )");

    my $file;
    foreach $file (("$FIG_Config::organisms/$genome/Features/$type/$type.links","$FIG_Config::organisms/$genome/Features/$type/links"))
    {
        if (-s $file)
        {
            my @links = `cat $file`;
            for ($i=0; ($i < @links) && (! (($links[$i] =~ /^(\S+)\t(\S.*\S)/) && ($1 eq $fid) && ($2 eq $link))); $i++) {}
            if (($i < @links) && open(TMP,">$file"))
            {
                splice(@links,$i,1);
                print TMP join("",@links);
                close(TMP);
            }
        }
    }
}

sub delete_all_peg_links {
    my($self,$peg) = @_;

    $self->delete_all_fid_links($peg);
}

sub delete_all_fid_links {
    my($self,$fid) = @_;
    my($i);

    if ($self->is_deleted_fid($fid)) { return }
    my $genome = $self->genome_of($fid);

    my $rdbH = $self->db_handle;
    $rdbH->SQL("DELETE FROM fid_links WHERE ( fid = \'$fid\' )");

    ($fid =~ /^fig\|\d+\.\d+\.([^.]+)\.\d+$/) || confess "bad fid $fid";
    my $type = $1;

    my $file;
    foreach $file (("$FIG_Config::organisms/$genome/Features/$type/$type.links","$FIG_Config::organisms/$genome/Features/$type/links"))
    {
        if (-s $file)
        {
            my @links = `cat $file`;
            my @links1 = grep { ! (($_ =~ /^(\S+)/) && ($1 eq $fid)) } @links;
            if ((@links1 < @links) && open(TMP,">$file"))
            {
                print TMP join("",@links1);
                close(TMP);
            }
        }
    }
}


=head3 Search Database

Searches the database for objects that match the query string in some way.

Returns a list of results if the query is ambiguous or an unique identifier
otherwise.

=cut

sub search_database {
  # get parameters
  my ($self, $query, $options) = @_;

  # get cgi
  my $cgi = new CGI;

  # turn query string into lower case
  $query = lc($query);
  my $ss_query = $query;
  $ss_query =~ s/ /_/g;
  my @tokenized = split(/ /, $query);

  # check for options, otherwise set default values
  my $limit = 100;
  if (defined($options)) {
      if ($options->{limit}) {
	  $limit = $options->{limit};
      }
  }

  # get database handle
  my $dbh = $self->db_handle();

  # check exact organism name and id
  my $result = $dbh->SQL("SELECT genome FROM genome WHERE LOWER(gname)='$query' OR genome='$query'");
  if (scalar(@$result) > 0) { return { type => 'organism', result => $result->[0]->[0] }; }

  # check exact subsystem
  $result = $dbh->SQL("SELECT subsystem FROM subsystem_index WHERE LOWER(subsystem)='$ss_query'");
  if (scalar(@$result) > 0) { return { type => 'subsystem', result => $result->[0]->[0] }; }

  # check fig-id
  $result = $dbh->SQL("SELECT id FROM features WHERE id='$query'");
  if (scalar(@$result) > 0) { return { type => 'feature', result => $result->[0]->[0] }; }

  # check unique alias
  $result = $dbh->SQL("SELECT id FROM ext_alias WHERE alias='$query'");
  if (scalar(@$result) > 0) { return { type => 'feature', result => $result->[0]->[0] }; }

  # exact search failed, sum up all the fuzzy searches
  my $return_value;

  # check functional role
  $result = $dbh->SQL("SELECT DISTINCT role, subsystem FROM subsystem_index WHERE LOWER(role) LIKE '%" . $query . "%'");
  if (scalar(@$result) > 0) { push(@$return_value, { type => 'functional_role', result => $result }); }

  # check organism name and domain
  $result = $dbh->SQL("SELECT DISTINCT genome, gname, maindomain FROM genome WHERE LOWER(gname) LIKE '%" . $query . "%' OR LOWER(maindomain)='$query'");
  if (scalar(@$result) > 0) { push(@$return_value, { type => 'organism', result => $result }); }

  # check aliases
  $result = $dbh->SQL("SELECT DISTINCT id, genome FROM features WHERE LOWER(aliases) like '%". $query . "%'");
  if (scalar(@$result)) { push(@$return_value, { type => 'feature', result => $result }); }

  # check subsystem
  $result = $dbh->SQL("SELECT DISTINCT subsystem FROM subsystem_index WHERE LOWER(subsystem) LIKE '%" . $ss_query . "%'");
  if (scalar(@$result) > 0) { push(@$return_value, { type => 'subsystem', result => $result }); }

  # check for extended search
  unless ($cgi->param('quick_search')) {
    my @tokens;
    foreach (@tokenized) {
      push(@tokens, "LOWER(role) LIKE '%" . $_ . "%'");
    }
    my $comp = join(' AND ', @tokens);
    $result = $dbh->SQL("SELECT DISTINCT prot, role, gname, org, maindomain FROM roles JOIN genome ON roles.org=genome.genome WHERE prot LIKE 'fig%' AND " . $comp . " LIMIT $limit");
    if (scalar(@$result) > 0) { 
	    # roles can contain a peg more than once, and we just need to get the last entry for each peg
	    my $flatten;
	    map {$flatten->{$_->[0]} = $_} grep {!$self->is_deleted_fid($_->[0])} @$result;
	    $result = [values %$flatten];
	    push(@$return_value, { type => 'proteins', result => $result }); 
    }
  }

  return $return_value;
}

sub flat {
  my ($in) = @_;

  my $out;

  foreach (@$in) { push(@$out, $_->[0]); }

  return $out;
}

###########
#
#

=head2 Peg Searches and Similarities

Some routines for dealing with peg search and similarities.

This is code lifted from pom.cgi and reformatted for more general use.

Find the given role in the given (via CGI params) organism.

We do this by finding a list of pegs that are annotated to have
this role in other organisms that are "close enough" to our organism

We then find pegs in this organism that are similar to
these pegs.

=cut

sub find_role_in_org
{
    my($self,$role, $org, $user, $sims_cutoff) = @_;

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

    $self->try_to_locate($org,$hits,[@cand[0..$how_many0]],$seen, $sims_cutoff);

    if (keys(%$hits) == 0)
    {
        splice(@cand,0,$how_many0+1);
        &try_to_locate($self,$org,$hits,\@cand,$seen, $sims_cutoff);
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

=head2 Utility Methods

=head3 is_ec

    my $flag = FIG::is_ec($role);

Return TRUE if the specified role is an EC number, else FALSE. This can be used to
determine whether a role is specified via a role ID or the role's EC number.

=over 4

=item role

Role ID or EC number to check.

=item RETURN

Returns TRUE if the specified role specification is an EC number, and FALSE if it
is a true role ID.

=back

=cut

sub is_ec {
    # Get the parameter.
    my ($role) = @_;
    # Check its structural format.
    my $retVal = ($role =~ /^(\d+|-)\.(\d+|-)\.(\d+|-)\.(\d+|-)$/);
    # Return the match indicator.
    return $retVal;

}

=head3 run_in_background

Background job support.

If one wants to turn a script into a background, invoke $fig->run_in_background($coderef). This will cause $coderef to be invoked as a background job. This means its output will be written to $FIG_Config::data/Global/background_jobs/<pid>, and that it shows up and is killable via the seed control panel.

=cut

sub run_in_background
{
    my($self, $coderef, $close_fds) = @_;

    if (ref($coderef) ne "CODE")
    {
        warn "FIG::run_in_background() invoked without a code reference\n";
        return;
    }

    my $job_dir = "$FIG_Config::data/Global/background_jobs";
    verify_dir($job_dir);

    my $child = fork;

    if (!defined($child))
    {
        die "run_in_background: fork failed: $!\n";
    }

    if ($child == 0)
    {
        #
        # In the child process.
        #

        POSIX::setsid();

        my $d = $self->db_handle();
        if ($d)
        {
            my $dbh = $d->{_dbh};
            $dbh->{InactiveDestroy} = 1;
        }

        if ($close_fds)
        {
            for (my $fd = 3; $fd < 32; $fd++)
            {
                POSIX::close($fd);
            }
        }

        my $my_job_dir = "$job_dir/$$";
        verify_dir($my_job_dir);

        open(my $info, ">$my_job_dir/INFO");
        my $now = asctime(localtime(time()));
        chomp($now);
        print $info "Background job $$ created from run_in_background by $> on $now\n";
        close($info);

        #
        # Redirect stdout/stderr to the OUTPUT file.
        #

        close(STDOUT);
        close(STDERR);

        open STDOUT, ">$my_job_dir/OUTPUT";
        open STDERR, ">&STDOUT";

        select(STDERR);
        $| = 1;
        select(STDOUT);
        $| = 1;

        #
        # Make stdin be from nowhere.
        #

        open STDIN, "</dev/null";

        #
        # Run the code.
        #

        open(my $stat0, ">$my_job_dir/STATUS");
        print $stat0 "Job started at $now\n";
        close($stat0);
        eval { &$coderef; };

        open(my $stat, ">$my_job_dir/STATUS");

        if ($@ eq "")
        {
            print $stat "Finished successfully\n";
        }
        else
        {
            print $stat "Code had an error:\n$@\n";
        }
        close($stat);

        #
        # We need to undef this, otherwise the DBrtns destructor
        # will do an explicit dbh->disconnect, which will undo any
        # effect of the InactiveDestroy set above.
        #

        $d = $self->db_handle();
        if ($d)
        {
            delete $d->{_dbh};
        }

        exit;
    }

    return $child;
}

sub get_all_jobs :List
{
    my($self) = @_;

    my $job_dir = "$FIG_Config::data/Global/background_jobs";

    opendir(my $jfh, $job_dir);
    my @jobs = grep { $_ =~ /^\d+$/ } readdir($jfh);
    closedir($jfh);

    return @jobs;
}


sub get_job :Scalar
{
    my($self, $job_id) = @_;

    my $job_dir = "$FIG_Config::data/Global/background_jobs/$job_id";

    if (-d $job_dir)
    {
        return FIG::Job->new($job_id, $job_dir);
    }
    else
    {
        return undef;
    }
}

sub get_current_arch :Scalar
{
    my $arch;

    open(FH, "<$FIG_Config::fig_disk/config/fig-user-env.sh");
    while (<FH>)
    {
        chomp;
        if (/^RTARCH=\"(.*)\"/)
        {
            $arch = $1;
            last;
        }
    }
    return $arch;
}


############################### Interfaces to Other Systems ######################################
#

=head2 External Interface Methods

This section contains the functionality introduced by the interface with GenDB.  The initial
two functions simply register when GenDB has a version of the genome (so we can set links
to it when displaying PEGS:

=head3 has_genome

usage: has_genome("GenDB",$genome)

Invoking this routine just records that GenDB has a copy of the genome designated by
$genome.

=cut

sub has_genome {
    my($system,$genome) = @_;

    print STDERR "$system has $genome\n";
}

=head3 dropped_genome

usage: dropped_genome("GenDB",$genome)

Invoking this routine just records that GenDB should no longer be viewed as
having a copy of the genome designated by $genome.

=cut

sub dropped_genome {
    my($system,$genome) = @_;

    print STDERR "$system dropped $genome\n";
}

=head3 link_to_system

usage: $url = link_to_system("GenDB",$fid) # usually $fid is a peg, but it can be other types of features, as well

This routine is used to get a URL that can be used to "flip" from one system to the other.  If
the feature is unknown to the system, undef should be returned.

=cut

sub link_to_system {
    my($system,$fid) = @_;

    return undef;
}

############################### Adding, Deleting, Altering Features  ####################

=head2 Feature Update Methods

The following routines support alteration of features

=head3 delete_feature

usage: $fig->delete_feature($user,$fid)

Invoking this routine deletes the feature designated by $fid.

=cut

sub delete_feature {
    my($self,$user,$fid) = @_;

    my $genome = &genome_of($fid);
    $self->log_update($user,$genome,$self->genus_species($genome),"Deleted Feature $fid");
    my $type   = &ftype($fid);
    my $dbh = $self->db_handle();
    my $file = $self->table_exists('deleted_fids') ? "$FIG_Config::organisms/$genome/Features/$type/deleted.features"
                                                   : "$FIG_Config::global/deleted.features";
    if (open(TMP,">>$file"))
    {
        flock(TMP,LOCK_EX) || confess "cannot lock deleted.features";
        print TMP "$fid\n";
        close(TMP);
        chmod 0777, $file;
    }
    if ($file eq "$FIG_Config::organisms/$genome/Features/$type/deleted.features")
    {
        $dbh->SQL("INSERT INTO deleted_fids (genome,fid) VALUES ('$genome','$fid')");
    }
    $self->{_deleted_fids}->{$fid} = 1;
}

sub undelete_feature {
    my($self,$user,$fid) = @_;

    my $genome = &genome_of($fid);
    $self->log_update($user,$genome,$self->genus_species($genome),"Undeleted Feature $fid");

    my $type   = &ftype($fid);
    my $dbh = $self->db_handle();
    &undelete_from_file($fid,"$FIG_Config::global/deleted.features");
    &undelete_from_file($fid,"$FIG_Config::organisms/$genome/Features/$type/deleted.features");

    if ($self->table_exists('deleted_fids'))
    {
        $dbh->SQL("DELETE FROM deleted_fids WHERE fid = '$fid'");
    }
    $self->{_deleted_fids}->{$fid} = 0;
}

# This is not done properly - the possibility of destructive overlap is obvious.  I doubt that
# it will be called 10 times during the lifetime of the SEED.  (RAO)

sub undelete_from_file {
    my($fid,$file) = @_;

    my $fidQ = quotemeta $fid;
    my @old = grep { $_ !~ /$fidQ$/ } `cat $file`;
    if (open(OLDDEL,">$file"))
    {
        foreach my $line (@old)
        {
            print OLDDEL $line;
        }
        close(OLDDEL);
    }
}


=head3 add_feature

    my $fid = $fig->add_feature($user,$genome,$type,$location,$aliases,$translation,$fid);

Invoking this routine adds the feature, returning a new (generated) $fid. It is
also possible to specify the feature ID, which is recommended if the feature is
to be permanent. (In order to do this the ID needs to be allocated from the
clearinghouse machine.) The translation is optional and only applies to PEGs.

=over 4

=item genome

ID of the genome to which the feature belongs.

=item type

Type of the feature (C<peg>, C<rna>, etc.)

=item location

Location of the feature, in the form of a comma-delimited list of location specifiers.
These are of the form I<contig>C<_>I<begin>C<_>I<end>, where I<contig> is the ID
of a contig, and I<begin> and I<end> are the starting and stopping offsets of the
location. These offsets are 1-based, and depending on the strand, the beginning
offset could be larger than the ending offset.

=item aliases

A comma-delimited list of alias names for the feature.

=item translation (optional)

The protein translation of the feature, if it is a peg.

=item fid (optional)

The ID to give to the new feature. If this parameter is omitted, an ID will be
generated automatically.

=item RETURN

Returns the new feature's ID if successful,or C<undef> if an error occurred.

=back

=cut

sub add_feature {
    my( $self, $user, $genome, $type, $location, $aliases, $sequence ) = @_;

    my $dbh = $self->db_handle();

    if ( $genome !~ /^\d+\.\d+$/ )
    {
        print STDERR "SEED error: add_feature failed due to bad genome id: $genome\n";
        return undef;
    }

    if ( $type !~ /^[0-9A-Za-z_]+$/ )
    {
        print STDERR "SEED error: add_feature failed due to bad type: $type\n";
        return undef;
    }

    if ( length ( $location ) > 5000 )
    {
        print STDERR "SEED error: add_feature failed because location is over 5000 char:\n";
        print STDERR "$location\n";
        return undef;
    }

    my @loc  = split( /,/, $location );
    my @loc2 = grep { defined($_->[0]) && $_->[1] && $_->[2] }
               map  { [ $_ =~ m/^(.+)_(\d+)_(\d+)$/ ] }
               @loc;

    if ( (@loc2 == 0) || ( @loc != @loc2 ) )
    {
        print STDERR "SEED error: add_feature failed because location is missing or malformed:\n";
        print STDERR "$location\n";
        return undef;
    }

    if ( my @bad_names = grep { length( $_->[0] ) > 96 } @loc2 )
    {
        print STDERR "SEED error: add_feature failed because location contains a contig name of over 96 char:\n";
        print STDERR join( ", ", @bad_names ) . "\n";
        return undef;
    }

    #  We should never recreate an existing feature:

    my ( $contig, $beg, $end );
    $contig = $loc2[0]->[0];
    $beg    = $loc2[0]->[1];
    my @same_contig = grep { $_->[0] eq $contig } @loc2;
    $end = $same_contig[-1]->[2];
    if ( $beg > $end )  { ( $beg, $end ) = ( $end, $beg ) }
    my ( $features, undef, undef ) = $self->genes_in_region( $genome, $contig, $beg, $end );

    my @same_loc = grep { scalar $self->feature_location( $_ ) eq $location }    # Same location
                   grep { /\.$type\.\d+$/ }                                      # Same type
                   @$features;                                                   # Near by features

    if ( @same_loc )
    {
        print STDERR "SEED Note: Attempt to recreate feature $same_loc[0]\n";
        return $same_loc[0];
    }

    my %seen = ();
    my @checksums = map { [ $_, $self->contig_md5sum( $genome, $_ ) ] }
                        grep { $_ && ( ! $seen{ $_ }++ ) }
                        map  { $_->[0] }
                        @loc2;
    my $fid = $self->fid_from_clearinghouse( $genome, $type, $location, \@checksums, $sequence );

    if ( ! $fid )
    {
        print STDERR "Failed to get a fid for $genome.$type at $location\n";
        return undef;
    }

    my ( $fidN ) = $fid =~ m/^fig\|\d+\.\d+\.[0-9A-Za-z_]+\.(\d+)$/;
    if ( ! $fidN || length( $fid ) > 32 )
    {
        print STDERR "SEED error: add_feature failed because the identifier is malformed or over 32 char: $fid\n";
        return undef;
    }

    $sequence ||= "";
    $aliases  ||= "";

    if ( 0 )   # GJO - Debug
    {
        print STDERR "SEED: Creating feature:\n"
                   . "   fid      = $fid\n"
                   . "   fidN     = $fidN\n"
                   . "   type     = $type\n"
                   . "   genome   = $genome\n"
                   . "   location = $location\n"
                   . "   contig   = $contig\n"
                   . "   minloc   = $beg\n"
                   . "   maxloc   = $end\n"
                   . "   aliases  = $aliases\n"
                   . "   sequence = $sequence\n";
    }

    if ($self->is_deleted_fid($fid))
    {
        $self->undelete_feature($user,$fid);
        $self->log_update($user,$genome,$self->genus_species($genome),"Undeleted Feature $fid");
	$self->broker_log('add_feature', {
	    user => $user,
	    genome => $genome,
	    fid => $fid,
	    type => 'undelete',
	    location => $location,
	    aliases => $aliases,
	    sequence => $sequence,
	});
        return $fid;
    }

    $self->broker_log('add_feature', {
	user => $user,
	genome => $genome,
	fid => $fid,
	type => 'new',
	location => $location,
	aliases => $aliases,
	sequence => $sequence,
    });

    $self->log_update($user,$genome,$self->genus_species($genome),"Added Feature $fid at $contig\_$beg\_$end");

    my $aliasesT = $aliases;
    $aliasesT =~ s/,\s*/\t/g;
    &add_tbl_entry( $fid, $location, $aliasesT );

    if ( $sequence )
    {
        $self->add_sequence( $fid, $sequence );
    }

    if ( ( $type eq "peg" ) and $sequence )
    {
        $self->enqueue_similarities([$fid]);

        my $md5 = Digest::MD5::md5_hex( uc $sequence );
        $dbh->SQL( "INSERT INTO protein_sequence_MD5 ( id, md5 ) VALUES ( '$fid', '$md5' )" );
    }

    my $rv = $dbh->SQL("INSERT INTO features (id,idN,type,genome,location,contig,minloc,maxloc,aliases)
                        VALUES ('$fid',$fidN,'$type','$genome','$location','$contig',$beg,$end,'$aliases')");

    my @aliases = split( /,\s*/, $aliases );
    my $alias;
    foreach $alias ( grep { /^(NP_|gi\||sp\||tr\||kegg\||uni\|)/ } @aliases )
    {
        $dbh->SQL("INSERT INTO ext_alias (id,alias,genome)
                   VALUES ('$fid','$alias','$genome')");
    }

    return $fid;
}


sub fid_from_clearinghouse
{
    my($self, $genome, $type, $location, $checksums, $translation) = @_;

    my $ch_url = "http://clearinghouse.theseed.org/Clearinghouse/clearinghouse_services.cgi";

    my $proxy = SOAP::Lite->uri("http://www.soaplite.com/Scripts")->proxy($ch_url);

    my $resp;
    eval {
        $resp = $proxy->add_feature($genome, $type, $location, $checksums, $translation);
    };
    if ($@)
    {
        warn "Error on proxy call: $@\n";
        return undef;
    }
    if ($resp->fault)
    {
        die "Failure on add_feature(genome=$genome,type=$type,location=$location): " .$resp->faultcode . ": " . $resp->faultstring . "\n";
    }

    return $resp->result;
}

=head3 clearinghouse_next_feature_id

    my $val = $fig->clearinghouse_next_feature_id($genome, $type)

Return the next feature ID that would be allocated by the clearinghouse for the given
genome and feature type.

=cut

sub clearinghouse_next_feature_id
{
    my($self, $genome, $type) = @_;

    my $ch_url = "http://clearinghouse.theseed.org/Clearinghouse/clearinghouse_services.cgi";
    my $proxy = SOAP::Lite->uri("http://www.soaplite.com/Scripts")->proxy($ch_url);

    my $resp;
    eval {
        $resp = $proxy->get_next_feature_id($genome, $type);
    };
    if ($@)
    {
        warn "Error on proxy call: $@\n";
        return undef;
    }
    if ($resp->fault)
    {
        warn "Failure on get_next_feature_id($genome, $type): " .$resp->faultcode . ": " . $resp->faultstring . "\n";
        return undef;
    }

    return $resp->result;
}

=head3 clearinghouse_register_metagenome_taxon_id

    my $tax = $fig->clearinghouse_register_metagenome_taxon_id($username, $genome_name)

Register a new taxon id for the MG-RAST metagenome server.

=cut

sub clearinghouse_register_metagenome_taxon_id
{
    my($self, $username, $genome_name) = @_;

    my $ch_url = "http://clearinghouse.theseed.org/Clearinghouse/clearinghouse_services.cgi";
    my $proxy = SOAP::Lite->uri("http://www.soaplite.com/Scripts")->proxy($ch_url);

    my $resp;
    eval {
        $resp = $proxy->register_metagenome_taxon_id($username, $genome_name);
    };
    if ($@)
    {
        warn "Error on proxy call: $@\n";
        return undef;
    }
    if ($resp->fault)
    {
        warn "Failure on register_metagenome_taxon_id($username, $genome_name): " .$resp->faultcode . ": " . $resp->faultstring . "\n";
        return undef;
    }

    return $resp->result;
}

=head3 clearinghouse_register_subsystem_id

    my $tax = $fig->clearinghouse_register_subsystem_id($ss_name);

Return a subsystem's short ID. Short IDs are maintained at a special
clearinghouse web site. If the subsystem does not yet have a short ID, a
new one will be assigned by the clearinghouse and returned.

=over 4

=item ss_name

Name of the subsystem whose ID is desired.

=item RETURN

ID of the desired subsystem.

=back

=cut

sub clearinghouse_register_subsystem_id {
    # Get the parameters.
    my ($self, $ss_name) = @_;
    # Return the result.
    return FIGRules::clearinghouse_register_subsystem_id($ss_name);
}


=head3 clearinghouse_lookup_subsystem_by_id

    my $tax = $fig->clearinghouse_lookup_subsystem_by_id($ss_name)

Register a subsystem id for the given subsystem name. Returns the existing id
if already present.

=cut

sub clearinghouse_lookup_subsystem_by_id
{
    my($self, $id) = @_;

    my $ch_url = "http://clearinghouse.theseed.org/Clearinghouse/clearinghouse_services.cgi";
    my $proxy = SOAP::Lite->uri("http://www.soaplite.com/Scripts")->proxy($ch_url);

    my $resp;
    eval {
        $resp = $proxy->lookup_subsystem_by_id($id);
    };
    if ($@)
    {
        warn "Error on proxy call: $@\n";
        return undef;
    }
    if ($resp->fault)
    {
        warn "Failure on lookup_subsystem_by_id($id): " .$resp->faultcode . ": " . $resp->faultstring . "\n";
        return undef;
    }

    return $resp->result;
}

=head3 clearinghouse_register_features

    my $val = $fig->clearinghouse_register_features($genome, $type, $num)

Register $num new features of type $type on genome $genome. Returns the starting index for the
new features.

=cut

sub clearinghouse_register_features
{
    my($self, $genome, $type, $num) = @_;

    my $ch_url = "http://clearinghouse.theseed.org/Clearinghouse/clearinghouse_services.cgi";
    my $proxy = SOAP::Lite->uri("http://www.soaplite.com/Scripts")->proxy($ch_url);

    my $resp;
    eval {
        $resp = $proxy->register_feature($genome, $type, $num);
    };
    if ($@)
    {
        warn "Error on proxy call: $@\n";
        return undef;
    }
    if ($resp->fault)
    {
        warn "Failure on register_feature($genome, $type, $num): " .$resp->faultcode . ": " . $resp->faultstring . "\n";
        return undef;
    }

    return $resp->result;
}


sub next_fid {
    my($self,$genome,$type) = @_;

    my $dbh = $self->db_handle();
    my $res = $dbh->SQL("select max(idN) from features where (genome = '$genome' and type = '$type')");
    return undef unless $res;

    my $fidN = $res->[0]->[0] + 1;
    while ($self->is_deleted_fid("fig\|$genome\.$type\.$fidN"))
    {
        $fidN++;
    }
    return "fig\|$genome\.$type\.$fidN";
}

sub replace_features_with {
    my($self,%args) = @_;

      my( $old_fids,  $user,   $genome,   $type,   $location,   $translation,   $function,  $fid) =
    @args{'old_fids', 'user',  'genome',  'type',  'location',  'translation',  'function', 'fid'};

    if ($old_fids =~ /^fig\|\d+\.\d+/) { $old_fids = [$old_fids] }
    if ((ref($old_fids) ne "ARRAY") || (@$old_fids < 1))  { return undef }

    if (! $user)     { return undef }
    if (! $genome)   { $genome = &FIG::genome_of($old_fids->[0]) }
    if (! $type  )   { $old_fids->[0] =~ /^fig\|\d+\.\d+\.([^\.]+)/; $type = $1 }
    if (! $location) { return undef }
    if (! $function) { $function = $self->function_of($old_fids->[0]) }

    my %aliases;
    foreach my $old_fid (@$old_fids)
    {
	if (($genome ne &genome_of($old_fid)) || ($type ne &ftype($old_fid))) { return undef }
	my @aliases = $self->feature_aliases($old_fid);
	foreach my $alias (@aliases)
	{
	    $aliases{$alias} = 1;
	}
    }
    my $all_aliases = join(",",sort keys(%aliases));

    my $new_fid = $self->add_feature($user,$genome,$type,$location,$all_aliases,$translation,$fid);

    if ($new_fid)
    {
	if ($function)
	{
	    $self->assign_function($new_fid,$user,$function);
	}
	foreach my $old_fid (@$old_fids)
	{
	    next if ($new_fid eq $old_fid);
	    $self->inherit_annnotations($old_fid,$new_fid);
	    $self->add_annotation($old_fid,$user,"Replaced by $new_fid");
	    $self->delete_feature($user,$old_fid);
	}
	my $all_old = join(",",@$old_fids);
	$self->log_update($user,$genome,$self->genus_species($genome),"Replaced Features $all_old with $new_fid");
    }
    return $new_fid;
}

sub inherit_annnotations {
    my($self,$old_fid,$new_fid) = @_;

    my @annotations = $self->feature_annotations($old_fid,"rawtime");
    foreach my $ann (@annotations)
    {
	my(undef, $timeStamp, $user, $annotation) = @$ann;
	$self->add_annotation($new_fid,$user,"Inherited from $old_fid\n\n$annotation",$timeStamp);
    }
}

=head3 is_deleted_fid_bulk

    my $hash = $fig->is_deleted_fid_bulk(@fids)

Returns a hash { $fid => 0|1  } based on the deleted status in the database.

=cut

sub is_deleted_fid_bulk
{
    my($self, @fids) = @_;

    my $del = {};
    $del->{$_} = 0 foreach @fids;
    my @fids2 = map { $_->[0] } grep { !(defined($_->[1] && $self->is_genome($_->[1]))) }
    		      map { /^fig\|(\d+\.\d+)\./ ? [$_, $1] : [$_] } @fids;

    $del->{$_} = 1 foreach @fids2;
    @fids = grep { !$del->{$_} } @fids;

    my @ok;

    my $dbh = $self->db_handle();

    while (@fids)
    {
	my @list = splice(@fids, 0, 1000);
	my $q = "?" . (",?" x (@list - 1));
	my $res = $dbh->{_dbh}->selectall_arrayref(qq(SELECT fid
						      FROM deleted_fids
						      WHERE fid IN ($q)),
						   undef, @list);
	$del->{$_->[0]} = 1 foreach @$res;
    }
    return $del;
}

sub is_deleted_fid {
    my($self,$fid) = @_;
    my($x,$y);

    #
    # If we have the FIG_Config deleted_fids_database_only flag turned on,
    # just use the database via valid_fids to handle this.
    #
    if ($FIG_Config::deleted_fids_database_only)
    {
	my $del = $self->is_deleted_fid_bulk($fid);
	return $del->{$fid};
    }

    if ($fid !~ /^fig\|\d+\.\d+\./) { return 0 }

    if (! ($x = $self->{_deleted_fids}))
    {
        $x = $self->{_deleted_fids} = {};
    }

    if (defined($y = $x->{$fid}))
    {
        return $y;
    }
    if (! $self->is_genome(&genome_of($fid)))
    {
        $x->{$fid} = 1;
        return 1;
    }

    #
    # If we've loaded the table, and it's not there, it's not deleted.
    #
    if ($self->{_deleted_fids_loaded})
    {
        return 0;
    }

    my $dbh = $self->db_handle();
    if (! $self->table_exists('deleted_fids'))
    {
        $dbh->create_table(tbl => 'deleted_fids',flds => 'genome varchar(16), fid varchar(32)');
        my $tmpfile = "$FIG_Config::temp/delfids$$";
        if ((-s "$FIG_Config::global/deleted.features") && open(TMP,">$tmpfile"))
        {
            open(GLOBDEL,"<$FIG_Config::global/deleted.features") || die "I could not open $FIG_Config::global/deleted.features";
            while (defined($y = <GLOBDEL>))
            {
                if ($y =~ /^fig\|(\d+\.\d+)/)
                {
                    print TMP "$1\t$y";
                }
            }
            close(GLOBDEL);
            close(TMP);
            $dbh->load_table(tbl => 'deleted_fids', file => $tmpfile, delim => "\t" );
            $dbh->create_index( tbl => 'deleted_fids', idx => 'deleted_fids_fid_idx', flds => 'fid');
            $dbh->create_index( tbl => 'deleted_fids', idx => 'deleted_fids_genome_idx', flds => 'genome');
            unlink($tmpfile);
        }
    }

    #
    # Cache the whole darn deleted table.
    #

    $self->{_deleted_fids_loaded} = 1;
    my $res = $dbh->SQL("SELECT fid FROM deleted_fids");
    map { $x->{$_->[0] } = 1 } @$res;

    return $x->{$fid};

    $res = $dbh->SQL("SELECT fid FROM deleted_fids WHERE fid = '$fid'");
    my $deleted = (@$res > 0);
    $x->{$fid} = $deleted;
    return $deleted;
}

sub fid_with_changed_location {
    my($self,$fid) = @_;
    my($x);

    if (! ($x = $self->{_changed_location_fids}))
    {
        $self->{_changed_location_fids} = {};
        if (open(TMP,"<$FIG_Config::global/changed.location.features"))
        {
            while ($_ = <TMP>)
            {
                if ($_ =~ /^(fig\|\d+\.\d+\.[a-zA-Z]+\.\d+)/)
                {
                    $self->{_changed_location_fids}->{$1}++;
                }
            }
            close(TMP);
        }
        $x = $self->{_changed_location_fids};
    }
    return $x->{$fid};
}


=head3 call_start

usage: $fig->call_start($genome,$loc,$translation,$against)

This routine can be invoked to produce an estimate of the correct start, given
a location in a genome believed to be a protein-encoding gene, along with a set
of PEGs that are believed to be orthologs.  If called in a list context,
it returns a list containing

    a string representing the estimated start location
    a confidence measure (better than 0.2 seems to be pretty solid)
    a new translation

If called in a scalar context, it returns its best prediction of the start.

=cut

sub call_start {
    my($self,$genome,$loc,$tran,$against) = @_;
    my($peg);

    my $orgdir;
    my $orgdir_flag = "";
    if (defined($orgdir = $self->{_orgdir})) {
	$orgdir_flag = "-orgdir=$orgdir";
    }

    open(TMP,"| start_data_for_set_of_pegs $orgdir_flag use-close > $FIG_Config::temp/tmp.objects$$")
        || die "could not set up pipe to start_data_for_set_of_pegs";
    print TMP "new|$genome\.peg\.1\t$loc\t$tran\n";
    foreach $peg (@$against)
    {
        print TMP "$peg\tno_recall\n";
    }
    close(TMP);

    &FIG::run("predict_starts $orgdir_flag $FIG_Config::temp/tmp.objects$$ $FIG_Config::temp/clear$$ $FIG_Config::temp/proposed$$ > /dev/null");

    if (-s "$FIG_Config::temp/proposed$$")
    {
        my @changes = `changed_starts $FIG_Config::temp/proposed$$ /dev/null`;
        if ((@changes == 1) && ($changes[0] =~ /^\S+\t\S+\t(\S+)\t(\S+)/))
        {
            my($new_loc,$conf) = ($1,$2);
            if (($ENV{FIG_VERBOSE}) && open(TMP,"<$FIG_Config::temp/proposed$$"))
            {
                while (defined($_ = <TMP>)) { print STDERR $_ }
                close(TMP);
            }
            my $proposed = wantarray ? join("",`cat $FIG_Config::temp/proposed$$`) : "";
            $proposed =~ s/^ID=[^\n]+\n//s;
            unlink("$FIG_Config::temp/tmp.objects$$","$FIG_Config::temp/clear$$","$FIG_Config::temp/proposed$$");
            return wantarray ? ($new_loc,$conf,$self->fixed_translation($tran,$genome,$loc,$new_loc),$proposed) : $new_loc;
        }
    }
    unlink("$FIG_Config::temp/tmp.objects$$","$FIG_Config::temp/clear$$","$FIG_Config::temp/proposed$$");
    return wantarray ? ($loc,0,$tran,"") : $loc;
}

sub fixed_translation {
    my($self,$old_tran,$genome,$old_loc,$new_loc) = @_;
    my($extra,$trimmed,$new_tran);

    if ($old_loc =~ /^(\S+)_(\d+)_(\d+)$/)
    {
        my($contigO,$begO,$endO) = ($1,$2,$3);

        if ($new_loc =~ /^(\S+)_(\d+)_(\d+)$/)
        {
            my($contigN,$begN,$endN) = ($1,$2,$3);
            if ($begO < $endO)
            {
                if ($begO < $begN)
                {
                    $trimmed = ($begN - $begO) / 3;
                    $new_tran = &translate($self->dna_seq($genome,join("_",($contigO,$begN,$begN+2))),undef,"start") .
                                substr($old_tran,$trimmed+1);
                }
                else
                {
                    $extra = ($begO - $begN) / 3;
                    $new_tran = &translate($self->dna_seq($genome,join("_",($contigO,$begN,$begO+2))),undef,"start") .
                                substr($old_tran,1);
                }
            }
            else
            {
                if ($begO > $begN)
                {
                    $trimmed = ($begO - $begN) / 3;
                    $new_tran = &translate($self->dna_seq($genome,join("_",($contigO,$begN,$begN-2))),undef,"start") .
                                substr($old_tran,$trimmed+1);
                }
                else
                {
                    $extra = ($begN - $begO) / 3;
                    $new_tran = &translate($self->dna_seq($genome,join("_",($contigO,$begN,$begO-2))),undef,"start") .
                                substr($old_tran,1);
                }
            }
            return $new_tran;
        }
    }
    return $old_tran;
}


=head3 pick_gene_boundaries

usage: $fig->pick_gene_boundaries($genome,$loc,$translation)

This routine can be invoked to expand a region of similarity to potential
gene boundaries.  It does not try to find the best start, but only the one that
is first after the beginning of the ORF.  It returns a list containing
the predicted location and the expanded translation.  Thus, you might use

($new_loc,$new_tran) = $fig->pick_gene_boundaries($genome,$loc,$tran);
$recalled            = $fig->call_start($genome,$new_loc,$new_tran,\@others);

to get the location of a recalled gene (in, for example, the process of correcting
a frameshift).

=cut

sub pick_gene_boundaries {
    my($self, $genome, $loc, $tran, $genetic_code, $search_region) = @_;

    if ($ENV{FIG_VERBOSE}) {
	print STDERR "Picking gene boundaries for org=$genome, loc=$loc";
	print STDERR ", tran=$tran" if $tran;
	print STDERR ", genetic_code=$genetic_code"   if $genetic_code;
	print STDERR ", search_region=$search_region" if $search_region;
	print STDERR "\n";
    }

    my($leftStop, $firstStart, $start, $end, $rightStop);

    my $full_loc = new FullLocation($self, $genome, $loc, $tran, $genetic_code);
    if ($full_loc->PickGeneBoundaries(-limit => $search_region)) {
	return ($full_loc->SeedString(), $full_loc->Translation());
    }

    return (undef, undef);
}



=head3 change_location_of_feature

usage: $fig->change_location_of_feature($fid,$location,$translation)

Invoking this routine changes the location of the feature.  The $translation argument
is optional (and applies only to PEGs).

The routine returns 1 on success and 0 on failure.

=cut

sub change_location_of_feature {
    my($self,$fid,$location,$translation) = @_;
    my($x);
    if ($self->is_deleted_fid($fid)) { return 0 }

    my $dbh = $self->db_handle();

    my $genome = &genome_of($fid);
    my $type   = &ftype($fid);

    my($got) = 0;
    my @loc = split(/,/,$location);
    my($contig,$beg,$end);
    if (($loc[0] =~ /^(\S+)_(\d+)_\d+$/) && (($contig,$beg) = ($1,$2)) && ($location =~ /(\d+)$/))
    {
        $end = $1;
        if ($beg > $end)  { ($beg,$end) = ($end,$beg) }
    }
    else
    {
        return 0;
    }

    my @tmp = grep { ($_ =~ /^(\S+)/) && ($1 eq $fid) } `grep '$fid' $FIG_Config::organisms/$genome/Features/$type/tbl`;
    if (@tmp > 0)
    {
        $x = $tmp[$#tmp];
        chop $x;
        my @flds = split(/\t/,$x);
        shift @flds;
        shift @flds;
        my $aliasesT =  (@flds > 0) ? join("\t",@flds) : "";
        &add_tbl_entry($fid,$location,$aliasesT);

        $dbh->SQL("UPDATE features SET location = '$location',
                                       contig = '$contig',
                                       minloc = $beg,
                                       maxloc = $end
                                   WHERE id = '$fid'");

        if (my $locations = $self->cached('_location'))
        {
            $locations->{$fid} = $location;
        }

        open(TMP,">>$FIG_Config::global/changed.location.features")
            || die "could not open $FIG_Config::global/changed.location.features";

        flock(TMP,LOCK_EX) || confess "cannot lock changed.location.features";
        print TMP "$fid\n";
        close(TMP);
        chmod 0777, "$FIG_Config::global/changed.location.features";
        $self->{_changed_location_fids} = undef;

        if (($type eq "peg") && defined($translation))
        {
            $self->add_sequence($fid,$translation);
        }
        $got = 1
    }
    return $got;
}

sub add_tbl_entry {
    my($fid,$location,$aliasesT) = @_;

    my $type;
    if ($fid =~ /^fig\|\d+\.\d+\.([a-zA-Z0-9_-]+)/)
    {
	$type = $1;
    }
    else
    {
	return "";
    }
    my $genome = &genome_of($fid);
    &verify_dir("$FIG_Config::organisms/$genome/Features/$type");
    my $file   = "$FIG_Config::organisms/$genome/Features/$type/tbl";
    open(TMP,">>$file")
        || die "could not open $file";
    flock(TMP,LOCK_EX) || confess "cannot lock $file";
    print TMP "$fid\t$location\t$aliasesT\n";
    close(TMP);
    chmod 0777, "$file";
}


sub add_sequence {
    my($self,$fid,$seq) = @_;

    my $type;
    if ($fid =~ /^fig\|\d+\.\d+\.([a-zA-Z0-9_-]+)/)
    {
	$type = $1;
    }
    else
    {
	return "";
    }
    my $genome = &genome_of($fid);
    &verify_dir("$FIG_Config::organisms/$genome/Features/$type");
    my $file   = "$FIG_Config::organisms/$genome/Features/$type/fasta";
    if (open(TMP,">>$file"))
    {
        flock(TMP,LOCK_EX) || confess "cannot lock $file";
        print TMP ">$fid\n";
        my $seek = tell TMP;
        my $ln   = length($seq);
        print TMP "$seq\n";
        close(TMP);
        chmod 0777, $file;
        my $fileno = $self->file2N($file);

	if ($type eq "peg")
	{
	    my $dbh = $self->db_handle();
	    $dbh->SQL("DELETE FROM protein_sequence_seeks where id = '$fid'");
	    $dbh->SQL("INSERT INTO protein_sequence_seeks (id,fileno,seek,len,slen)
                              VALUES ('$fid',$fileno,$seek,$ln+1,$ln)");
	}
    }
}

sub peg_in_gendb
{
    my($self, $peg) = @_;
    my $genome = &genome_of($peg);

    return $self->genome_in_gendb($genome);
}

sub genome_in_gendb
{
    my($self, $genome) = @_;
    return 0;
}

### Some rendering stuff
#

=head2 genome_to_gg

Render a genome's contig as GenoGraphics objects.

=cut

sub genome_to_gg
{
    my($self, $genome, $contig, $width) = @_;

    my $gg = [];

    my $len = $self->contig_ln($genome, $contig);

    my $next_color = 0;
    my %sub_color;

    for (my $start = 0; $start + $width < $len; $start += $width)
    {
        my $label = $start;
        my $end = $start + $width;

        my($genes, $g_beg, $g_end) = $self->genes_in_region($genome, $contig, $start, $end);

        my $map = [];

        for my $gene (@$genes)
        {
            my $loc = $self->feature_location($gene);
            my($c, $b, $e) = $self->boundaries_of($loc);

            my $shape;

            if ($b < $e)
            {
                $shape = "rightArrow";
            }
            else
            {
                $shape = "leftArrow";
                ($b, $e) = ($e, $b);
            }

            my($type, $peg_n) = ($gene =~ /fig\|\d+\.\d+\.(\w+)\.(\d+)$/);

            my $color = "red";
            if ($type eq 'rna')
            {
                $color= 'black';
            }

            my @a = $self->feature_aliases($gene);
            my @gene_names = grep { /^[a-zA-Z]{4}$/ } @a;
            if (@gene_names)
            {
                $peg_n = $gene_names[0];
            }

            my @subs = $self->peg_to_subsystems($gene);
            if (@subs)
            {
                my $sub = $subs[0];
                if (not exists $sub_color{$sub})
                {
                    my $c = $next_color + 1;
                    $next_color = ($next_color + 1) % 20;
                    $sub_color{$sub} = "color$c";
                }
                $color = $sub_color{$sub};
            }

            $b = $start if $b < $start;
            $e = $end if $e > $end;

            push(@$map, [$b - $start, $e - $start, $shape, $color, $peg_n, '', '']);
        }

        push(@$gg, [$label, 0, $width, $map]);
    }

    for my $sub (sort keys %sub_color)
    {
        my $map = [3000,  $width - 10,  'rect', $sub_color{$sub}, $sub, '', ''];
        push(@$gg, ['', 0, $width, $map]);
    }
    return $gg;
}

=head2 Markup Helper Methods

This section contains the methods used to read and write Markup data. The
markup data associates labels with sections of a feature's translation.

In the SEED, Markup data is stored in a separate file for each marked feature
in the the feature type subdirectory for an organism. So, for example, the
PEG markups for C<fig|83333.1.peg.4> would be in the file

    FIG/Data/Organisms/83333.1/peg/markup4.tbl

The file is stored in tab-separated form. Each line contains the following
fields

=over 4

=item start

1-based offset into the translation of the first amino acid to mark

=item len

number of amino acids to mark

=item label

label identifying the type of markup

=back

Reading and writing these tiny files is extremely fast, but they do have more
overhead than would be expected if the data were stored in a single flat file
managed by pointers from the FIG database. If that apprach becomes desirable,
then only this section of FIG.pm needs to be changed.

=cut

#

=head3 ReadMarkups

    my $marks = $fig->ReadMarkups($fid);

Read the markup data for the specified feature. The markings are returned as a
list of triples. Each triple contains the start location of a markup, the
length of the markup, and the label.

=over 4

=item fid

ID of the feature whose markups are to be read.

=item RETURN

Returns a reference to list of 3-tuples. Each list element will consist of the
starting offset of the markup (1-based), the length of the markup, and the label.
All values are expressed in terms of distance into the protein translation of the
feature.

=back

=cut
#: Return Type $@@;
sub ReadMarkups {
    # Get the parameters.
    my ($self, $fid) = @_;
    # Declare the return variable.
    my $retVal = [];
    # Get the name of the markup file.
    my $fname = _MarkupFileName($fid);
    # Get the contents of the file and parse it.
    if (-e $fname) {
        push @{$retVal}, map { [ Tracer::ParseRecord($_) ] } Tracer::GetFile($fname);
    }
    # Return the result.
    return $retVal;
}

=head3 WriteMarkups

    $fig->WriteMarkups($fid, \@marks);

Write out the markups for the specified feature. If the markup file for the
specified feature does not exist, it will be created. If it does exist, it
will be completely overwritten.

=over 4

=item fid

ID of the feature whose markups are to be written

=item marks

Reference to a list of markups. Each markup is in the form of a 3-tuple consisting
of the 1-based offset to the start of the markup, the length of the markup, and
the markup label. The offset and length are specified in terms of the protein
translation string.

=back

=cut
#: Return Type ;
sub WriteMarkups {
    # Get the parameters.
    my ($self, $fid, $marks) = @_;
    # Locate the output file.
    my $fname = _MarkupFileName($fid);
    # Open it for output.
    Open(\*OUTMARKS, ">$fname");
    # Write out the mark data.
    for my $mark (@{$marks}) {
        print OUTMARKS join("\t", @{$mark}) . "\n";
    }
    # Close the output file.
    close OUTMARKS;
}

=head3 _MarkupFileName

    my $name = FIG::_MarkupFileName($fid);

Return the name of the file containing the markup data for the specified feature.

=over 4

=item fid

ID of the feature whose markup file is desired.

=item RETURN

Returns the full path of the file containing the feature markups for the feature desired.

=back

=cut
#: Return Type $;
sub _MarkupFileName {
    # Get the parameters.
    my ($fid) = @_;
    # Declare the return variable. We prime it with the organism directory.
    my $retVal = $FIG_Config::organisms;
    # Parse the feature ID.
    my ($genome, $type, $idx);
    if ($fid =~ /fig\|(\d+\.\d+)\.([a-z]+)\.(\d+)/) {
        ($genome, $type, $idx) = ($1,$2,$3);
    } else {
        Confess("Invalid feature ID $fid specified in Markup call.");
    }
    # Compute the file name from the pieces of the feature ID.
    $retVal .= "/$genome/Features/$type/markup$idx.tbl";
    # Return the result.
    return $retVal;
}

=head2 UserData Helper Methods

This section contains the methods used to implement UserData access. User data is
stored in a subdirectory given by the user's name under the C<Users> directory
in the Global directory tree. In other words, the data for the default user
C<basic> would be at C<$FIG_Config::global/Users/basic>.

In each directory, the C<capabilities.tbl> file contains the capability data and
the C<preferences.tbl> file contains the preferences. Currently, preferences are
stored in a single file, but if performance becomes a problem we may split them
by category.

Each of these files has two columns of data-- a key and a value. In the preferences
file the key is a hierarchical construct with the pieces separated by colons, and
the value is essentially a free-format string understood only by the application. In
the capabilities file the key is a group name, and the value is an access level--
C<RW> (full access), C<RO> (read-only access), or C<NO> (no access).

Group names and key names are not allowed to contain white space. Tabs are used to
separate them from the value strings or access levels. The value strings for
preferences cannot contain tabs or new-lines. A backslash escape mechanism
will be used to allow tabs and new-lines to be specified in the preference values.

The files are sorted by key, to make updates easier.

The special C<Security_Default> subdirectory is used to track the default security
options for each secure object. The object's security group and default level
are specified in a file whose name is formed by appending the object ID to the
object type with an extension of "tbl". So, for example, the file containing the
security default information for Genome 83333.1 would be

    $FIG_Config::global/Users/Security_Default/Genome_83333.1.tbl

Each of these is a tiny file with the group name and default access level for that
organism or subsystem. The two fields of the file are tab-separated, and any new-line
character at the end is ignored.

=head3 GetDefault

    my ($group, $level) = $fig->GetDefault($objectID, $objectType);

Return the group name and default access level for the specified object.

=over 4

=item objectID

ID of the object whose capabilities data is desired.

=item objectType

Type of the object whose capabilities data is desired. This should be expressed
as a Sprout entity name. Currently, the only types supported are C<Genome>
and C<Subsystem>.

=item RETURN

Returns a two-element list. The first element is the name of the group
to witch the object belongs; the second is the default access level
(C<RW>, C<RO>, or C<NO>). If the object is not found, an empty list
should be returned.

=back

=cut

sub GetDefault {
    # Get the parameters.
    my ($self, $objectID, $objectType) = @_;
    # Declare the return variable.
    my @retVal = ();
    # Compute the file name for this object.
    my $fileName = _GetObjectCapabilityFile($objectType, $objectID);
    # Only proceed if the file exists and has data.
    if ($fileName && -e $fileName) {
        # Open the file and read the first line.
        Open(\*DEFAULTIN, "<$fileName");
        # Read the first (and only) line of the file.
        @retVal = _GetInputKVRecord(\*DEFAULTIN);
        # Close the file.
        close DEFAULTIN;
    }
    # Return the result.
    return @retVal;
}

=head3 GetPreferences

    my $preferences = $fig->GetPreferences($userID, $category);

Return a map of preference keys to values for the specified user in the
specified category.

=over 4

=item userID

ID of the user whose preferences are desired.

=item category (optional)

Name of the category whose preferences are desired. If omitted, all
preferences should be returned.

=item RETURN

Returns a reference to a hash mapping each preference key to a value. The
keys are fully-qualified; in other words, the category name is included.
It is acceptable for the hash to contain key-value pairs outside the
category. In other words, if it's easier for you to read the entire
preference set into memory, you can return that one set every time
this method is called without worrying about the extra keys.

=back

=cut

sub GetPreferences {
    # Get the parameters.
    my ($self, $userID, $category) = @_;
    # Get the preferences. Note we use the category name followed by a colon
    # (the official separator character) to restrict the preferences to the
    # ones we want.
    my %retVal = _GetUserDataFile($userID, 'preferences', "$category:");
    # Return the data.
    return \%retVal;
}

=head3 GetCapabilities

    my $level = $fig->GetCapabilities($userID);

Return a map of group names to access levels (C<RW>, C<RO>, or C<NO>) for the
specified user.

=over 4

=item userID

ID of the user whose access level is desired.

=item RETURN

Returns a reference to a hash mapping group names to the user's access level
for that group.

=back

=cut

sub GetCapabilities {
    # Get the parameters.
    my ($self, $userID, $category) = @_;
    # Get the complete list of capabilities.
    my %retVal = _GetUserDataFile($userID, 'capabilities');
    # Return the data.
    return \%retVal;
}

=head3 AllowsUpdates

    my $flag = $fig->AllowsUpdates();

Return TRUE if this access object supports updates, else FALSE. If the access object
does not support updates, none of the B<SetXXXX> methods will be called.

=cut

sub AllowsUpdates {
    return 1;
}

=head3 SetDefault

    $fig->SetDefault($objectID, $objectType, $group, $level);

Set the group and default access level for the specified object.

=over 4

=item objectID

ID of the object whose access level and group are to be set.

=item objectType

Type of the relevant object. This should be expressed as a Sprout entity name.
Currently, only C<Genome> and C<Subsystem> are supported.

=item group

Name of the group to which the object will belong. A user's access level for
this group will override the default access level.

=item level

Default access level. This is the access level used for user's who do not have
an explicit capability specified for the object's group.

=back

=cut

sub SetDefault {
    # Get the parameters.
    my ($self, $objectID, $objectType, $group, $level) = @_;
    # Find the target file.
    my $fileName = _GetObjectCapabilityFile($objectType, $objectID);
    if (! $fileName) {
        Confess("Invalid object $objectType ($objectID) specified in SetDefault.");
    } else {
        # Write out the new default data.
        Open(\*DEFAULTOUT, ">$fileName");
        _PutOutputKVRecord(\*DEFAULTOUT, $group, $level);
        close DEFAULTOUT;
    }
}

=head3 SetCapabilities

    $fig->SetCapabilities($userID, \%groupLevelMap);

Set the access levels by the specified user for the specified groups.

=over 4

=item userID

ID of the user whose capabilities are to be updated.

=item groupLevelMap

Reference to a hash that maps group names to access levels. The legal
access levels are C<RW> (read-write), C<RO> (read-only), and C<NO> (no
access). An undefined value for the access level indicates the default
level should be used for that group. The map will not replace all of
the user's capability date; instead, it overrides existing data, with
the undefined values indicating the specified group should be deleted
from the list.

=back

=cut

sub SetCapabilities {
    # Get the parameters.
    my ($self, $userID, $groupLevelMap) = @_;
    # Get the relevant file name.
    my $fileName = _GetUserDataDirectory($userID);
    # Insure this used is real.
    if (! $fileName) {
        Confess("Invalid user $userID specified when updating capabilities.");
    } else {
        # Process the updates.
        _ProcessUpdates("$fileName/capabilities.tbl", $groupLevelMap);
    }
}

=head3 SetPreferences

    $fig->SetPreferences($userID, \%preferenceMap);

Set the preferences for the specified user.

=over 4

=item userID

ID of the user whose preferences are to be udpated.

=item preferenceMap

Reference to a hash that maps each preference key to its value. The
keys should be fully-qualified (that is, they should include the
category name). A preference key mapped to an undefined value will
use the default preference value for that key. The map will not
replace all of the user's preference data; instead, it overrides
existing data, with the undefined values indicating the specified
preference should be deleted from the list.

=back

=cut

sub SetPreferences {
    # Get the parameters.
    my ($self, $userID, $preferencesMap) = @_;
    # Get the relevant file name.
    my $fileName = _GetUserDataDirectory($userID);
    # Insure this user is real.
    if (! $fileName) {
        Confess("Invalid user $userID specified when updating capabilities.");
    } else {
        # Process the updates.
        _ProcessUpdates("$fileName/preferences.tbl", $preferencesMap);
    }
}

=head3 CleanupUserData

    $fig->CleanupUserData();

Release any data being held in memory for use by the UserData object.

=cut

sub CleanupUserData {
    # There is no data to clean up.
}

=head2 UserData Utilities

=head3 GetObjectCapabilityFile

    my $fileName = FIG::_GetObjectCapabilityFile($objectType, $objectID);

This is an internal method that computes the name of the file containing the
default group and access data for a specified object. It returns the file
name.

=cut

sub _GetObjectCapabilityFile {
    # Get the parameters.
    my ($objectType, $objectID) = @_;
    # Clean name to insure it's valid.
    my $cleanObject = $objectID;
    $cleanObject =~ tr/: /__/;
    # Insure the security default directory exists.
    my $directory = "$FIG_Config::global/Users/Security_Default";
    Insure($directory);
    # Form the file name.
    my $retVal = "$directory/${objectType}_$cleanObject.tbl";
    # Return the result.
    return $retVal;
}

=head3 to_structured_english

    my ($ev_code_list, $subsys_list, $english_string) = $fig->to_structured_english($fig, $peg, $escape_flag);

Create a structured English description of the evidence codes for a PEG,
in either HTML or text format. In addition to the structured text, we
also return the subsystems and evidence codes for the PEG in list form.

=over 4

=item peg

ID of the protein or feature whose evidence is desired.

=item escape_flag

TRUE if the output text should be HTML, else FALSE

=item RETURN

Returns a three-element list. The first element is a reference to a list of evidence codes,
the second is a list of the subsystem containing the peg, and the third is the readable
text description of the evidence.

=back

=cut

sub to_structured_english {
    my ($self, $peg, $escape_flag, %options) = @_;
    return FIGRules::to_structured_english($self, $peg, $escape_flag, %options);
}

=head3 GetUserDataDirectory

    my $directoryName = FIG::_GetUserDataDirectory($userName);

Return the name of the directory containing the user's preference and capability
data. If the user does not have a directory, return C<undef>.

=over 4

=item userName

Name of the user whose directory is desired.

=item RETURN

Returns the name of the user's preference/capability directory. If the user does
not exist, the directory will be created automatically. If this policy is changed,
return C<undef> to indicate an invalid user name.

=back

=cut

sub _GetUserDataDirectory {
    # Get the parameters.
    my ($userName) = @_;
    # Compute the directory name.
    my $retVal = "$FIG_Config::global/Users/$userName";
    # Insure it exists.
    Insure($retVal);
    return($retVal);
}

=head3 GetUserDataFile

    my %userData = FIG::_GetUserDataFile($userID, $type, $prefix);

Create a hash from the user data file of the specified type. The user data file
contains two tab-delimited fields. The first field will be read in as the key
of the hash and the second as the data value. The file must be sorted, and
only records beginning with the character string in I<$prefix> will be put
in the hash.

=over 4

=item userID

Name of the user whose preference or capability data is desired.

=item type

Type of file desired: C<preferences> or C<capabilities>.

=item RETURN

Returns a hash containing all the key/value pairs in the user file of the
specified type. If the file is not found, will return an empty hash.

=back

=cut

sub _GetUserDataFile {
    # Get the parameters.
    my ($userID, $type, $prefix) = @_;
    # Declare the return value.
    my %retVal = ();
    # Try to find the user's directory.
    my $directory = _GetUserDataDirectory($userID);
    # Only proceed if it exists.
    if ($directory) {
        # Create the input file name.
        my $fileName = "$directory/$type.tbl";
        # If the file exists, we open it.
        if (-e $fileName) {
            Open(\*USERDATA, "<$fileName");
            # Use a null string for an undefined prefix, then compute the
            # minimum and maximum permissible key values. The EOF trick
            # works because keys should not contain non-ASCII characters.
            my $minKey = (defined $prefix ? $prefix : "");
            my $maxKey = $minKey . Tracer::EOF;
            # Read until we're done.
            my $done = 0;
            while (! $done) {
                # Get the next record.
                my ($key, $value) = _GetInputKVRecord(\*USERDATA);
                # Process according to the nature of the data on the line.
                if (! defined $key || $key ge $maxKey) {
                    # Here we're done. We've either hit end-of-file or
                    # the current line's key is too big.
                    $done = 1;
                } elsif ($key ge $minKey) {
                    # Here we want to keep the line.
                    $retVal{$key} = $value;
                }
            }
            # Close the file.
            close USERDATA;
        }
    }
    # Return the hash.
    return %retVal;
}

=head3 ProcessUpdates

    FIG::_ProcessUpdates($fileName, \%map);

Apply the specified updates to a key-value file. The records in the key-value file must
be sorted. If a key in the map matches a key in the file, the file's key value is
replaced. If a key in the map is not found in the file, it is added. If a key in the
map is found in the file and it has an undefined value in the map, then the key
is deleted.

=over 4

=item fileName

Name of the file to be updated.

=item map

Reference to a hash mapping keys to values. The keys may not contain any whitespace.
The value will be escaped before it is written.

=back

=cut

sub _ProcessUpdates {
    # Get the parameters.
    my ($fileName, $map) = @_;
    # Create a temporary file for the update.
    my $tmpFileName = "$fileName$$.tmp";
    # Get the map keys in lexical order.
    my @keys = sort keys %{$map};
    # Push on the EOF constant.
    push @keys, Tracer::EOF;
    # These variable will contain the key and value fields of the current
    # record of the input file.
    my ($lineKey, $lineValue) = (Tracer::EOF, undef);
    # If the input file does not exist, we pretend it's empty. Otherwise,
    # we read the first line.
    if (-e $fileName) {
        Open(\*USERDATAIN, "<$fileName");
        ($lineKey, $lineValue) = _GetInputKVRecord(\*USERDATAIN);
    }
    # Finally, we open the temp file for output.
    Open(\*USERDATAOUT, ">$tmpFileName");
    # Get the first key.
    my $key = shift @keys;
    # Loop until we reach the end of both lists.
    while ($key lt Tracer::EOF || $lineKey lt Tracer::EOF) {
        # Compare the keys to determine what to do next.
        if ($lineKey lt $key) {
            # Here we must read the next record. First we have to write
            # the previous one. Note that if $lineValue is undefined,
            # the record is discarded automatically.
            _PutOutputKVRecord(\*USERDATAOUT, $lineKey, $lineValue);
            ($lineKey, $lineValue) = _GetInputKVRecord(\*USERDATAIN);
        } elsif ($lineKey eq $key) {
            # Here we have a match. We select the new key's value as the
            # value of the line key and let the loop spin. When the key
            # is written to the output file, the new value will be used.
            # if the new value is undefined, the record is thrown away,
            # which is exactly what we want.
            $lineValue = $map->{$key};
            $key = shift @keys;
        } else {
            # Here the key in the map is new, so we write it to the
            # output file and get the next key.
            _PutOutputKVRecord(\*USERDATAOUT, $key, $map->{$key});
            $key = shift @keys;
        }
    }
    # Close the files.
    close USERDATAOUT;
    close USERDATAIN;
    # Replace the old file with the temporary. We delete the old file first so
    # that a rename is used for the move, which is safer.
    unlink $fileName;
    move($tmpFileName, $fileName);
}

=head3 GetInputKVRecord

    my ($key, $value) = FIG::_GetInputKVRecord($handle);

Read a key/value pair from the specified input file. If we are at end-of-file
the key returned will be the C<Tracer::EOF> constant. The key and value are
separated by a tab. The value will be unescaped if it exists.

=over 4

=item handle

Open handle for the input file.

=item RETURN

Returns a two-element list. The first element will be the first field of the
input record; the second element will be the second field. If we are at
end-of-file, the first element will be the C<Tracer::EOF> constant.

=back

=cut

sub _GetInputKVRecord {
    # Get the parameters.
    my ($handle) = @_;
    # Declare the return variables.
    my ($key, $value);
    # Read from the file.
    my $line = <$handle>;
    # Check to see if we got something.
    if (defined $line) {
        # Parse and return what we got. Note we strip the line terminator first.
        my $stripped = Tracer::Strip($line);
        ($key, $value) = split /\t/, $stripped, 2;
        # Insure the value is defined. If it is, we un-escape it.
        if (! defined $value) {
            $value = "";
        } else {
            $value = Tracer::UnEscape($value);
        }
    } else {
        # Here we've hit end-of-file, so we stuff in a trailer.
        ($key, $value) = (Tracer::EOF, "");
    }
    # Return the key and value.
    return ($key, $value);
}

=head3 PutOutputKVRecord

    FIG::_PutOutputKVRecord($handle, $key, $value);

Write a key-value pair to the output file. The value will automatically be
escaped. A tab will be used to separate the fields.

=over 4

=item handle

Open output file handle.

=item key

First field to put in the output record.

=item value

Value field to put in the output record. It will automatically be escaped. If it
is undefined, the method will have no effect. An undefined value therefore serves
as a deleted-line marker.

=back

=cut

sub _PutOutputKVRecord {
    # Get the parameters.
    my ($handle, $key, $value) = @_;
    # Only proceed if we have a value.
    if (defined $value) {
        # Escape the value.
        my $trueValue = Tracer::Escape($value);
        # Write the output record.
        print $handle "$key\t$trueValue\n";
    }
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

  if ($organism eq 'All') {
      return "$FIG_Config::global/Models/All/Scenarios";
  }
  else {
      return "$FIG_Config::organisms/$organism/Scenarios";
  }
}

=head3 get_scenario_info

    FIG->get_scenario_info(@subsystem_names)

Returns a reference to a hash containing the scenario information for the specified
subsystems.  The hash keys are subsystem names, the hash values are hashes keyed
by "scenarios" and "reactions".  "reactions" keys a hash of functional roles to lists
of reactions.  "scenarios" keys a hash of scenario names and with yet more hashes as values.
The keys to these hashes are the strings "input_compounds", "output_compounds", "map_ids",
"additional_reactions" and "ignore_reactions", values are references to lists of KEGG ids.
If a subsystem has no scenarios, no hash entry is created for that subsystem.

=over 4

=item @subsystem_names

A list of subsystem names.

=back

=cut

sub get_scenario_info{
  my ($self, @subsystem_names) = @_;

  my $return = {};

  foreach my $ssa (@subsystem_names)
  {
      $ssa =~ s/[ \/]/_/g;

      # using get_subsystem is really slow, and so we are going to cat the file and return that
      if (open(SSA,"<$FIG_Config::data/Subsystems/$ssa/hope_kegg_info")) {
        my @lines = <SSA>;
	chomp @lines;

	for (my $i = 0; $i < scalar @lines; $i += 6)
	{
	    my $scenario_name = $lines[$i];
	    my $input_compounds = $lines[$i+1];
	    my $output_compounds = $lines[$i+2];
	    my $map_ids = $lines[$i+3];
	    my $additional_reactions = $lines[$i+4];
	    my $ignore_reactions = $lines[$i+5];

	    $return->{$ssa}->{"scenarios"}->{$scenario_name}->{input_compounds} = [split(/,/,$input_compounds)];

	    my @output_compounds_lists;
	    my @inner_list;

	    foreach my $cpd (split(/,/,$output_compounds))
	    {
		if ($cpd =~ /\(/)
		{
		    $cpd =~ s/\(//g;
		    push @inner_list, $cpd;
		}
		elsif (scalar @inner_list > 0 && $cpd =~ /\)/)
		{
		    $cpd =~ s/\)//g;
		    push @inner_list, $cpd;
		    my @new_inner_list = @inner_list;
		    push @output_compounds_lists, \@new_inner_list;
		    @inner_list = ();
		}
		elsif (scalar @inner_list > 0)
		{
		    push @inner_list, $cpd;
		}
		else
		{
		    push @output_compounds_lists, [$cpd];
		}
	    }

	    $return->{$ssa}->{"scenarios"}->{$scenario_name}->{output_compounds} = \@output_compounds_lists;
	    $return->{$ssa}->{"scenarios"}->{$scenario_name}->{map_ids} = [split(/,/,$map_ids)];
	    $return->{$ssa}->{"scenarios"}->{$scenario_name}->{additional_reactions} = [split(/,/,$additional_reactions)];
	    $return->{$ssa}->{"scenarios"}->{$scenario_name}->{ignore_reactions} = [split(/,/,$ignore_reactions)];
	}
	close(SSA);
    }

      if (open(HR,"<$FIG_Config::data/Subsystems/$ssa/hope_reactions")) {
	  my @lines = <HR>;
	  chomp @lines;
	  foreach my $line (@lines) {
	      my ($fr, $reactions) = split "\t", $line;
	      $return->{$ssa}->{"reactions"}->{$fr} = [split(/,/,$reactions)];
	  }
	  close(HR);
      }
  }
  return $return;
}


#############
# Atomic regulon stuff. Temporary, should be replaced by
# Sapling calls.
#

=head3 in_atomic_regulon

    @membership = $fig->in_atomic_regulon($fid);

@membership is a list of triples [genome-id, regulon-num, regulon-size].

=cut

sub in_atomic_regulon
{
    my($self, $fid) = @_;

    my $res = $self->db_handle->SQL(qq(SELECT genome, regulon, size
				       FROM atomic_regulon
				       WHERE fid = ?), undef, $fid);
    return $res;
}

sub features_in_atomic_regulon
{
    my($self, $fid_list) = @_;

    my $dbh = $self->db_handle->{_dbh};
    return undef unless ref($fid_list) eq 'ARRAY';
    my $cond = join(",", map { $dbh->quote($_) } @$fid_list);
    my $res = $dbh->selectall_hashref(qq(SELECT fid, genome, regulon, size
					 FROM atomic_regulon
					 WHERE fid IN ($cond)), 'fid');
    return $res;
}


########################################################


#
# RabbitMQ message broker methods.
#

sub broker_connect
{
    my($self) = @_;
    my $conn = eval { Net::RabbitMQ->new(); };
    if ($@ =~ /Can\'t locate object/)
    {
	warn "RabbitMQ code not loaded";
	return;
    }
    my $sockfd = $conn->connect($FIG_Config::rabbitmq_host, {
	user => $FIG_Config::rabbitmq_user,
	password => $FIG_Config::rabbitmq_password,
    });
    if (!$sockfd)
    {
	warn "Cannot connect to RabbitMQ on $FIG_Config::rabbitmq_host";
	return;
    }

    my $chan = 1;
    $conn->channel_open($chan);

    $self->{_broker}->{connection} = $conn;
    $self->{_broker}->{channel} = $chan;
    $self->{_broker}->{sockfd} = $sockfd;

    my $exch = $FIG_Config::rabbitmq_exchange;

    $self->{_broker}->{exchange} = $exch;

    $conn->exchange_declare($chan, $exch, { exchange_type => "topic", durable => 1, auto_delete => 0 });

    $self->{_broker}->{connected} = 1;

    return 1;
}

sub broker_log
{
    my($self, $msg, $data) = @_;

    return if (!defined $FIG_Config::rabbitmq_host) || $FIG_Config::rabbitmq_host eq '';

    if (!$self->{_broker}->{connected})
    {
	#
	# If too many failures, don't even try.
	#
	return if ($self->{_broker}->{error_msgs} > 10);

	my $ok = eval { $self->broker_connect(); };
	if (!$ok || $@)
	{
	    if ($self->{_broker}->{error_msgs}++ < 10)
	    {
		carp "broker_log: could not connect to broker";
	    }
	    return;
	}
    }

    my $seed_id = $FIG_Config::seed_id || "unspecified";

    my $conn = $self->{_broker}->{connection};
    my $chan = $self->{_broker}->{channel};
    my $exch = $self->{_broker}->{exchange};

    my $enc_data = encode_json($data);
    $conn->publish($chan, "seed_log.$seed_id.$msg", $enc_data, { exchange => $exch },
	       { content_type => "application/json",
		 delivery_mode => 2, # persistent
		 });

}

sub get_uuid
{
    if ($haveUUID)
    {
	my($u, $str);
	UUID::generate($u);
	UUID::unparse($u, $str);
	return $str;
    }
    else
    {
	if (open(my $fh, "-|", "uuidgen"))
	{
	    my $str = <$fh>;
	    chomp $str;
	    close($fh);
	    return $str;
	}
	else
	{
	    return "INVALID";
	}
    }

}

sub is_refseq_id {
    my($id) = @_;

}


=head2 FIG::Job module

=cut

### Begin FIG::Job module

package FIG::Job;

use FIGAttributes;
use base 'FIGAttributes';

sub new
{
    my($class, $job_id, $job_dir) = @_;

    my $self = {
        id => $job_id,
        dir => $job_dir,
    };
    return bless $self, $class;
}

sub status :Scalar
{
    my($self) = @_;

    return &FIG::file_read("$self->{dir}/STATUS");
}

sub running :Scalar
{
    my($self) = @_;
    my $rc;
    warn "running test on $self->{id}\n";
    if (kill(0, $self->{id}) > 0)
    {
        $rc = 1;
    }
    else
    {
        $rc = 0;
    }
    warn "running returns $rc\n";

    return $rc;
}

sub info :Scalar :List
{
    my($self) = @_;
    return &FIG::file_read("$self->{dir}/INFO");
}

sub output :Scalar :List
{
    my($self) = @_;
    return &FIG::file_read("$self->{dir}/OUTPUT");
}



######### End FIG::Job ##

package FIG;

#
# DAS support.
#

=head3 init_das

Initialize a DAS data query object.

=cut

sub init_das
{
    my($self, $url, $dsn) = @_;

    my $das_data_dir = "$FIG_Config::var/DAS";

    if (-d $das_data_dir)
    {
        return new SeedDas($self,$das_data_dir, $url, $dsn);
    }
    else
    {
        return undef;
    }
}

package FIG::SimLock;

#
# Little package to implement a lock for sims work.
#

use strict;
use Fcntl qw/:flock/;  # import LOCK_* constants

sub new
{
    my($class) = @_;

    my $pool_dir = "$FIG_Config::global/sim_pools";
    &FIG::verify_dir($pool_dir);

    #
    # Lock the pool directory.
    #
    open(my $lock, ">$pool_dir/lockfile");

    flock($lock, LOCK_EX);

    my $self = {
        lock_fh => $lock,
    };

    return bless($self, $class);
}

sub DESTROY
{
    my($self) = @_;

    warn "$$ unlocking sims lock\n";
    $self->{lock_fh}->close();
}

package FIG;

{
    package GenomeDataCache;
    use strict;
    use Data::Dumper;
    use Carp;
    use Scalar::Util qw(weaken);

    sub new
    {
	my($class, $fig) = @_;

	my $self = {
	    fig => $fig,
	    cache => undef,
	};
	weaken($self->{fig});
	return bless $self, $class;
    }


    sub genomes
    {
	my($self, $complete, $restrictions, $domain) = @_;

	$self->load_cache();
	my @out;
	for my $ent (values %{$self->{cache}})
	{
	    next unless $ent->{genome};
	    next if $complete && !$ent->{complete};
	    next if $restrictions && !$ent->{restrictions};
	    next if $domain && $domain ne $ent->{maindomain};
	    push(@out, $ent->{genome});
	}

	return sort { FIG::by_genome_id($a,$b) } @out;
    }


    sub genome_list
    {
	my($self) = @_;

	$self->load_cache();
	my $out = [];
	for my $ent (values %{$self->{cache}})
	{
	    next unless $ent->{genome};
	    push(@$out, [$ent->{genome}, $ent->{gname}, $ent->{maindomain}]);
	}

    @$out = sort { lc $a->[1] cmp lc $b->[1] } @$out;
	return $out;
    }


    sub genome_info
    {
	my($self) = @_;

	$self->load_cache();
	my $out = [];
	for my $ent ( sort { FIG::by_genome_id( $a->{genome}, $b->{genome} ) }
	              grep { $_->{genome} }
	              values %{$self->{cache}}
	            )
	{
	    push(@$out, [ @$ent{qw(genome gname szdna maindomain pegs rnas complete taxonomy)} ]);
	}

    @$out = sort { lc $a->[1] cmp lc $b->[1] } @$out;
	return $out;
    }


    sub get_results_for_genome_counts
    {
	my($self, $complete) = @_;

	$self->load_cache();
	my $out = [];
	for my $ent ( values %{$self->{cache}} )
	{
	    next unless $ent->{genome};
	    next if $complete && !$ent->{complete};
	    push(@$out, [$ent->{genome}, $ent->{maindomain}]);
	}

	return $out;
    }


    sub AUTOLOAD
    {
	my($self, @params) = @_;
	our $AUTOLOAD;
	my $name = $AUTOLOAD;
	$name =~ s/.*:://;

	if (defined($self->{cache}))
	{
	    confess "Method not found: $name";
	}

	$self->load_cache();

	return $self->$name(@params);

    }

    sub DESTROY
    {
	my($self) = @_;
	print STDERR "Destroy cache $self $$\n";
    }

    sub load_cache
    {
	my($self) = @_;

	return if defined($self->{cache});

	# print STDERR "Loading genome cache\n";

	my $sth= $self->{fig}->db_handle->dbh->prepare("SELECT * FROM genome");
	$sth->execute;

	for my $field (@{$sth->{NAME}})
	{
	    no strict 'refs';
	    *{$field} = sub {
		my($self, $id) = @_;

		#
		# If we have forked, the AUTOLOAD will not trigger in the child,
		# so we must check for the cache existence here.
		#
		if (!$self->{cache})
		{
		    # print STDERR "Loading genome cache in child\n";
		    $self->load_cache();
		}
		return $self->{cache}->{$id}->{$field};
	    };
	}
	$self->{cache} = $sth->fetchall_hashref('genome');
    }
}

{
    package SubsystemDataCache;
    use Carp;
    use Data::Dumper;
    use Scalar::Util qw(weaken);
    sub new
    {
	my($class, $fig) = @_;

	my $self = {
	    fig => $fig,
	    cache => undef,
	};
	weaken($self->{fig});
	return bless $self, $class;
    }

    sub cache_ok
    {
	my($self) = @_;
	return 0 if !defined($self->{cache});
	my @s = stat($FIG_Config::var . "/ss_modified");
	return 0 unless @s;
	return ($s[9] < $self->{cache_time});
    }

    sub flush_cache
    {
	my($self) = @_;
	$self->{cache} = undef;
    }

    sub all_subsystems
    {
	my($self, $only_usable) = @_;

	$self->load_cache();

	my @out;
	for my $ent (values %{$self->{cache}})
	{
	    next if ($only_usable && ! $self->{fig}->usable_subsystem_classification([$ent->{class_1}, $ent->{class_2}]));
	    push(@out, $ent->{subsystem});
	}

	return @out;
    }

    sub all_subsystems_detailed
    {
	my($self, $only_usable) = @_;

	$self->load_cache();

	my @out;
	for my $ent (values %{$self->{cache}})
	{
	    next if ($only_usable && ! $self->{fig}->usable_subsystem_classification([$ent->{class_1}, $ent->{class_2}]));
	    push(@out, { %$ent } );
	}

	return @out;
    }

    sub subsystem_classification
    {
	my($self, $ssa) = @_;
	$self->load_cache();
	my $c1 = $self->class_1($ssa);
	my $c2 = $self->class_2($ssa);
	return [$c1, $c2];
    }

    sub AUTOLOAD
    {
	my($self, @params) = @_;
	our $AUTOLOAD;
	my $name = $AUTOLOAD;
	$name =~ s/.*:://;

	if (defined($self->{cache}))
	{
	    confess "Method not found: $name";
	}

	$self->load_cache();

	return $self->$name(@params);

    }

    sub DESTROY
    {
	my($self) = @_;
#	print "Destroy cache $self\n";
    }

    sub load_cache
    {
	my($self) = @_;
	return if defined($self->{cache}) && $self->cache_ok();

	$self->{cache_time} = time;
	my $expt = "";
	if ($FIG_Config::exclude_experimental_subsystems)
	{
	    $expt = " WHERE class_1 NOT LIKE 'Experimental%'";
	}
	my $sth= $self->{fig}->db_handle->dbh->prepare("SELECT * FROM subsystem_metadata $expt");
	$sth->execute;

	for my $field (@{$sth->{NAME}})
	{
	    no strict 'refs';
	    *{$field} = sub {
		my($self, $id) = @_;
		#
		# If we have forked, the AUTOLOAD will not trigger in the child,
		# so we must check for the cache existence here.
		#
		if (!$self->{cache} || !$self->cache_ok())
		{
		    $self->load_cache();
		}
		return $self->{cache}->{$id}->{$field};
	    };
	}
	$self->{cache} = $sth->fetchall_hashref('subsystem');

	$_->{curator} =~ s/^master:// for values %{$self->{cache}};
    }
}

1;
