#!/usr/bin/perl

use strict;
use warnings;

use DBD::mysql;
use Digest::MD5::File qw(file_md5_hex);
use File::stat;
use Data::Dumper;

=head1 NAME

bindex.pl - gathers metadata about bio-files into a database, now bam files

=head1 About

Yeah.

=cut

main();

sub main {

    my $dbh = DBI->connect( 'DBI:mysql:biofiles', 'biofiles', '' )
      || die "Error connection to biofiles: $@";

# We have no state, so we could easily reconnect, but for now half hour should be enough.
    my $sth = $dbh->prepare('SET wait_timeout = 1800')
      || die "prepare: $dbh->errstr()";
    $sth->execute() || die "execute: $dbh->errstr()";
    $sth->finish();

# Use mlocate to find files of interest; report symlinks as files, check that file still exists and is readable and use regex.
    my @bamfiles = `/usr/bin/locate -Per '\\.bam\$'`;

    # Retrieve file identification data for found file.
    my $sth = $dbh->prepare(
'SELECT md5sum, filesize, unix_timestamp(filetime) FROM bamfile WHERE filename = ?'
    ) || die "prepare: $dbh->errstr()";

    # Update filesize and timestamp of a file in db.
    my $sth2 = $dbh->prepare(
'UPDATE bamfile SET filesize = ?, filetime = from_unixtime(?) WHERE md5sum = ?'
    ) || die "prepare: $dbh->errstr()";

    foreach my $bamfile (@bamfiles) {
        chomp($bamfile);

        # What if it no longer exists?
        next if ( -d $bamfile );
        next if ( not -e $bamfile );

        #    print "Processing $bamfile\n";
        my $bam = Bamfile->new( dbh => $dbh, filename => $bamfile );

# It's possible to load in pre-calculated MD5SUM's for example with:
# LOAD DATA INFILE '/tmp/md5sum.lst' IGNORE INTO TABLE bamfile COLUMNS TERMINATED BY '  ' (md5sum,filename);
# In which case this code will fill in the filesize & time.

        # Retrieve file identification data for found file.
        $sth->execute($bamfile) || die "execute: $dbh->errstr()";

        # If file exists in database, make sure it hasn't changed.
        my ( $md5sum, $filesize, $filetime, $filename );
        if ( ( $md5sum, $filesize, $filetime ) = $sth->fetchrow() ) {
            if ( $bam->changed ) {
                if (   ( 0 == $filesize )
                    || ( 0 == $filetime ) )
                {
                    # Get file properties for updating them in the database.
                    my $fstat = stat($bamfile) || die "File $bamfile stat: $@";
                    $sth2->execute( $fstat->size, $fstat->mtime, $md5sum )
                      || die "execute: $dbh->errstr()";
                }
            }
        }
    }

    $sth2->finish();
    $sth->finish();

    $dbh->disconnect();
}

package Bamfile;

use warnings;
use strict;

use DBD::mysql;
use Digest::MD5::File qw(file_md5_hex);
use File::stat;
use Data::Dumper;
use Cwd qw(abs_path);

=head1 NAME

Bamfile - Binary Alignment Map from sequencing

=head1 About

This object is used for getting and tracking metadata about BAM files in the system.

=head1 Constructor

=head2 new

  my $bamfile = Bamfile->new(dbh => $dbh, filename => 'bamfile.bam');

=over

=item dbh

DBI DBD::mysql database handle with database set up.

=back

=cut

sub new {
    my $class   = shift;
    my %options = @_;
    my $self    = {
        dbh      => undef,
        filename => undef,
        md5sum   => undef,
        filesize => undef,
        filetime => undef,
        %options,
    };
    bless( $self, $class );

    $self->load();

    return ($self);
}

=head1 Attributes

=head2 dbh

Database handle

=head2 filename

Absolute path in filesystem

=head2 filesize

size of complete file in bytes

=head2 filetime

UNIX epoch timestamp of file modification time

=head2 md5sum

md5sum of whole file in lowe-case hex; if it doesn't exist in database it's calculated

=cut

sub md5sum {
    my $self = shift;

    my $md5sum   = $self->{md5sum};
    my $filename = $self->{filename};

    while ( -l $filename ) { $filename = abs_path($filename) }

    if ( ( !defined $md5sum ) && ( defined $filename ) ) {
        $md5sum = file_md5_hex($filename);
    }

    return $md5sum;
}

=head2 changed

Indicates whether the filesystem file matches one in database.

=cut

sub changed {
    my $self = shift;

    my $filename = $self->{filename};
    my $filetime = $self->{filetime};
    my $filesize = $self->{filesize};

    # Get file properties for determining if it's still the same size.
    my $fstat = stat($filename) || die "File stat: $@";

    if ( ( $filesize != $fstat->size ) || ( $filetime != $fstat->mtime ) ) {
        print(  "different $filename: time $filetime new "
              . $fstat->mtime
              . ", size $filesize new "
              . $fstat->size
              . "\n" );
        return -1;
    }
    return 0;

}

=head1 Methods

=head2 load

Attempts to load identifying information from database.

=cut

sub load {
    my $self = shift;

    my $dbh      = $self->{dbh};
    my $filename = $self->{filename};

    if ( defined $filename ) {

        # Retrieve file identification data for the file, considering aliases

        while ( -l $filename ) { $filename = abs_path($filename) }

        my $sth = $dbh->prepare(
            <<END
SELECT b.md5sum, b.filesize, unix_timestamp(b.filetime)
    FROM alias a
        RIGHT OUTER JOIN bamfile b
            ON a.md5sum = b.md5sum
    WHERE ? in (a.filename, b.filename)
END
        ) || die "prepare: $dbh->errstr()";
        $sth->execute($filename) || die "execute: $dbh->errstr()";

# If file path doesn't exist in database, check if it exists elsewhere by same md5sum.
        if (
            !(
                ( $self->{md5sum}, $self->{filesize}, $self->{filetime} ) =
                $sth->fetchrow()
            )
          )
        {
            print "Checking new file: $self->{filename}\n";
            $self->{md5sum} = $self->md5sum;
            print "$self->{md5sum}  $self->{filename}\n";

            my $sth2 = $dbh->prepare(
'SELECT filesize, unix_timestamp(filetime), filename FROM bamfile WHERE md5sum = ?'
            ) || die "prepare: $dbh->errstr()";
            $sth2->execute( $self->{md5sum} ) || die "execute: $dbh->errstr()";
            if ( my ( $filesize, $filetime, $filename ) = $sth2->fetchrow() ) {
                print(
"Found by md5sum: filename $filename, filesize $filesize, filetime $filetime\n"
                );
                $self->{filesize} = $filesize;
                $self->{filetime} = $filetime;
            }
            else {

                # It doesn't exist in database, so for now add it.
                my $fstat = stat( $self->{filename} )
                  || die "File $filename stat: $@";
                $self->{filesize} = $fstat->size;
                $self->{filetime} = $fstat->mtime;

                $self->store();
            }
            $sth2->finish();
        }
        $sth->finish();

  # Update the path alias whether it existed or not, so we can spot stale names.
        $self->update_alias();

        return -1;
    }

    return 0;
}

=head2 store_alias

Store current object into the file name alias list.

=cut

sub update_alias {
    my $self = shift;

    my $dbh      = $self->{dbh};
    my $filename = $self->{filename};
    my $md5sum   = $self->{md5sum};
    my $filetime = $self->{filetime};

    my $realfile;
    if ( -l $self->{filename} ) {
        $realfile = $filename;
        while ( -l $realfile ) { $realfile = abs_path($realfile) }
    }

    # Profiling shows REPLACE and ON DUPLICATE KEY are as fast in this case.
    my $sth =
      $dbh->prepare(
'INSERT INTO alias (md5sum, filename, filetime, symlink) VALUES ( ?, ?, from_unixtime(?), ? ) ON DUPLICATE KEY UPDATE modified=NOW()'
      ) || die "prepare: $dbh->errstr()";
    $sth->execute( $md5sum, $filename, $filetime, $realfile )
      || die "execute: $dbh->errstr()";
    $sth->finish();

    return -1;
}

=head2 store

Store the bamfile object into database; currently insert only to notice errors.

=cut

sub store {
    my $self = shift;

    my $dbh      = $self->{dbh};
    my $filename = $self->{filename};
    my $md5sum   = $self->{md5sum};
    my $filesize = $self->{filesize};
    my $filetime = $self->{filetime};

    my $sth = $dbh->prepare(
'INSERT INTO bamfile ( md5sum, filesize, filetime, filename) VALUES ( ?, ?, FROM_UNIXTIME(?), ?)'
    ) || die "prepare: $dbh->errstr()";
    $sth->execute( $md5sum, $filesize, $filetime, $filename )
      || die "execute: $dbh->errstr()";
    $sth->finish();

    return -1;
}

1;
