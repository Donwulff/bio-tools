#!/usr/bin/perl

use strict;
use warnings;

# Generates list of minimal changes between two dbSNP versions
# Contigs are listed without prefix and should be in same order in input files
# + new SNP, = changes SNP, - removed SNP

open( my $old, '<', $ARGV[0] ) || die "Opening $ARGV[0]: $!\n";
open( my $new, '<', $ARGV[1] ) || die "Opening $ARGV[1]: $!\n";

my ( $oldrow, $newrow, $oldchr );

# Skip header, which according to spec can only be at beginning. We don't care about content right now.
while ( defined( $oldrow = <$old> ) && ( substr( $oldrow, 0, 1 ) eq '#' ) ) { }
while ( defined( $newrow = <$new> ) && ( substr( $newrow, 0, 1 ) eq '#' ) ) { }
MAIN: while ( defined $oldrow || defined $newrow ) {
    if ( !defined $oldrow ) {
        my @newcols1 = split( "\t", $newrow );
        if ( !defined($oldchr) || !( $oldchr eq $newcols1[0] ) ) {
            print("$newcols1[0]\n");
            $oldchr = $newcols1[0];
        }
        print "+" . join( "\t", @newcols1[1..4] ) . "\n";
        $newrow = <$new>;
        next;
    }
    if ( !defined $newrow ) {
        my @oldcols1 = split( "\t", $oldrow );
        if ( !defined($oldchr) || !( $oldchr eq $oldcols1[0] ) ) {
            print("$oldcols1[0]\n");
            $oldchr = $oldcols1[0];
        }
        print "-" . $oldcols1[2] . "\n";
        $oldrow = <$old>;
        next;
    }

    # Process single genomic location at a time.
    my @oldrows  = ($oldrow);
    my @oldcols1 = split( "\t", $oldrow );
    my @newrows  = ($newrow);
    my @newcols1 = split( "\t", $newrow );

    if ( !defined($oldchr) || !( $oldchr eq $newcols1[0] ) ) {
        print("$newcols1[0]\n");
        $oldchr = $newcols1[0];
    }

    if ( !( $oldcols1[0] eq $newcols1[0] ) || ( $oldcols1[1] <= $newcols1[1] ) )
    {
        $oldrow = <$old>;
        my @oldcols = split( "\t", $oldrow // '' );
        while ( @oldcols && ( $oldcols1[0] eq $oldcols[0] )
            && ( $oldcols1[1] eq $oldcols[1] ) )
        {
            push( @oldrows, $oldrow );
            $oldrow = <$old>;
            @oldcols = split( "\t", $oldrow // '' );
        }

        @oldrows = sort {
            substr( ( split( "\t", $a ) )[2], 2 ) <=>
              substr( ( split( "\t", $b ) )[2], 2 )
        } @oldrows;
        if (   ( $oldcols1[0] eq $newcols1[0] )
            && ( $oldcols1[1] < $newcols1[1] ) )
        {
            print map { "-" . ( split( "\t", $_ ) )[2] . "\n" } @oldrows;
        }
    }

    if ( !( $oldcols1[0] eq $newcols1[0] ) || ( $oldcols1[1] >= $newcols1[1] ) )
    {
        $newrow = <$new>;
        my @newcols = split( "\t", $newrow // '' );
        while ( @newcols && ( $newcols1[0] eq $newcols[0] )
            && ( $newcols1[1] eq $newcols[1] ) )
        {
            push( @newrows, $newrow );
            $newrow = <$new>;
            @newcols = split( "\t", $newrow // '' );
        }

        @newrows = sort {
            substr( ( split( "\t", $a ) )[2], 2 ) <=>
              substr( ( split( "\t", $b ) )[2], 2 )
        } @newrows;
        if (   ( $oldcols1[0] eq $newcols1[0] )
            && ( $oldcols1[1] > $newcols1[1] ) )
        {
            print
              map { "+" . join( "\t", ( split( "\t", $_ ) )[ 1 .. 4 ] ) . "\n" }
              @newrows;
        }

    }

 # If genomic location matches, we need to make sure the variants match as well.
    if (   ( $oldcols1[0] eq $newcols1[0] )
        && ( $oldcols1[1] == $newcols1[1] ) )
    {
        my ( $oldrnum, $newrnum ) = ( 0, 0 );
        while (1) {
            if ( !defined( $oldrows[$oldrnum] ) ) {
                while ( defined( $newrows[$newrnum] ) ) {
                    print "+"
                      . join( "\t",
                        ( split( "\t", $newrows[$newrnum] ) )[ 1 .. 4 ] )
                      . "\n";
                    $newrnum++;
                }
                last;
            }
            if ( !defined( $newrows[$newrnum] ) ) {
                while ( defined( $oldrows[$oldrnum] ) ) {
                    my @oldcols = split( "\t", $oldrows[$oldrnum] );
                    print "-$oldcols[2]\n";
                    $oldrnum++;
                }
                last;
            }
            my @oldcols = split( "\t", $oldrows[$oldrnum] );
            my @newcols = split( "\t", $newrows[$newrnum] );
            if ( $oldcols[2] eq $newcols[2] ) {
                if (   !( $oldcols[3] eq $newcols[3] )
                    || !( $oldcols[4] eq $newcols[4] ) )
                {
                    print "="
                      . join( "\t",
                        ( split( "\t", $newrows[$newrnum] ) )[ 1 .. 4 ] )
                      . "\n";
                }
                $oldrnum++;
                $newrnum++;
                next;
            }
            else {
                if ( substr( $oldcols[2], 2 ) < substr( $newcols[2], 2 ) ) {
                    print "+"
                      . join( "\t",
                        ( split( "\t", $newrows[$newrnum] ) )[ 1 .. 4 ] )
                      . "\n";
                    $newrnum++;
                }
                else {
                    print "-$oldcols[2]\n";
                    $oldrnum++;
                }
            }
        }
    }
}
