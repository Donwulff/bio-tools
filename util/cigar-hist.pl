#!/usr/bin/perl
use strict;
use warnings;

# Experimental quality control histogram, checking CIGAR and MD tags.
# May want to run it like:
# samtools view -F 0xF00 bwa-mapped-duplicate-marked.bam | perl cigar-hist.pm > mapped-duplicate-marked.bam.cigar

my @histogram;
my @trimming;
my @trimcount;
my @deletion;
my @insertion;

my @nucleotides;
my @variants;
my @varlength;

my @exact;
my @exactnum;

while (<STDIN>) {
    my @col = split("\t");
    my ( $cigar, $md ) = ( $col[5], grep { /^MD:Z:/ } @col );
    my ( $maxexact, $numexact, $nucu, $nucc, $match, $trimmed, $trimc ) =
      ( 0, 0, 0, 0, 0, 0, 0 );
    if ( defined($md) ) {
        $md = substr( $md, 4 );
        while ( $md =~ /([0-9]+)/g ) {
            my $count = $1;
            $maxexact = $count if $count > $maxexact;
            $numexact++;
            if ( $md =~ /[\^]?([^0-9]*)/g ) {
                my $chr = $1;
                my $len = length($chr);
                if ($len) {
                    $nucc++;
                    $nucu += $len;
                    $varlength[$len]++;
                }
            }
        }
        $exact[$maxexact]++;
        $exactnum[$numexact]++;
        $nucleotides[$nucu]++;
        $variants[$nucc]++;
    }
    while ( $cigar =~ /([0-9]+)/g ) {
        my $count = $1;
        if ( $cigar =~ /(.)/g ) {
            my $chr = $1;
            if ( 'M' eq $chr ) {
                $match = $count if $count > $match;
            }
            elsif ( 'S' eq $chr or 'H' eq $chr ) {
                $trimmed += $count;
                $trimc++;
            }

            # Insertion and deletion can only be inside matches.
            elsif ( 'D' eq $chr ) {
                $deletion[$count]++;
            }
            elsif ( 'I' eq $chr ) {
                $insertion[$count]++;
            }
            else { print }
        }
    }
    $histogram[$match]++;
    $trimming[$trimmed]++;
    $trimcount[$trimc]++;
}

# Not pass by reference for simplicity
sub histogram {
    my @histogram = @_;
    foreach ( sort { $a <=> $b } keys @histogram ) {
        print "$_ : $histogram[$_]\n" if defined $histogram[$_];
    }
}

print("Lengths of reads matching reference\n");
histogram(@histogram);
print("Amount of trimming per read\n");
histogram(@trimming);
print("Number of ends trimmed per read\n");
histogram(@trimcount);
print("Lengths of deletions\n");
histogram(@deletion);
print("Lengths of insertions\n");
histogram(@insertion);

print("Edit distance per read\n");
histogram(@nucleotides);
print("Number of variants per read\n");
histogram(@variants);
print("Length of variants\n");
histogram(@varlength);

print("Maximum exact match per read\n");
histogram(@exact);
print("Number of matching segments per read\n");
histogram(@exactnum);
