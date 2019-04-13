#!/bin/sh

# Converts Yoruba Sequence (L3e2b1a1) to rCRS (H2a2a1)
# Some DNA analysis use UCSC hg19 human reference, which uses Yoruban mtDNA sequence.
# Coordinate conversion between the two is easy enough, but some of the reference alleles are also different.
# The 'diff' array holds differences from rCRS to the Yoruban haplogroup.
# If the sample has a call for a reverse of one of the Yoruban variants (Ie. towards rCRS/H2a2a1), skip the variant.
# If the sample has no call for one of those locations, print the Yoruban reference variant instead.
# This way we remove variants which match the rCRS reference, and add missing ones which match the Yoruban reference.
# Note that at least some of the early vcf files didn't have calls for all locations of mtDNA,
# in which case spurious Yoruban variants will be added instead!

# Reference: https://www.mitomap.org/MITOMAP/MitoSeqs
# Convert and annotate online: https://mseqdr.org/mvtool.php
# Variant haplogroup analysis: https://dna.jameslick.com/mthap/

# Use tabix if it exists, creating index if needed.
PARSE="zgrep '^chrM[[:space:]]' $1"
if which tabix > /dev/null;
then
  if [ ! -e $1.tbi ];
  then
    if tabix $1 2> /dev/null;
    then
      PARSE="tabix $1 chrM"
    fi
  else
    PARSE="tabix $1 chrM"
  fi
fi

eval $PARSE | awk -v ORS=' ' ' \
BEGIN { \
  split("A73G C150T T195C A263G \
         A750G A1438G T2352C T2483C A2706G A4769G T5580C C7028T A8701G A8860G A9377G \
         T9540C A10398G A10819G T10873C G11719A C12705T T14212C C14766T G14905A G15301A A15326G \
         T16172C T16189C C16223T C16320T", diff); \
} \
/^chrM[[:space:]]/ \
{ \
  if( $2>=16192 ) $2-=2; \
    else if( $2>=3109 ) $2-=1; \
    else if( $2>=318 ) $2-=2; \
    else if( $2>=311 ) $2-=1; \
  for(snp in diff) \
    if ($5$2$4==diff[snp]) { \
      delete diff[snp]; \
      next; \
    } \
  print "*"$2$5; \
} \
END { \
  for(snp in diff) print substr(diff[snp],2) \
}' | sort -n
