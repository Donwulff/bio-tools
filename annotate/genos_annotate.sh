#!/bin/bash
# Genos Research provided exome VCF's have weird format which cointains NT=Not Targeted Regions and NC=No Call
# This script is intended to merge thhe Genos Research VCF with the variants in dbSNP to fill in ref-calls
# Different bcftools have had different behavior for merging SNP's, so this needs to be tested
# Two different merges are attempted, one for multiallelic and one for separate calls
# After merging, output is annotated with CADD and MCAP pathogenity values
# This script has some additional bits to try to check concordance against DNA.Land imputed genomes to see if we got merging right
set -x
# Set or give on command line; Genos sample is assumed to be gVCF without REF lines.
SAMPLE=$1
DBSNP=/mnt/GenomicData/All_20170710_GRCh37p13.vcf.gz
REF=/mnt/FDA/hs37d5.fa
# http://krishna.gs.washington.edu/download/CADD/v1.3/whole_genome_SNVs_inclAnno.tsv.gz; 344GB so not just auto-downloading it.
CADD1=/mnt/CADD_v1.3/prescored/whole_genome_SNVs_inclAnno.tsv.gz
CADD2=/mnt/CADD_v1.3/prescored/InDels_inclAnno.tsv.gz
CADD3=/mnt/CADD_v1.3/local.tab.gz
# Should be SAMPLE, but it isn't now.
SAMPLENAME=${SAMPLE}

/bin/mkdir -p data

if [ ! -e ${SAMPLE}_chrnum.vcf.gz ];
then
  echo "## Renaming chromosomes to numeric..."
  bcftools annotate --rename-chrs /mnt/GenomicData/chr_to_num $SAMPLE.vcf.gz \
    | bgzip -c > ${SAMPLE}_chrnum.vcf.gz
#  tabix -p vcf ${SAMPLE}_chrnum.vcf.gz
fi

if [ ! -e ${SAMPLE}_multi.vcf.gz.tbi ];
then
  bcftools norm -cw -f ${REF} -m -any ${SAMPLE}_chrnum.vcf.gz | bgzip -c > ${SAMPLE}_norm.vcf.gz
  tabix -p vcf ${SAMPLE}_norm.vcf.gz
  bcftools norm -cw -f ${REF} -m +any ${SAMPLE}_chrnum.vcf.gz | bgzip -c > ${SAMPLE}_multi.vcf.gz
  tabix -p vcf ${SAMPLE}_multi.vcf.gz
fi

if [ ! -e ${SAMPLE}_multi_dbSNP.vcf.gz.tbi ];
then
  echo "## Annotating with dbSNP annotations..."
  tabix -p vcf ${DBSNP}
  if [ ! -e ${SAMPLE}_multi_dbSNP.vcf.gz ];
  then
    /usr/bin/bcftools annotate -a $DBSNP -c ID,^INFO/RS,^INFO/RSPOS,^INFO/RV ${SAMPLE}_norm.vcf.gz \
      | bgzip -c > ${SAMPLE}_dbSNP.vcf.gz
    /usr/bin/bcftools annotate -a $DBSNP -c ID,^INFO/RS,^INFO/RSPOS,^INFO/RV ${SAMPLE}_multi.vcf.gz \
      | bgzip -c > ${SAMPLE}_multi_dbSNP.vcf.gz
  fi
  LASTPOS=`/bin/zgrep "^1[[:space:]]" ${SAMPLE}_dbSNP.vcf.gz | tail -1 | cut -f2`
  if [ $LASTPOS -lt 249212624 ];
  then
    echo "-- bcftools has a bug causing it to cut chromosomes, try to install working version"
    exit
  fi
  tabix -p vcf ${SAMPLE}_dbSNP.vcf.gz
  tabix -p vcf ${SAMPLE}_multi_dbSNP.vcf.gz
fi

if [ ! -e data/All_GRCh37_reference.vcf.gz ];
then
  echo "## Building a fake reference-only sample..."
  bcftools norm -cw -f ${REF} $DBSNP | sed '
s/INFO$/INFO\tFORMAT\tsample-id/
/^#/!s/$/\tGT:GQ:DP\t0\/0:0:0/
8i##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
' | bgzip -c > data/All_GRCh37_reference.vcf.gz
fi;

if [ ! -e data/All_GRCh37_reference_nomult.vcf.gz ];
then
  echo "## Building a fake reference-only sample..."
  bcftools norm -cw -f ${REF} -m -any $DBSNP | sed '
s/INFO$/INFO\tFORMAT\tsample-id/
/^#/!s/$/\tGT:GQ:DP\t0\/0:0:0/
8i##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
' | bgzip -c > data/All_GRCh37_reference_nomult.vcf.gz
fi;

if [ ! -e data/All_GRCh37_reference_multi.vcf.gz ];
then
  echo "## Building a fake reference-only sample..."
  bcftools norm -cw -f ${REF} -m +any $DBSNP | sed '
s/INFO$/INFO\tFORMAT\tsample-id/
/^#/!s/$/\tGT:GQ:DP\t0\/0:0:0/
8i##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
' | bgzip -c > data/All_GRCh37_reference_multi.vcf.gz
fi;

# Perhaps check timestamp? Check broken on purpose now.
if [ ! -e data/All_GRCh37_${SAMPLE}_multi.vcf.gz.tbi.not ];
then
  echo "## Renaming sample in reference file..."
  echo "$SAMPLENAME" > ${SAMPLE}_samplename
  bcftools reheader -s ${SAMPLE}_samplename data/All_GRCh37_reference.vcf.gz -o data/All_GRCh37_${SAMPLE}.vcf.gz
  tabix -p vcf data/All_GRCh37_${SAMPLE}.vcf.gz
  bcftools reheader -s ${SAMPLE}_samplename data/All_GRCh37_reference_nomult.vcf.gz -o data/All_GRCh37_${SAMPLE}_nomult.vcf.gz
  tabix -p vcf data/All_GRCh37_${SAMPLE}_nomult.vcf.gz
  bcftools reheader -s ${SAMPLE}_samplename data/All_GRCh37_reference_multi.vcf.gz -o data/All_GRCh37_${SAMPLE}_multi.vcf.gz
  tabix -p vcf data/All_GRCh37_${SAMPLE}_multi.vcf.gz
fi

# Documentation on bcftools concat -D looks wrong, instead of -d none that does't even exist it's -d all
if [ ! -e ${SAMPLE}_nomult_dbSNP_concat.vcf.gz ];
then
  echo "## Expanding REF and NC blocks... (non-annotated version)"
  /usr/bin/bcftools concat -D -a ${SAMPLE}_dbSNP.vcf.gz data/All_GRCh37_${SAMPLE}_nomult.vcf.gz | gawk -e '
  BEGIN { OFS="\t"; }
  /END=/ { chrom = $1; split( $8, info, ";" ); end = int( substr( info[1], 5 ) ); filter = $7; }
  !/^#/ { if( $1==chrom && $2>end ) print; if ( $1==chrom && $2<=end && filter=="NC" && $7!="NC" ) { $7 = filter; $10 = "./."; print; } }
  /^#/ { print; }' | bgzip -c > ${SAMPLE}_nomult_dbSNP_concat.vcf.gz
fi

if [ ! -e ${SAMPLE}_dbSNP_concat.vcf.gz ];
then
  echo "## Expanding REF and NC blocks... (non-annotated version)"
  /usr/bin/bcftools concat -D -a ${SAMPLE}_dbSNP.vcf.gz data/All_GRCh37_${SAMPLE}.vcf.gz | gawk -e '
  BEGIN { OFS="\t"; }
  /END=/ { chrom = $1; split( $8, info, ";" ); end = int( substr( info[1], 5 ) ); filter = $7; }
  !/^#/ { if( $1==chrom && $2>end ) print; if ( $1==chrom && $2<=end && filter=="NC" && $7!="NC" ) { $7 = filter; $10 = "./."; print; } }
  /^#/ { print; }' | bgzip -c > ${SAMPLE}_dbSNP_concat.vcf.gz
fi

if [ ! -e ${SAMPLE}_multi_dbSNP_concat.vcf.gz ];
then
  echo "## Expanding REF and NC blocks... (non-annotated multiallelic version)"
  /usr/bin/bcftools concat -D -a ${SAMPLE}_multi_dbSNP.vcf.gz data/All_GRCh37_${SAMPLE}_multi.vcf.gz | gawk -e '
  BEGIN { OFS="\t"; }
  /END=/ { chrom = $1; split( $8, info, ";" ); end = int( substr( info[1], 5 ) ); filter = $7; }
  !/^#/ { if( $1==chrom && $2>end ) print; if ( $1==chrom && $2<=end && filter=="NC" && $7!="NC" ) { $7 = filter; $10 = "./."; print; } }
  /^#/ { print; }' | bgzip -c > ${SAMPLE}_multi_dbSNP_concat.vcf.gz
fi

# Now we drop the multiallelic version...
if [ ! -e ${SAMPLE}_dbSNP_CADD1.3.vcf.gz ];
then
  echo "## Annotating with CADD annotations, this'll take a while..."
  # Construct headers for CADD run; generic way left just in case
  /bin/zcat ${CADD1} | head -3 | grep "^#Chrom" | cut -c2- | tr "\\t" "\\n" | grep -Ev "^(Chrom|Pos|Ref|Alt)$" \
    | gawk -e '{ print "##INFO=<ID="$1",Number=1,Type=String,Description=\"CADD\">" }' > data/CADD_INFO_lines.vcf
  # Requires matching columns; bool = Flag, num = Float, int = Integer can't be handled by bcftools...
  cat data/CADD_columns.txt | grep -Ev "^(Chrom|Pos|Ref|Alt)$" | gawk -e '
  {
    gsub( /[\(\)]/, "", $2 );
    switch($3) {
      case "bool":   $3="String"; break;
      case "string":
      case "factor": $3="String"; break;
      case "num":    $3="String"; break;
      case "int":    $3="String";
    };
    desc=substr( $0, index( $0, " " )+1);
    desc=substr( desc, index( desc, " " )+1);
    desc=substr( desc, index( desc, " " )+1);
    gsub( /"/, "", desc );
    print "##INFO=<ID="$2",Number=1,Type="$3",Description=\""desc"\">"
  }' > data/CADD_INFO_lines.vcf

  COLUMNS=`zcat ${CADD1} | head -3 | grep "^#Chrom" | cut -c2- | tr "\\t" ,`
  /bin/zcat ${SAMPLE}_dbSNP.vcf.gz \
    | bcftools annotate -a ${CADD1} -h data/CADD_INFO_lines.vcf -c ${COLUMNS} -Ou \
    | bcftools annotate -a ${CADD2} -h data/CADD_INFO_lines.vcf -c ${COLUMNS} -Ou \
    | bcftools annotate -a ${CADD3} -h data/CADD_INFO_lines.vcf -c ${COLUMNS} \
    | bgzip -c > ${SAMPLE}_dbSNP_CADD1.3.vcf.gz
fi

if [ ! -e data/mcap_v1_0.tab.gz.tbi ];
then
  echo "## Constructing MCAP 1.0 annotation..."
  wget http://bejerano.stanford.edu/MCAP/downloads/dat/mcap_v1_0.txt.gz -o data
  zcat data/mcap_v1_0.txt.gz | bgzip -c > data/mcap_v1_0.tab.gz
  tabix -b2 -e2 data/mcap_v1_0.tab.gz
  echo '##INFO=<ID=MCAP,Number=1,Type=Float,Description="M-CAP">' > data/mcap_v1_0.tab.gz.head
fi

if [ ! -e ${SAMPLE}_dbSNP_CADD1.3_MCAP1.0.vcf.gz ];
then
  bcftools annotate -a data/mcap_v1_0.tab.gz -h data/mcap_v1_0.tab.gz.head -c "chrom,pos,ref,alt,INFO/MCAP" ${SAMPLE}_dbSNP_CADD1.3.vcf.gz \
    | bgzip -c > ${SAMPLE}_dbSNP_CADD1.3_MCAP1.0.vcf.gz
  tabix -p vcf ${SAMPLE}_dbSNP_CADD1.3_MCAP1.0.vcf.gz
fi;

# New version uses "-d all" and -D doesn't seem to work, but cuts the sample file off at ~200k as well!
if [ ! -e ${SAMPLE}_dbSNP_CADD1.3_MCAP1.0_concat.vcf.gz ];
then
  echo "## Expanding REF and NC blocks..."
  /usr/bin/bcftools concat -D -a ${SAMPLE}_dbSNP_CADD1.3_MCAP1.0.vcf.gz data/All_GRCh37_${SAMPLE}.vcf.gz | gawk -e '
  BEGIN { OFS="\t"; }
  /END=/ { chrom = $1; split( $8, info, ";" ); end = int( substr( info[1], 5 ) ); filter = $7; }
  !/^#/ { if( $1==chrom && $2>end ) print; if ( $1==chrom && $2<=end && filter=="NC" && $7!="NC" ) { $7 = filter; $10 = "./."; print; } }
  /^#/ { print; }' | bgzip -c > ${SAMPLE}_dbSNP_CADD1.3_MCAP1.0_concat.vcf.gz
fi

if [ ! -e ${SAMPLE}_dbSNP_CADD1.3_MCAP1.0_concat.types ];
then
  echo "## Report of site types found in sample:"
  /bin/zgrep -v -E "(^#|^Y|^M)" ${SAMPLE}_dbSNP_CADD1.3_MCAP1.0_concat.vcf.gz | gawk -e '
  {
    if(substr($3,1,2)==".") novel="Novel"; else novel="dbSNP";
    if($7==".") $7="REF"; print novel " " $7
  }' | sort | uniq -c | tee ${SAMPLE}_dbSNP_CADD1.3_MCAP1.0_concat.types
fi

if [ ! -e ${SAMPLE}_nonc.vcf.gz ];
then
  echo "## Creating clean calls, PASS (non-REF), REF and all non-NC subsets"
  /bin/zgrep -v -E "[[:space:]](NC|Indel_QC|SNP_QC)[[:space:]]" ${SAMPLE}_dbSNP_CADD1.3_MCAP1.0_concat.vcf.gz | bgzip -c > ${SAMPLE}_clean.vcf.gz
  /bin/zgrep -E "(^#|[[:space:]]PASS[[:space:]])" ${SAMPLE}_dbSNP_CADD1.3_MCAP1.0_concat.vcf.gz | bgzip -c > ${SAMPLE}_pass.vcf.gz
  /bin/zgrep -v -E "[[:space:]](PASS|NC|Indel_QC|SNP_QC)[[:space:]]" ${SAMPLE}_dbSNP_CADD1.3_MCAP1.0_concat.vcf.gz | bgzip -c > ${SAMPLE}_refs.vcf.gz
  /bin/zgrep -v "[[:space:]]NC[[:space:]]" ${SAMPLE}_dbSNP_CADD1.3_MCAP1.0_concat.vcf.gz | bgzip -c > ${SAMPLE}_nonc.vcf.gz
fi
  vcftools --gzvcf ${SAMPLE}_clean.vcf.gz --gzdiff ../DNALand/dnl.imported.vcf.gz --diff-indv-discordance --diff-indv-map data/genos_dnl.map --out ${SAMPLE}_clean
  vcftools --gzvcf ${SAMPLE}_pass.vcf.gz --gzdiff ../DNALand/dnl.imported.vcf.gz --diff-indv-discordance --diff-indv-map data/genos_dnl.map --out ${SAMPLE}_pass
  vcftools --gzvcf ${SAMPLE}_refs.vcf.gz --gzdiff ../DNALand/dnl.imported.vcf.gz --diff-indv-discordance --diff-indv-map data/genos_dnl.map --out ${SAMPLE}_refs
  vcftools --gzvcf ${SAMPLE}_nonc.vcf.gz --gzdiff ../DNALand/dnl.imported.vcf.gz --diff-indv-discordance --diff-indv-map data/genos_dnl.map --out ${SAMPLE}_nonc
  vcftools --gzvcf ${SAMPLE}_clean.vcf.gz --gzdiff ../DNALand/dnl.imported.vcf.gz --diff-site-discordance --diff-indv-map data/genos_dnl.map --out ${SAMPLE}_clean
  vcftools --gzvcf ${SAMPLE}_pass.vcf.gz --gzdiff ../DNALand/dnl.imported.vcf.gz --diff-site-discordance --diff-indv-map data/genos_dnl.map --out ${SAMPLE}_pass
  vcftools --gzvcf ${SAMPLE}_refs.vcf.gz --gzdiff ../DNALand/dnl.imported.vcf.gz --diff-site-discordance --diff-indv-map data/genos_dnl.map --out ${SAMPLE}_refs
  vcftools --gzvcf ${SAMPLE}_nonc.vcf.gz --gzdiff ../DNALand/dnl.imported.vcf.gz --diff-site-discordance --diff-indv-map data/genos_dnl.map --out ${SAMPLE}_nonc

if [ ! -e ${SAMPLE}.23andMe ];
then
  echo "## Producing 23andMe style output..."
  cp /mnt/GenomicData/23andMe.header ${SAMPLE}.23andMe
  /bin/zgrep -v "^Y" ${SAMPLE}_clean.vcf.gz | /usr/bin/gawk -e '
    BEGIN { OFS = "\t"; }
    !/^#/ { a[1] = substr($10,1,1)=="0"?$4:(substr($10,1,1)=="1"?$5:"-");
            if(length($10)>=3) { a[2] = substr($10,3,1)=="0"?$4:(substr($10,3,1)=="1"?$5:"-") };
            asort(a);
            if(length($4) == 1 && length($5) == 1 && $3 != ".")
                print $3,$1,$2,a[1]a[2] }
  ' | /usr/bin/uniq >> ${SAMPLE}.23andMe
fi

#vcftools --gzvcf annotated_dbSNP_CADD1.3_MCAP1.0_concat.vcf.gz --gzdiff DNALand/dnl.imported.vcf.gz --diff-indv-discordance --diff-indv-map genos_dnl.map


#vcftools --diff-indv-map genos_dnl.map --diff-site-discordance --gzvcf annotated_dbSNP_CADD1.3_concat.vcf.gz --gzdiff DNALand/dnl.imported.vcf.gz

#grep -v -- "-nan" out.diff.sites | less | grep -v "0$" | cut -f1,2 | tail -n+2 | xargs -IX bash -c "zcat annotated_nonc.vcf.gz DNALand/dnl87_frj.imported.vcf.gz | grep '^X[[:space:]]' ; echo ----" > discordant.txt

#grep -v -- "-nan" sample_clean.diff.sites | grep -v "0$" | cut -f1,2 | tail -n+2 | xargs -IX bash -c "echo --- X ---; zcat 13211510242132_dbSNP_CADD1.3_MCAP1.0_concat.vcf.gz ../DNALand/dnl87_frj.imported.vcf.gz | grep '^X[[:space:]]'"

#grep -v -- "-nan" sample_clean.diff.sites | grep -v "0$" | cut -f1,2 | tail -n+2 | xargs -IX bash -c "echo --- X ---; zcat 13211510242132_clean.vcf.gz ../DNALand/dnl87_frj.imported.vcf.gz | grep '^X[[:space:]]'" > clean.txt
