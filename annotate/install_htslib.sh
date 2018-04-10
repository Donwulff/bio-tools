#!/bin/bash
# Installs prerequisites, gets and runs Thomas Krahn's BigY2 (hg38) annotation script
# Tested on Windows 10 with Windows Store installed WSL Ubuntu
# This script assumes you have only one BigY VCF zip in your Windows download folder(s)
sudo apt update
sudo apt upgrade
sudo apt-get install unzip
# Recommended by linuxbrew, don't think it's needed for this
#sudo apt-get install build-essential
sed -i "/\/home\/linuxbrew\/\.linuxbrew/d" ~/.profile
. ~/.profile
sh -c "$(curl -fsSL https://raw.githubusercontent.com/Linuxbrew/install/master/install.sh)"
echo 'export PATH="/home/linuxbrew/.linuxbrew/bin:$PATH"' >>~/.profile
echo 'export MANPATH="/home/linuxbrew/.linuxbrew/share/man:$MANPATH"' >>~/.profile
echo 'export INFOPATH="/home/linuxbrew/.linuxbrew/share/info:$INFOPATH"' >>~/.profile
. ~/.profile
brew install htslib
brew install samtools
brew install bcftools
wget -nc https://gist.github.com/tkrahn/283462028c61cd213399ba7f6b773893/raw/38f9e5c2448247e27d35023fd09e7c7923b9000b/bigY_hg38_pipeline.sh
sed -i "s/^wget http:/wget -N http:/g" bigY_hg38_pipeline.sh
sed -i "/snps_hg38.vcf.gz.tbi$/a \
# The ; in middle of INFO data confuses bcftools 1.8, swap the downloaded files around and remove the char\n\
mv -f snps_hg38.vcf.gz snps_hg38.vcf.gz.bak\n\
mv -f snps_hg38.vcf.gz.tbi snps_hg38.vcf.gz.tbi.bak\n\
\n\
zcat snps_hg38.vcf.gz.bak | sed 's/R1b-_;L1066/R1b-_L1066/g' | bgzip -c > snps_hg38.vcf.gz\n\
tabix -f snps_hg38.vcf.gz\n\
\n\
touch snps_hg38.vcf.gz -f snps_hg38.vcf.gz.bak\n\
touch snps_hg38.vcf.gz.tbi -f snps_hg38.vcf.gz.tbi.bak\n\
echo \"modifying annotation file complete\"" bigY_hg38_pipeline.sh
chmod a+x bigY_hg38_pipeline.sh
unzip /mnt/c/Users/*/Downloads/bigy2-*.zip
./bigY_hg38_pipeline.sh
