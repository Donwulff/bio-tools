#!/bin/bash

# Function to clone or pull a git repository
git_clone_or_pull() {
    REPO_URL=$1
    TARGET_DIR=$2

    if [ -d "$TARGET_DIR" ]; then
        echo "Directory $TARGET_DIR exists. Pulling latest changes..."
        cd $TARGET_DIR
        git pull
        cd ..
    else
        echo "Cloning repository $REPO_URL into $TARGET_DIR..."
        git clone $REPO_URL $TARGET_DIR
    fi
}

# Set environment variables
export CFLAGS="-march=native -flto=8 -O3"
export CPPFLAGS=$CFLAGS
export LDFLAGS=$CFLAGS
export LIBS="../zlib/libz.a $CFLAGS"
export AR="gcc-ar"
export RANLIB="gcc-ranlib"

# Update package list
sudo apt update

# Install necessary packages
sudo apt install -y build-essential autoconf

# Build zlib
git_clone_or_pull https://github.com/cloudflare/zlib.git zlib
cd zlib
./configure
make -j8
sudo make install
cd ..

# Install additional packages for htslib
sudo apt install -y libbz2-dev liblzma-dev

# Build htslib
git_clone_or_pull https://github.com/samtools/htslib.git htslib
cd htslib
git submodule update --init --recursive
autoreconf -i
./configure
sed -i 's/ -lz / /g' Makefile *.mk
make -j8
sudo make install
cd ..

# Install additional package for samtools
sudo apt install -y libncurses-dev

# Build samtools
git_clone_or_pull https://github.com/samtools/samtools.git samtools
cd samtools
autoreconf -i
./configure
sed -i 's/ -lz / /g' Makefile
make -j8
sudo make install
cd ..

# Install additional package for bcftools
sudo apt install -y libgsl-dev

# Build bcftools
git_clone_or_pull https://github.com/samtools/bcftools.git bcftools
cd bcftools
autoreconf -i
sed -i 's/ -lz / /g' Makefile
./configure --enable-libgsl
make -j8
sudo make install
cd ..

# Build bwa
git_clone_or_pull https://github.com/lh3/bwa.git bwa
cd bwa
# Modify Makefile to remove -lz, then append external CFLAGS and LIBS
sed -i 's/-lz//g' Makefile
sed -i 's/^CFLAGS= /CFLAGS := /' Makefile
sed -i 's/^LIBS= /LIBS := /' Makefile
sed -i "/^CFLAGS :=/ s/$/ \${CFLAGS}/" Makefile
sed -i "/^LIBS :=/ s/$/ \${LIBS}/" Makefile
make -j8
sudo cp -a bwa /usr/local/bin/
cd ..

