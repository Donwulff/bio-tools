#!/bin/bash

set -e

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
export CXXFLAGS=$CFLAGS
export LDFLAGS=$CFLAGS
export LIBS="../zlib/libz.a $CFLAGS"
export AR="gcc-ar"
export RANLIB="gcc-ranlib"

# Update package list
sudo apt update
sudo apt upgrade -y

# Install necessary packages
sudo apt install -y build-essential autoconf

# Install additional package for libdeflate
sudo apt install -y cmake

git_clone_or_pull https://github.com/ebiggers/libdeflate.git libdeflate
cd libdeflate
cmake -DCMAKE_BUILD_TYPE="Release" -DCMAKE_C_FLAGS="$CFLAGS"
make -j8
sudo -E make install
cd ..

# Build zlib
git_clone_or_pull https://github.com/cloudflare/zlib.git zlib
cd zlib
./configure
make -j8
set +e
sudo -E make install
set -e
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
sudo -E make install
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
sudo -E make install
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
sudo -E make install
cd ..

# Build bedtools
git_clone_or_pull https://github.com/arq5x/bedtools2.git bedtools2
cd bedtools2
sed -i 's/ -lz / /g' Makefile
sed -i 's/^BT_CXXFLAGS = /BT_CXXFLAGS := /' Makefile
sed -i 's/^BT_LIBS    = /BT_LIBS := /' Makefile
sed -i "/^BT_CXXFLAGS :=/ s/$/ \${CPPFLAGS}/" Makefile
sed -i "/^BT_LIBS :=/ s/$/ \${LIBS}/" Makefile
make -j8
sudo -E make install
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
sudo cp bwa /usr/local/bin/
cd ..

# Install additional package for isa-l
sudo apt install -y libtool nasm

# Build isa-l
git_clone_or_pull https://github.com/intel/isa-l.git isa-l
cd isa-l
./autogen.sh
./configure
make -j8
sudo -E make install
cd ..

# Build fastp - use v0.22.0 or this until the memleak is fixed
git_clone_or_pull https://github.com/Donwulff/fastp fastp
cd fastp
make -j8
sudo -E make install
cd ..

# Get and install precompiled k8-0.2.5
wget https://github.com/attractivechaos/k8/releases/download/0.2.5/k8-0.2.5.tar.bz2
tar jxf k8-0.2.5.tar.bz2
sudo cp k8-0.2.5/k8-Linux /usr/local/bin

exit # bwa-postalts.js doesn't look compatible with 1.0, huge memory leak.

# Build nodejs
wget -O- https://nodejs.org/dist/v18.17.0/node-v18.17.0.tar.gz | tar -zxf -
cd node-v18.17.0
export CFLAGS="-march=native -O3"
export CPPFLAGS=$CFLAGS
export CXXFLAGS=$CFLAGS
export LDFLAGS=$CFLAGS
./configure
make -j8
# Build k8 under nodejs
git_clone_or_pull https://github.com/attractivechaos/k8 k8
cd k8
export CFLAGS="-march=native -flto=8 -O3"
export CPPFLAGS=$CFLAGS
export CXXFLAGS=$CFLAGS
export LDFLAGS=$CFLAGS
sed -i 's/^CFLAGS= /CFLAGS := /' Makefile
sed -i 's/^LIBS= /LIBS := /' Makefile
sed -i "/^CFLAGS :=/ s/$/ \${CFLAGS}/" Makefile
sed -i "/^LIBS :=/ s/$/ \${LIBS}/" Makefile
make
sudo cp k8 /usr/local/bin
cd ..
cd ..
