#! /bin/bash
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# Install SExtractor for the Dockerized build.

set -xeuo pipefail

version=2.28.0

cd /tmp
curl -fsSL "https://github.com/astromatic/sextractor/archive/refs/tags/${version}.tar.gz" |tar xz
pushd sextractor-${version}
patch -p1 -i ../sextractor-fixes.patch
./autogen.sh
./configure --prefix=/dasch --with-cfitsio-libdir=/dasch/lib --with-cfitsio-incdir=/dasch/include LIBS="-lz"
make
make install
popd
rm -rf sextractor-${version}
