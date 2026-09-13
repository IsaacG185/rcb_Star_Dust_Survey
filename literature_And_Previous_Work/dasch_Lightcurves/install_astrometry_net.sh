#! /bin/bash
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# Install Astrometry.Net for the Dockerized build.

set -xeuo pipefail

version=0.93

cd /tmp
curl -fsSL "http://astrometry.net/downloads/astrometry.net-${version}.tar.gz" |tar xz
pushd astrometry.net-${version}
make CFITS_INC="-I/dasch/include" CFITS_LIB="/dasch/lib/libcfitsio.a -lz"
make install INSTALL_DIR=/dasch CFITS_INC="-I/dasch/include" CFITS_LIB="/dasch/lib/libcfitsio.a -lz"
popd
rm -rf astrometry.net-${version}
