#! /bin/bash
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# Install SCAMP for the Dockerized build.

set -xeuo pipefail

version=2.10.0

cd /tmp
curl -fsSL "https://github.com/astromatic/scamp/archive/refs/tags/v${version}.tar.gz" |tar xz
pushd scamp-${version}
patch -p1 -i ../scamp-fixes.patch
./autogen.sh
./configure --prefix=/dasch
make
make install
popd
rm -rf scamp-${version}
