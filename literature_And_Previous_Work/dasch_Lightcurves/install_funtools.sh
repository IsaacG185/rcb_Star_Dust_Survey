#! /bin/bash
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# Install Funtools for the Dockerized build.

set -xeuo pipefail

version=1.4.8

cd /tmp
curl -fsSL "https://github.com/ericmandel/funtools/archive/refs/tags/v${version}.tar.gz" |tar xz
pushd funtools-${version}
./mkconfigure
./configure --prefix=/dasch
make
make install
popd
rm -rf funtools-${version}
