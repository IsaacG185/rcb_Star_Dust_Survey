#! /bin/bash
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# Install Starbase for the Dockerized build.

commit=815abc2bfc9481d2c1ad8d7f9d0ef04baa2467f3

set -xeuo pipefail

cd /tmp
curl -fsSL "https://github.com/jbroll/starbase/archive/$commit.tar.gz" |tar xz
mv "starbase-$commit" starbase
pushd starbase
patch -p1 -i ../starbase-fixes.patch

./configure --prefix=/dasch
make
make install

pushd tawk
./configure --prefix=/dasch
make
make install
popd

popd
rm -rf starbase
