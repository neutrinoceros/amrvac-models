#!/usr/bin/env bash


pushd .

mkdir build
cd build
ln -s ../mod_usr.t mod_usr.t
$AMRVAC_DIR/setup.pl -d 1 -arch default

set -eu

make -j 4
cd ~1

cd out/
cp ../build/amrvac .
mpirun -n 4 amrvac -i ../conf1D.nml ../first.par
mpirun -n 4 amrvac -i ../conf1D.nml ../second.par
