#!/usr/bin/env bash
set -eu

pushd .

cd build1D/
make
cd ~1


cd out/

mpirun -n 4 ../build1D/amrvac -i ../conf1D.nml ../first.par
mpirun -n 4 ../build1D/amrvac -i ../conf1D.nml ../second.par
