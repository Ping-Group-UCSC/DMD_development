#!/bin/bash

mkdir build
cp -r src_v4.5.7 ./build/src_v4.5.7
cp Makefile ./build/
cp make.inc ./build/
cd ./build
mkdir bin

module purge
module load openmpi/5.0.3-gcc-13.2.0 
module load cmake
module load gcc/13.2.0
module list

make -j4
