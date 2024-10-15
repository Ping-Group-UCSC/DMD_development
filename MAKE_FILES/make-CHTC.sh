#!/bin/bash

module purge
module load openmpi/5.0.3-gcc-13.2.0 
module load cmake
module load gcc/13.2.0
module list

make -j4
