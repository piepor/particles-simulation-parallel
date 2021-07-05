#!/bin/sh
rm *.ppm *.jpg *.dmp *.sta
rm particles_f.exe
module load gnu
rm *.exe
gfortran -o particles_f.exe -pg -O3 -fbounds-check particles_f.f90
./particles_f.exe
