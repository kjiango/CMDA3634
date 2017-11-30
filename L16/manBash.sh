#! /bin/bash
#
#PBS -l walltime=00:05:00
#PBS -l nodes=1:ppn=12
#PBS -W group_list=hokiespeed
#PBS -q normal_q
#PBS -j oe
cd $PBS_O_WORKDIR

module purge
module load gcc
gcc -O3 -o mandelbrot mandelbrot.c png_util.c -I. -lpng -lm -fopenmp

for n in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20}
do
./mandelbrot 4096 4096 $n 
done
