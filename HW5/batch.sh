#! /bin/bash
#
#PBS -l walltime=00:05:00
#PBS -l nodes=1:ppn=1:gpus=1
#PBS -W group_list=newriver
#PBS -q p100_dev_q
#PBS -A CMDA3634
#
cd $PBS_O_WORKDIR
#
module purge
module load cuda
#
nvcc -O3 -o mandelbrot -arch=sm_60 mandelbrot.cu png_util.c -I. -lpng -lm
#
./mandelbrot 4096 4096 32
./mandelbrot 4096 4096 64
./mandelbrot 4096 4096 128
./mandelbrot 4096 4096 256
./mandelbrot 4096 4096 512
#
rm mandelbrot
#
echo "Exited script normally"
