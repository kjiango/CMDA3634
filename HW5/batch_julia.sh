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
nvcc -O3 -o julia -arch=sm_60 cudaJulia.cu png_util.c -I. -lpng -lm
#
./julia 4096 4096 64
#
rm julia
#
echo "Exited script normally
