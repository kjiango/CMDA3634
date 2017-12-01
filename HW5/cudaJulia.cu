/* 

To compile:

   gcc -O3 -o mandelbrot mandelbrot.c png_util.c -I. -lpng -lm

To create an image with 4096 x 4096 pixels (last argument will be used to set number of threads):

    ./mandelbrot 4096 4096 1

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "png_util.h"

// Q2a: add include for CUDA header file here:
#include "cuda.h"

#define MXITER 1000


typedef struct {
  
  double r; // real
  double i; // imaginary
  
}complex_t;


// return iterations before z leaves mandelbrot set for given c
__device__ int testpoint(complex_t c, complex_t z){

  int iter;
  double temp;

  z = c;
  
  for(iter = 0; iter < MXITER; iter++){  
    temp = (z.r*z.r) - (z.i*z.i) + c.r;
    
    z.i = z.r*z.i*2. + c.i;
    z.r = temp;
    
    if((z.r*z.r+z.i*z.i) > 4.0){
      return iter;
    }
  }
  return iter; 
}


// perform Julia iteration on a grid of numbers in the complex plane
// record the  iteration counts in the count array
__global__ void julia(int Nre, int Nim, complex_t zmin, complex_t dz, complex_t c, float *count){ 

  // Q2c: replace this loop with a CUDA kernel
  complex_t z;
  int thread = threadIdx.x;
  int block = blockIdx.x;
  int blockSize = blockDim.x;
  int id = block*blockSize + thread;

  int m = id%Nre; // real axis
  int n = id%Nim; // imag axis

  z.r = zmin.r + dz.r*m;
  z.i = zmin.i + dz.i*n;
     
  count[m + n*Nre] = (float) testpoint(c, z);
}


/**
Main method
*/
int main(int argc, char **argv){

  // to create a 4096x4096 pixel image [ last argument is placeholder for number of threads ] 
  // usage: ./mandelbrot 4096 4096 32 

  int Nre = atoi(argv[1]);
  int Nim = atoi(argv[2]);
  int Nthreads = atoi(argv[3]);

  // Q2b: set the number of threads per block and the number of blocks here:
  int Nblocks = (Nre*Nim + Nthreads-1)/Nthreads;

  // storage for the iteration counts
  float *count;
  float *device_count;
  count = (float*)malloc(Nre*Nim*sizeof(float));
  cudaMalloc(&device_count, Nre*Nim*sizeof(float));

  // Parameters for a bounding box for "c" that generates an interesting image
  const float centRe = -1.6, centIm = 0.312;
  const float diam = 3.14;

  complex_t zmin; 
  complex_t zmax;
  complex_t dz;
  complex_t c;

  zmin.r = centRe - 0.5*diam;
  zmax.r = centRe + 0.5*diam;
  zmin.i = centIm - 0.5*diam;
  zmax.i = centIm + 0.5*diam;

  //set step sizes
  dz.r = (zmax.r-zmin.r)/(Nre-1);
  dz.i = (zmax.i-zmin.i)/(Nim-1);

  c.i = 0.1560;
  c.r = -0.8;

  clock_t start = clock(); //start time in CPU cycles

  // compute julia set
  julia <<<Nthreads, Nblocks>>> (Nre, Nim, zmin, dz, c, count); 
  
  // copy from the GPU back to the host here
  cudaMemcpy(count, device_count, Nre*Nim*sizeof(float), cudaMemcpyDeviceToHost);

  clock_t end = clock(); //start time in CPU cycles
  
  // print elapsed time
  printf("elapsed = %f\n", ((double)(end-start))/CLOCKS_PER_SEC);


  // output julia to png format image
  FILE *fp = fopen("julia.png", "w");

  printf("Printing julia.png...");
  write_hot_png(fp, Nre, Nim, count, 0, 80);
  printf("done.\n");

  free(count);

  exit(0);
  return 0;
}  