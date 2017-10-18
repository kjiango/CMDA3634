#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Q2a: add OpenMP header include here
#include "omp.h"

int main(int argc, char **argv){

  long long int test = 0;
  long long int Ninside = 0; // number of random points inside 1/4 circle

  double newPi = 0, estPi = 0;
  double tol = 1e-8;

  // Q2b: add OpenMP API code to set number of threads to 10 here
  int Nthreads = 10;
  omp_set_num_threads(Nthreads); // start with Nthreads threads
  
  struct drand48_data *drandData; 
  drandData = (struct drand48_data*) malloc(Nthreads*sizeof(struct drand48_data));

  // Q2c: add an OpenMP parallel region here, wherein each thread initializes 
  //      one entry in drandData using srand48_r and seed based on thread number
#pragma omp parallel {
  long int seed = 0;
  srand48_r(seed, drandData+0);
}
  
  do{
    int Ninnertests=10000;

    estPi = newPi;
    // Q2d: add an OpenMP parallel for directive here
    //      to split this loop amongst the 10 threads
    //      [ add a reduction clauss for the Ninside variable ]
    for(int n=0;n<Ninnertests;++n){
      double x, y;
      ++test;
      
      // call threadsafe reentrant version of drand48
      drand48_r(drandData+0, &x);
      drand48_r(drandData+0, &y);
      
      if(x*x+y*y<1){
        ++Ninside;
      }
    }

    newPi = Ninside/(double)test;
    printf("newPi = %lf\n", 4.*newPi);
  }while(fabs(newPi-estPi)>tol);

  printf("\n");
  printf("estPi = %lf\n", 4.*newPi);

  free(drandData);

  return 0;
}
