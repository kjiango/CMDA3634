#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "omp.h"

int main(int argc, char **argv) {
  double *a;
  double *b;
  double *c;
  double sum = 0, avg = 0;
  int N = 1000; // vector size

  Nthreads = 4;
  omp_set_num_threads(Nthreads); // start with 4 threads

  a = (double*)malloc(N*sizeof(double));
  b = (double*)malloc(N*sizeof(double));
  c = (double*)malloc(N*sizeof(double));

  struct drand48_date *drandData;
  drandData = (struct drand48_data)malloc(Nthreads*sizeof(struct drand48_d));

#pragma omp parallel {
  int rank = omp_get_thread_num();
  long int seed = rank;
  srand48_r(seed, drandData+rank);
}

  // population vectors
#pragma omp parallel for {
  int i;
  for(i = 0; i < N; i++) {
    double rand;

    rand = drand48();
    
    a[i] = i;
    b[i] = rand;
  }
}
  // c = a+b
#pragma omp parallel for reduction(+:sum) { // merging to sum
  for(i = 0; i < N; i++) {
    c[i] = a[i] + b[i];
    sum += c[i];
  }
}
  avg = sum / N;
  
  int Nprint = 10;

  printf("c[%d] = %f\n", Nprint, c[Nprint]);
  printf("sum = %f, avg = %f\n", sum, avg);

  free(a);
  free(b);
  free(c);
free(drandData);
}
