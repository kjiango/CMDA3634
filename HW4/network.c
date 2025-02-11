#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "omp.h"

typedef struct {

  // Start your code
  int nodes;
  int maxConnects;
  int totalConnects;
  double d;
  int* connectCount;
  int* source;
  int* dest;
  int* offset;
  double* oldPageRank;
  double* newPageRank;
  // End your code

}Network;

typedef struct {
  int dest;
  int source;

}Link;


int destCompare(const void * u, const void * v) {
  Link *linkU = (Link*)u;
  Link *linkV = (Link*)v;

  return linkU->dest - linkV->dest;
}


Network networkReader(char* filename) {

  // Start your code
  Network g;
  int i, j, temp;
  FILE *fp;
  int source;
  int dest;
  int count=0;

  fp = fopen(filename, "r");
  fscanf(fp, "%d,", &g.nodes);
  temp = fscanf(fp, "%d,", &g.maxConnects);

  g.connectCount = (int*)calloc(g.nodes, sizeof(int)); 
  
  int* sourceArray = (int*)calloc(g.nodes*g.maxConnects, sizeof(int));
  int* destArray = (int*)calloc(g.nodes*g.maxConnects, sizeof(int));

  g.oldPageRank = (double*)calloc(g.nodes, sizeof(double));
  g.newPageRank = (double*)calloc(g.nodes, sizeof(double));
 
  // Check if fscanf has read to the End of File
  do{
    temp = fscanf(fp, "%d,%d,", &source, &dest);     
    if(temp == EOF) break;

    sourceArray[count] = source;
    destArray[count] = dest;

    g.connectCount[sourceArray[count]]++;
    count++;
  }while(temp != EOF);
  
  g.totalConnects = count;
  fclose(fp);
  
  Link links[g.totalConnects];
  
  for(i = 0; i < g.totalConnects; i++) {
    links[i].source = sourceArray[i];
    links[i].dest = destArray[i];
  }
  free(sourceArray);
  free(destArray);

  qsort(&links, g.totalConnects, sizeof(Link), destCompare);
  g.source = (int*)calloc(g.totalConnects, sizeof(int));
  g.dest = (int*)calloc(g.totalConnects, sizeof(int));
  g.offset = (int*)calloc(g.nodes+1, sizeof(int));
 
  for(i = 0; i < g.totalConnects; i++) {
    g.source[i] = links[i].source;
    g.dest[i] = links[i].dest;
    g.offset[g.dest[i]+1]++;
  }
  for(i = 1; i < g.nodes+1; i++) {
    g.offset[i] += g.offset[i-1];
  }
  return g;
  // End your code

}


/**
Computes the difference norm between the old and new page ranks
*/
double computeDiff(Network net) {
  int i; 
  double diff = 0; 

  for(i = 0; i < net.nodes; i++) {
    diff = pow(net.oldPageRank[i] - net.newPageRank[i], 2);
  }
  return sqrt(diff);
}


/**
Updates the page rank values using the PR algorithm
*/
double updatePageRank(Network net) {
  int j, i;  
  int source;
  int dest;
  double L;
  double PR;

#pragma omp parallel for
  for(i = 0; i < net.nodes; i++) {
    for(j = net.offset[i]; j < net.offset[i+1]; j++) {
      dest = net.dest[j];
      source = net.source[j];   
      L = (double)net.connectCount[source];
      PR = net.oldPageRank[source];

      net.newPageRank[dest] += net.d*(PR/L);
    }
  }
  return computeDiff(net);
}


/**
Computes the page rank values for each node
*/
void computePageRank(Network net) {
  
  double tol = 1e-6;
  double diff = 1;
  double* temp;
  int k;
  double nodes = (double)net.nodes;
  double constant = (1-net.d)/nodes;

  /* Initialize old PR vlaues */
  for(k = 0; k < net.nodes; k++) {
    net.oldPageRank[k] = 1/nodes;
  }
  while(diff > tol) {
    /* Set all new PR values to constant */
  #pragma omp parallel for reduction(+:diff)
    for(k = 0; k < net.nodes; k++) {
      net.newPageRank[k] = constant;
    }

    /* Update PR based on connections */
    diff = updatePageRank(net);

    /* Switch pointers updated values in oldPR array
       Gets us ready to run again */
    temp = net.oldPageRank;
    net.oldPageRank = net.newPageRank;
    net.newPageRank = temp;
  }

}


/**
Frees memory in the Network struct
*/
void networkDestructor(Network net) {

  free(net.connectCount);
  free(net.source);
  free(net.dest);
  free(net.offset);
  free(net.oldPageRank);
  free(net.newPageRank);

}


/**
Main method 
*/
int main(int argc, char** argv) {

  // Reads the filename of the data file
  char* filename = argv[1];
  int n_threads = 4;

  omp_set_num_threads(n_threads);

  //Start your code
  Network net = networkReader(filename);

  int j;
  int k = 0;
  net.d = 0.85;

  computePageRank(net);

  for(j = 0; j < net.nodes; j++) {
    printf("PageRank of node %d is: %f\n", j, net.oldPageRank[j]);
  }
  networkDestructor(net);
  
  return 0;

}
