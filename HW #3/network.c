#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mpi.h"

typedef struct{
  // Network Properties

  // Number of Nodes in whole network
  int nodes;

  // Maximum number of connections each node can have
  int maxConnects;

  // Number of connections stored in struct
  int totalConnects;

  // Number of nodes this processor is responsible for
  int locNodes;

  // Relaxation constant
  double d;

  // Array pointers

  // Number of outbound connections array pointer
  int* connectCount;
  // Sources array pointers
  int* source;
  // Destinations array pointers
  int* dest;

  // Book keeping for send and recieve buffers
  int* sendCount; // Size of send
  int* recvCount; // nuber of recvs the process is to expect
  int* sendOffset; // Managing one array for all outgoing data
  int* recvOffset; // Managing one array for all incoming data
  int* recvLoc; // Companion of recvBuffer. Contins local node destinantion of received data

  int totalSends; // Total sends this processor needs to make
  int totalRecvs; // Total recvs this processor has to make

  // PageRank array pointers
  double* oldPageRank;
  double* newPageRank;

}Network;

/**
Network reader with parallel computing
 */
Network mpi_networkReader(char* filename){
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  Network g;
  int i, j, temp;
  FILE *fp;
  int source;
  int dest;

  /*locNodes = g.nodes/size;
  if((nodes%size) - 1 >= rank) {
    locNodes++;   // determines number of nodes per processor
    } */
  if(rank == 0){
    /* If root, read from data and distribute to all processors.
       Each processor should have global variables like number of nodes/max connects,
       every connection where one of its nodes are involved
       and the number of outbound connections of each of its nodes. See assignment
       for splitting up node ownership and local indexing information.
    */
    int count = 0;

    fp = fopen(filename, "r");
    fscanf(fp, "%d,", &g.nodes);
    temp = fscanf(fp, "%d,", &g.maxConnects);
    //MPI_Bcast(&g.nodes,1,MPI_INT,0,MPI_COMM_WORLD);
    //MPI_Bcast(&g.maxConnects,1,MPI_INT,0,MPI_COMM_WORLD);
    
    g.connectCount = (int*)calloc(g.nodes, sizeof(int));
    g.source = (int*)calloc(g.nodes*g.maxConnects, sizeof(int));
    g.dest = (int*)calloc(g.nodes*g.maxConnects, sizeof(int));
    int *proc_source = (int*)calloc(g.nodes*g.maxConnects, sizeof(int));
    int *proc_dest = (int*)calloc(g.nodes*g.maxConnects, sizeof(int));
    // double *sendBuffer = 

    g.sendOffset = (int*)malloc(sizeof(int));
    g.recvOffset = (int*)malloc(sizeof(int));

    do{
      temp = fscanf(fp, "%d,%d,", &source, &dest);

      if(temp == EOF) break;

      g.source[count] = source;
      g.dest[count] = dest;
      g.connectCount[g.source[count]]++;

      count++;
    }while(temp != EOF);

    for(i = 0; i < g.nodes; i++) {
      int *sent = (int*)malloc(sizeof(int));
      sent[0] = g.connectCount[i];
      MPI_Send(sent,1,MPI_INT,i%size,1,MPI_COMM_WORLD);
    }
    for(i = 0; i < count; i++) {
	proc_source[g.source[i]%size]++;
	proc_dest[g.dest[i]%size]++;
    }
    //MPI_Bcast(proc_source,size,MPI_INT,0,MPI_COMM_WORLD);
    //MPI_Bcast(proc_dest,size,MPI_INT,0,MPI_COMM_WORLD);

    /*
    for(i = 0; i < size; i++) {
      int *sourceBuffer = (int*)malloc(proc_source[i] + sizeof(int));
      int *destBuffer = (int*)malloc(proc_dest[i] + sizeof(int));

      int temp1 = 0;
      int temp2 = 0;

      } */
    fclose(fp);

    if(rank != 0) {
      MPI_Send(g.source, count, MPI_INT, rank, 1, MPI_COMM_WORLD);
      MPI_Send(g.dest, count, MPI_INT, rank, 2, MPI_COMM_WORLD);
      MPI_Send(g.connectCount, count, MPI_INT, rank, 3, MPI_COMM_WORLD);
    }
    int *local_source = (int*)calloc(g.nodes*g.maxConnects, sizeof(int));
    int *local_dest = (int*)calloc(g.nodes*g.maxConnects, sizeof(int));

    for(i = 0; i < count; i++) {
      if(g.source[i]%size == g.dest[i]%size) {
	local_source[i] = g.source[i];
	local_dest[i] = g.dest[i];
      }
    }
    
    return g;
  }
  else{
    /* Receive data from root */
    Network g;
    int i;
    MPI_Status status;
    g.connectCount = (int*)calloc(g.nodes, sizeof(int));
    g.source = (int*)calloc(g.nodes*g.maxConnects, sizeof(int));
    g.dest = (int*)calloc(g.nodes*g.maxConnects, sizeof(int));
    
    MPI_Recv(&g.source, 1, MPI_INT, rank, 1, MPI_COMM_WORLD, &status);
    MPI_Recv(&g.dest, 1, MPI_INT, rank, 2, MPI_COMM_WORLD, &status);
    MPI_Recv(&g.connectCount, 1, MPI_INT, rank, 3, MPI_COMM_WORLD, &status);
  }
}

/**
Book keeping method. Tabulates counts and sets up data transfers 
 */
// Some order assumed: recv buffer contains chunks of recv data in ascending rank order.                   
// That is, buffer looks like [Data from 0, Data from 1, ..., Data from n].                                 
// It is assumed that the order source/destination data is stored is the same on all processors.            
// This assumption saves us having to send companion destination data

Network mpi_tally(Network net){
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int source, dest;
  int sourceRank, destRank;
  int j;

  // Send count arrays so each process knows
  // How many sends/recvs needed with each other node
  net.sendCount = (int*)calloc(size, sizeof(int));
  net.recvCount = (int*)calloc(size, sizeof(int));
  
  // Tabulate size of sends/recvs this process needs to make each update
  for(j=0; j<net.totalConnects;j++){
    source = net.source[j];
    dest = net.dest[j];
    // Remember: (global node number) mod (number of processors) gives 
    // the processor the node belongs to. 
    sourceRank = source%size;
    // Similar to compute which rank destination belongs to
    destRank = dest%size;
    if(sourceRank==rank && destRank!=rank){
      net.sendCount[destRank]++;
    }
    if(destRank==rank && sourceRank!=rank){
      net.recvCount[sourceRank]++;
    }   
  }
  // Count total number of send/recvs this processor needs to make
  net.totalSends = 0;
  net.totalRecvs = 0;
  
  for(j=0; j<size; j++) {
    net.totalSends += net.sendCount[j];
    net.totalRecvs += net.recvCount[j];
  }
  // Book keeping for sends and recevs
  // Sets offset markers to divide send/recv buffers
  // Effectivly reserves enough space for a buffer
  // to/from each processor
  net.sendOffset = (int*) calloc(size, sizeof(int));
  net.recvOffset = (int*) calloc(size, sizeof(int));
  for(j=0;j<size-1;j++) {
    net.sendOffset[j+1] =  net.sendOffset[j]+net.sendCount[j];
    net.recvOffset[j+1] =  net.recvOffset[j]+net.recvCount[j];
  }
  // Book keeping for predetermining node destination of 
  // incoming data 
  net.recvLoc = (int*) calloc(net.totalRecvs, sizeof(int));
  int *count = (int*) calloc(size, sizeof(int));
  // Iterate through all connections once again.
  for(j=0; j<net.totalConnects;j++){
    source = net.source[j];
    dest = net.dest[j];
    destRank = dest%size;
    sourceRank = source%size;
    if(destRank==rank && sourceRank!=rank){
      *(net.recvLoc+(net.recvOffset[sourceRank]+count[sourceRank])) =  dest/size; 
      // store local node index of where data is going. We already know the data is on the right node
      count[sourceRank]++;
    }
  }
  return net;
  
}

/**
Updates the PageRank of each node using the PageRank algorithm
 */
Network  mpi_updatePageRank(Network net){

  // MPI variables
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // Variables
  int j, i, k;
  MPI_Status status;
  int tag = 999;
  int source; // current source
  int dest; // current destination
  int destRank; // Prosessor rank we need to send to
  int sourceRank; // Processor that owns source
  int marker;
  double PR;
  double L;

  // Buffers and counter needed for data transfer
  int *count = (int*) calloc(size, sizeof(int));
  double *recvBuffer = (double*) calloc(net.totalRecvs, sizeof(double));
  double *sendBuffer = (double*) calloc(net.totalSends, sizeof(double));

  /* For all connections */
  for(j=0;j<net.totalConnects;j++){
    source = net.source[j];
    dest = net.dest[j];
    destRank = dest%size;
    sourceRank = source%size;

    // If source belongs to this node, do something with the data
    if(sourceRank==rank){
      PR = net.oldPageRank[source/size];
      L = (double) net.connectCount[source/size];
      if(destRank == rank){
	// If dest also belongs to this node, update its PageRank
	net.newPageRank[dest/size] += net.d*(PR/L);
      }
      else{
	// If dest does not belong to this processor, prepare data to be sent
	*(sendBuffer+(net.sendOffset[destRank]+count[destRank])) = net.d*(PR/L);
	count[destRank]++;	
      }
    }
  }
  marker = 0;
  // Send / Rcv
  // For all processors, loops over whose turn it is to send
  for(j=0; j<size; j++){
    if(j==rank){
      // Send to everyone
      for(i=0; i < size; i++){
	if(i!=j && net.sendCount[i] > 0){ // sendCount to self will be zero, no worries
	  MPI_Send(sendBuffer+net.sendOffset[i], net.sendCount[i], MPI_DOUBLE, 
		   i, tag, MPI_COMM_WORLD); 
	}	
      }
    }
    if(j!=rank && net.recvCount[j] > 0){
      // Recieve from j
      MPI_Recv(recvBuffer+marker, net.recvCount[j], MPI_DOUBLE, j, tag,
	       MPI_COMM_WORLD, &status);
      marker+=net.recvCount[j]; // move marker so we dont overwrite data on next recv
    }   
  }
  int loc; // local node number of where recv data is going
  // Read through recved data and update PageRank of nodes the data belongs to
  for(j=0; j < net.totalRecvs; j++){
    loc = net.recvLoc[j];
    net.newPageRank[loc] += recvBuffer[j];  
  }
  
  /* Update complete. Free buffers */ 
  free(sendBuffer);
  free(recvBuffer);
  free(count);

  return net;
}

/**
Computes the difference norm between the old and new PageRank arrays
 */
double mpi_computeDiff(Network net){
  
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  /* Compute local sum of square differences and perform a reduction
     to sum over all processors, then take square root. 
  */
  int i;
  double local_diff = 0;
  double global_diff;

  for(i = 0; i < net.nodes; i++) {
    local_diff += pow(net.newPageRank[i] - net.oldPageRank[i], 2);
  }
  // if(rank == 0) // MAYBE??
  MPI_Allreduce(&local_diff, &global_diff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return sqrt(global_diff);
}

/**
Updates the PageRank of each node using the PageRank algorithm
 */
Network mpi_computePageRank(Network net){

  double tol = 1e-6;
  double diff = 1;
  double* temp;
  int k;
  double nodes = (double) net.nodes;
  double constant = (1-net.d)/nodes;
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  /* Initialize old PR vlaues */
  for(k=0;k<net.locNodes;k++){
    net.oldPageRank[k] = 1/nodes;
  }

  while(diff > tol){
    /* Set nodes new PR values to constant */
    for(k=0;k<net.locNodes;k++){
      net.newPageRank[k] = constant;
    }

    /* Update PR based on connections */
    net = mpi_updatePageRank(net);
    
    /* Compute norm of difference  */ 
    diff = mpi_computeDiff(net);
    
    /* Have root node print diff */
    if(rank==0){
      printf("diff = %f\n", diff);
    }
    
    /* Switch pointers so updated values in oldPR array             
       Gets us ready to run again */
    temp = net.oldPageRank;
    net.oldPageRank = net.newPageRank;
    net.newPageRank = temp;
  }
  return net;
  
}

/**
Free all (heap) arrays allocated. Make Valgrind happy 
 */
void networkDestructor(Network net){
  free(net.connectCount);
  free(net.source);
  free(net.dest);
  free(net.sendCount);
  free(net.recvCount);
  free(net.sendOffset);
  free(net.recvOffset);
  free(net.oldPageRank);
  free(net.newPageRank);
  free(net.recvLoc);
}

/**
Main Method
 */
int main(int argc, char** argv){
  
  MPI_Init(&argc, &argv);
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // Reads the filename of the data file
  char* filename = argv[1];

  //Start your code
  int j;
  int k = 0;
  
  // Read in data, distribute to all processors
  Network net = mpi_networkReader(filename);
  net.d = 0.85;

  // Tabulation and book keeping for data exchange
  net = mpi_tally(net);

  net = mpi_computePageRank(net);
 
 for(j=0; j< net.nodes;j++){
   // If node belongs to this processor, print its PageRank  
   if(rank == j%size){ 
     printf("PageRank of node %d is: %f\n", j, net.oldPageRank[j/size]);
   }
 }
 // Program complete. Deconstruct network and finalize
 networkDestructor(net);
 MPI_Finalize();

 return 0;
 
}
