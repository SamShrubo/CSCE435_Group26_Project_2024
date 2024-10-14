// Kaitlyn Griffin
// CSCE 435 - Sample Sort Implementation
// create pseudocode for sample sort using MPI

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>

#include <caliper/cali.h>
#include <caliper/cali-manager.h>
#include <adiak.hpp>

#define MASTER 0  /* taskid of first task */

boolean correctness_check(int *array, int size) {
    for(int i = 0; i < size-1; i++) {
        if(array[i] > array[i + 1]) {
            return false;
        }
        return true;
    }
}

// need a function to create array depending on the type of input and size
// Input sizes: 2^16, 2^18, 2^20, 2^22, 2^24, 2^26, 2^28
// Input types: Sorted, Random, Reverse sorted, 1% perturbed
void create_input(int input_type){}

// Processor count: 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024

// perform sample sort in parallel
int main(int argc, char* argv[]) {
    CALI_CXX_MARK_FUNCTION;

    double total_time, time_per_process = 0.0;
    int numtasks, taskid, array_size, N, inputType;
    int* globalArray = NULL;

    // file input should be size of array, input type, and processor count

    if (argc == 3) {
        // take in the size of the array x with 2^x elements
        array_size = atoi(argv[1]);  
        // exponentiate 2 to the power of x to get final array size
        N = exp(2, array_size);
        // get the input type, 0 for sorted, 1 for random, 2 for reverse sorted, 3 for 1% perturbed
        inputType = atoi(argv[2]);  // Get the input type
    } else {
        printf("\n Please provide the size of the array and the input type\n");
        return 0;
    }

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

    int localArraySize = N / numtasks;
    // TODO: should also handle when numtasks does not cleanly divide N

    int *localArray = (int*) malloc(localArraySize * sizeof(int));

    // generate arrays in each process

    // START OF PARALLEL SECTION
    double total_time_start = MPI_Wtime();
    double process_time_start = MPI_Wtime();
    
    /* 
    from slides:
        Sample sort is used when uniformly distributed assumption is not true
            Distributed to m buckets and sort each with quicksort
            Draw sample of size s
            Sort samples and choose m-1 elements to be splitters
            Split into m buckets and proceed with bucket sort
    */

    // distribute elements of array to m buckets
   
    // sort each local array with quicksort
    
    // draw sample of size s
    
    // sort sample and select pivots

    // globally broadcast pivots to all processes
    
    // partition local arrays into buckets based on pivots
    
    // exchange buckets between processes
    
    // sort received buckets

    // END OF PARALLEL SECTION
    double process_time_end = MPI_Wtime();
    
    // gather sorted subarrays back to master process
    if(taskid == MASTER) {
        // gather all subarrays
        int* sortedArray = (int*)malloc(N * sizeof(int));
        MPI_Gather(localArray, localArraySize, MPI_INT, sortedArray, localArraySize, MPI_INT, MASTER, MPI_COMM_WORLD);
    } else {
        MPI_Gather(localArray, localArraySize, MPI_INT, NULL, 0, MPI_INT, MASTER, MPI_COMM_WORLD);
    }

    // Record total time taken
    double total_time_end = MPI_Wtime();
    total_time = total_time_end - total_time_start;
    time_per_process = process_time_end - process_time_start;

    // use MPI_Send and MPI_Recv to gather all of the process times
    // (similarly to how it was done in previous labs)

    // correctness check (to ensure sorting worked)

    MPI_Finalize();
}
