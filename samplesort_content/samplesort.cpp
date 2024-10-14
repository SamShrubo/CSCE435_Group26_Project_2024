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

#define MASTER 0               /* taskid of first task */
#define FROM_MASTER 1          /* setting a message type */
#define FROM_WORKER 2          /* setting a message type */

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
void create_input(int input_size, int input_type /* also need pointer to local array or begin and end pointers to entire array*/){
    if(input_type == 0) { // sorted

    } else if (input_type == 1) { // random

    } else if (input_type == 2) { // reverse sorted

    } else if (input_type == 3) { // 1% perturbed

    }
}

// Processor counts: 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024

// perform sample sort in parallel
//// required region annotation
// main
int main(int argc, char* argv[]) {
    CALI_CXX_MARK_FUNCTION;

    // file input should be size of array, input type, and processor count
    int exponent, array_size, inputType;
    if (argc == 3) {
        // take in the size of the array x with 2^x elements
        exponent = atoi(argv[1]);  
        // exponentiate 2 to the power of x to get final array size
        array_size = exp(2, exponent);
        // get the input type, 0 for sorted, 1 for random, 2 for reverse sorted, 3 for 1% perturbed
        inputType = atoi(argv[2]);  // Get the input type
    } else {
        printf("\n Please provide the size of the array and the input type\n");
        return 0;
    }

    // input validation
    if(exponent % 2 != 0 && exponent < 16 && exponent > 28) {
        printf("\n Please provide a valid size of the array as an even number between 16 and 28, inclusive\n");
        return 0;
    } else if(inputType < 0 && inputType > 3) {
        printf("\n Please provide a valid input type between 0 and 3\n");
        printf("0 for sorted, 1 for random, 2 for reverse sorted, 3 for 1%% perturbed\n");
        return 0;
    }

    // variable creation
    double 
        total_time, 
        time_per_process = 0.0;
    int 
        numtasks, 
        taskid;
    int* 
        globalArray = NULL;

    // adiak variables
    string 
        algorithm = "samplesort",
        programming_model = "mpi",
        datatype = "int",
        input_type, // choices: "Sorted", "Random", "Reverse sorted", "1% perturbed"
        scalability, // choices: "weak", "strong"
        implementation_source; // choices: ("online", "ai", "handwritten")
    int 
        group_num = 26,
        size_of_data_type = sizeof(int),
        input_size = array_size,
        num_procs; // number of processors (MPI ranks)

    // offical start of program
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

    //// required region annotation
    // data_init_runtime

    int localArraySize = N / numtasks;
    // TODO: should also handle when numtasks does not cleanly divide N

    int *localArray = (int*) malloc(localArraySize * sizeof(int));

    // generate arrays in each process

    // end data_init_runtime

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

    //// required region annotation:
    // correctness_check (to ensure sorting worked)
        // TODO: call correctness_check
        // if false, sorting was not successful
        // if true, sorting worked
    // end correctness_check

    // required adiak code
    adiak::init(NULL);
    adiak::launchdate();    // launch date of the job
    adiak::libraries();     // Libraries used
    adiak::cmdline();       // Command line used to launch the job
    adiak::clustername();   // Name of the cluster
    adiak::value("algorithm", algorithm); // The name of the algorithm you are using (e.g., "merge", "bitonic")
    adiak::value("programming_model", programming_model); // e.g. "mpi"
    adiak::value("data_type", data_type); // The datatype of input elements (e.g., double, int, float)
    adiak::value("size_of_data_type", size_of_data_type); // sizeof(datatype) of input elements in bytes (e.g., 1, 2, 4)
    adiak::value("input_size", input_size); // The number of elements in input dataset (1000)
    adiak::value("input_type", input_type); // For sorting, this would be choices: ("Sorted", "ReverseSorted", "Random", "1_perc_perturbed")
    adiak::value("num_procs", num_procs); // The number of processors (MPI ranks)
    adiak::value("scalability", scalability); // The scalability of your algorithm. choices: ("strong", "weak")
    adiak::value("group_num", group_number); // The number of your group (integer, e.g., 1, 10)
    adiak::value("implementation_source", implementation_source); // Where you got the source code of your algorithm. choices: ("online", "ai", "handwritten").

    // required performance metrics
    double 
        minimum_time_per_rank,
        maximum_time_per_rank,
        average_time_per_rank,
        total_time,
        variance_time_per_rank;

    /* USE MPI_Reduce to calculate minimum, maximum and the average times for processes.
    MPI_Reduce (&sendbuf,&recvbuf,count,datatype,op,root,comm). https://hpc-tutorials.llnl.gov/mpi/collective_communication_routines/ */
    // use MPI_Reduce to 

    if(taskid == 0) {
        // have the master process conglomerate expected performance metrics

    } else if (taskid == 1) {
        // the worker process might do something? not sure on this yet
    }

    MPI_Finalize();
}
// end main
