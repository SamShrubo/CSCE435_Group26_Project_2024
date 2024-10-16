// Kaitlyn Griffin
// CSCE 435 - Sample Sort Implementation
// create pseudocode for sample sort using MPI

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <algorithm>
#include <string>

#include <caliper/cali.h>
#include <caliper/cali-manager.h>
#include <adiak.hpp>

using std::vector;
using std::swap;
using std::string;

#define MASTER 0               /* taskid of first task */
#define FROM_MASTER 1          /* setting a message type */
#define FROM_WORKER 2          /* setting a message type */

bool correctness_check(int *array, int size) {
    for(int i = 0; i < size-1; i++) {
        if(array[i] > array[i + 1]) {
            return false;
        }
    }
    return true;
}

void create_sorted_array(int *localArray, int rank, int local_size){
    int start = rank * local_size;
    for(int i = 0; i < local_size; i++) {
        localArray[i] = start + i;
    }
}

void create_random_array(int *localArray, int local_size){
    for(int i = 0; i < local_size; i++) {
        localArray[i] = std::rand();
    }
}

void create_reverse_sorted(int *localArray, int rank, int local_size){
    int start = rank * local_size;
    for(int i = 0; i < local_size; i++) {
        localArray[i] = start + (local_size - i - 1);
    }
}

void create_one_percent_perturbed(int *localArray, int rank, int local_size){
    create_sorted_array(localArray, rank, local_size);
    // perturb the array
    int num_to_perturb = (int)(std::round(local_size * 0.01));
    if(num_to_perturb == 0) { // ensure array will always be perturbed in some way
        num_to_perturb = 1;
    }
    for(int i = 0; i < num_to_perturb; i++) {
        std::swap(localArray[rand() % local_size], localArray[rand() % local_size]);
    }
}

// need a function to create array depending on the type of input and size
// Input sizes: 2^16, 2^18, 2^20, 2^22, 2^24, 2^26, 2^28
// Input types: Sorted, Random, Reverse sorted, 1% perturbed
void create_input(int *localArray, int local_size, int input_type, int rank){
    if(input_type == 0) { // sorted
        create_sorted_array(localArray, rank, local_size);
    } else if (input_type == 1) { // random
        create_random_array(localArray, local_size);
    } else if (input_type == 2) { // reverse sorted
        create_reverse_sorted(localArray, rank, local_size);
    } else if (input_type == 3) { // 1% perturbed
        create_one_percent_perturbed(localArray, rank, local_size);
    }
}

// referenced from https://www.geeksforgeeks.org/cpp-program-for-quicksort/
// with slight adjustments to use an int array instead of vector<int>
int partition(int* arr, int low, int high){
    int pivot = arr[high];
    int i = (low - 1);
    for (int j = low; j <= high - 1; j++) {
        if (arr[j] <= pivot) {
            i++;
            std::swap(arr[i], arr[j]);
        }
    }
    std::swap(arr[i + 1], arr[high]);
    return (i + 1);
}

// referenced from https://www.geeksforgeeks.org/cpp-program-for-quicksort/
// with slight adjustments to use an int array instead of vector<int>
void quicksort(int *localArray, int low, int high) {
    if(low < high) {
        int pi = partition(localArray, low, high);

        quicksort(localArray, low, pi-1);
        quicksort(localArray, pi+1, high);
    }
}

void printArray(int* arr, int size, int rank) {
    printf("Processor %d: \n", rank);
    printf("size: %d\n", size);
    for(int i = 0; i < size; i++) {
        printf("%d ", arr[i]);
    }
    printf("\n");
}

// perform sample sort in parallel
int main(int argc, char* argv[]) {
    CALI_CXX_MARK_FUNCTION;

    // main start
    CALI_MARK_BEGIN("main");

    // file input should be size of array, input type, and processor count
    int exponent, array_size, inputType;
    if (argc == 3) {
        // take in the size of the array x with 2^x elements
        exponent = atoi(argv[1]);  
        // exponentiate 2 to the power of x to get final array size
        array_size = pow(2, exponent);
        // get the input type, 0 for sorted, 1 for random, 2 for reverse sorted, 3 for 1% perturbed
        inputType = atoi(argv[2]);  // Get the input type
    } else {
        printf("\n Please provide the size of the array and the input type\n");
        return 0;
    }

    // input validation
    if(exponent % 2 != 0 /*|| exponent < 16*/ || exponent > 28) {
        printf("\n Please provide a valid size of the array as an even number between 16 and 28, inclusive\n");
        return 0;
    } else if(inputType < 0 || inputType > 3) {
        printf("\n Please provide a valid input type between 0 and 3\n");
        printf("0 for sorted, 1 for random, 2 for reverse sorted, 3 for 1%% perturbed\n");
        return 0;
    }

    // variable creation
    int 
        numtasks, 
        taskid;

    // adiak variables
    string 
        algorithm = "samplesort",
        programming_model = "mpi",
        data_type = "int",
        input_type, // choices: "Sorted", "Random", "Reverse sorted", "1% perturbed"
        scalability = "strong", // choices: "weak", "strong"
        implementation_source = "temp //TOCHANGE"; // choices: ("online", "ai", "handwritten")
    int 
        group_number = 26,
        size_of_data_type = sizeof(int),
        input_size = array_size,
        num_procs; // number of processors (MPI ranks)

    switch(inputType) {
        case 0:
            input_type = "Sorted";
            break;
        case 1:
            input_type = "Random";
            break;
        case 2:
            input_type = "Reverse sorted";
            break;
        case 3:
            input_type = "1%% perturbed";
            break;
    }

    // Processor counts: 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024

    // offical start of program
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

    cali::ConfigManager mgr;
    mgr.start();

    // START OF PARALLEL SECTION

    // data_init_runtime start
    CALI_MARK_BEGIN("data_init_runtime");
    int local_size = array_size / numtasks;
    int* localArray = new int[local_size];

    // distribute elements of array to m buckets
    // generate arrays in each process
    create_input(localArray, local_size, inputType, taskid);

    // end data_init_runtime
    CALI_MARK_END("data_init_runtime");

    // sort each local array with quicksort
    CALI_MARK_BEGIN("comp");
    CALI_MARK_BEGIN("comp_small"); 
    // using small computation because only the local arrays are involved
    quicksort(localArray, 0, local_size-1);
    CALI_MARK_END("comp_small");
    CALI_MARK_END("comp");

    printf("After first quicksort on processor %d\n", taskid);
    
    // draw sample of size s
    int sample_size = num_procs - 1;
    int* local_sample = new int[sample_size];
    for(int i = 0; i < sample_size; i++) {
        int index = ((i+1) * local_size) / (sample_size + 1);
        // choosing evenly spaced samples from throughout the local array
        local_sample[i] = localArray[index];
    }

    printf("After sample selection on processor %d\n", taskid);
    
    
    printf("at barrier");
    MPI_Barrier(MPI_COMM_WORLD);
    // gather all samples in master 
    int* gathered_sample = NULL;
    if(taskid == MASTER){
        gathered_sample = new int[sample_size * num_procs];
    }
    CALI_MARK_BEGIN("comm");
    CALI_MARK_BEGIN("comm_small");
    MPI_Gather(local_sample, sample_size, MPI_INT, gathered_sample, sample_size, MPI_INT, MASTER, MPI_COMM_WORLD);
    CALI_MARK_BEGIN("comm");
    CALI_MARK_BEGIN("comm_small");


    printf("before master sorts sample \n");

    // master sorts sample and selects m-1 elements to be splitters
    int* splitters = new int[sample_size];
    if(taskid == MASTER){
        quicksort(gathered_sample, 0, sample_size * num_procs - 1);
        for(int i = 1; i < num_procs; i++) {
            splitters[i-1] = gathered_sample[i * sample_size];
        }

        printf("After splitter selection in master\n");
        printArray(gathered_sample, sample_size * num_procs, taskid);
    }

    printf("After master sorts sample on processor %d\n", taskid);

    // globally broadcast splitters to all processes
    CALI_MARK_BEGIN("comm");
    CALI_MARK_BEGIN("comm_small");
    // using small computation because only the local arrays are involved
    MPI_Bcast(splitters, sample_size, MPI_INT, MASTER, MPI_COMM_WORLD);
    CALI_MARK_END("comm_small");
    CALI_MARK_END("comm");
    
    /* 
    from slides:
        Sample sort is used when uniformly distributed assumption is not true
            Distributed to m buckets and sort each with quicksort
            Draw sample of size s
            Sort samples and choose m-1 elements to be splitters
            Split into m buckets and proceed with bucket sort
    */

   printf("partition local arrays into buckets based on splitters on processor %d\n", taskid);

    // TODO: partition local arrays into buckets based on splitters
    // Allocate send_counts and send_displs
    int* send_counts = new int[num_procs]();
    for(int i = 0; i < local_size; i++) {
        int bucket = 0;
        while(bucket < sample_size && localArray[i] > splitters[bucket]){
            bucket++;
        }
        send_counts[bucket]++;
    }

    // Calculate send_displs
    int* send_displs = new int[num_procs];
    send_displs[0] = 0;
    for(int i = 1; i < num_procs; i++) {
        send_displs[i] = send_displs[i-1] + send_counts[i-1];
    }

    // Allocate send_buffer
    int total_send = 0;
    for(int i = 0; i < num_procs; i++) {
        total_send += send_counts[i];
    }
    int* send_buffer = new int[total_send];

    // Initialize current positions for each bucket
    int* current_positions = new int[num_procs];
    for(int i = 0; i < num_procs; i++) {
        current_positions[i] = send_displs[i];
    }

    // Populate send_buffer
    for(int i = 0; i < local_size; i++) {
        int bucket = 0;
        while(bucket < sample_size && localArray[i] > splitters[bucket]){
            bucket++;
        }
        send_buffer[current_positions[bucket]++] = localArray[i];
    }

    printf("Right before exchange buckets between processes in process %d\n", taskid);

    // TODO: exchange buckets between processes
    // All-to-all communication to get recv_counts
    int* recv_counts_array = new int[num_procs];
    CALI_MARK_BEGIN("comm");
    CALI_MARK_BEGIN("comm_small");
    MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts_array, 1, MPI_INT, MPI_COMM_WORLD);
    CALI_MARK_END("comm_small");
    CALI_MARK_END("comm");

    // Calculate recv_displs and total_recv
    int* recv_displs = new int[num_procs];
    recv_displs[0] = 0;
    for(int i = 1; i < num_procs; i++) {
        recv_displs[i] = recv_displs[i-1] + recv_counts_array[i-1];
    }
    int total_recv = 0;
    for(int i = 0; i < num_procs; i++) {
        total_recv += recv_counts_array[i];
    }

    // Allocate recv_buffer
    int* recv_buffer = new int[total_recv];

    // Perform all-to-all variable communication
    CALI_MARK_BEGIN("comm");
    CALI_MARK_BEGIN("comm_large");
    MPI_Alltoallv(send_buffer, send_counts, send_displs, MPI_INT,
                 recv_buffer, recv_counts_array, recv_displs, MPI_INT,
                 MPI_COMM_WORLD);
    CALI_MARK_END("comm_large");
    CALI_MARK_END("comm");

    printf("Sort received buckets in process %d\n", taskid);
    
    // TODO: sort received buckets
    CALI_MARK_BEGIN("comp");
    CALI_MARK_BEGIN("comp_large");
    quicksort(recv_buffer, 0, total_recv - 1);
    CALI_MARK_END("comp_large");
    CALI_MARK_END("comp");

    printf("Gather sorted subarrays back to master process %d\n", taskid);
    
    // gather sorted subarrays back to master process
    int* sortedArray = NULL;
    if(taskid == MASTER) {
        sortedArray = new int[array_size];
    }

    CALI_MARK_BEGIN("comm");
    CALI_MARK_BEGIN("comm_large");
    // gather all subarrays
    MPI_Gather(localArray, local_size, MPI_INT, sortedArray, local_size, MPI_INT, MASTER, MPI_COMM_WORLD);
    CALI_MARK_END("comm");
    CALI_MARK_END("comm_large");

    printf("Subarrays gathered in process %d\n", taskid);
    
    bool correct_check = true;
    if(taskid == MASTER) {
        // correctness_check start
        CALI_MARK_BEGIN("correctness_check");
        // call correctness_check
        correct_check = correctness_check(sortedArray, array_size);  
        if(correct_check) { // if true, sorting worked
            printf("\nSamplesort was successful.\n");
        } else { // if false, sorting was not successful
            printf("\nSamplesort was NOT successful.\n");
        }
        // end correctness_check
        CALI_MARK_END("correctness_check");
        printArray(sortedArray, array_size, taskid);
    }
    
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
    /*double 
        minimum_time_per_rank,
        maximum_time_per_rank,
        average_time_per_rank,
        total_time,
        variance_time_per_rank;
    */

   printf("End of main in process %d\n", taskid);

    CALI_MARK_END("main");

    mgr.stop();
    mgr.flush();

    MPI_Finalize();
}
// end main
