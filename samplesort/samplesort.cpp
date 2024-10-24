// Kaitlyn Griffin
// CSCE 435 - Sample Sort Implementation

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <algorithm>
#include <string>
#include <ctime>
#include <random>

#include <caliper/cali.h>
#include <caliper/cali-manager.h>
#include <adiak.hpp>

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

void create_random_array(int *localArray, int rank, int local_size){
    std::srand(time(0) + rank);
    for(int i = 0; i < local_size; i++) {
        localArray[i] = std::rand() % RAND_MAX;
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

// function to create array depending on the type of input and size
// Input sizes: 2^16, 2^18, 2^20, 2^22, 2^24, 2^26, 2^28
// Input types: Sorted, Random, Reverse sorted, 1% perturbed
void create_input(int *localArray, int local_size, const std::string& input_type, int rank){
    if(input_type == "Sorted") { // sorted
        create_sorted_array(localArray, rank, local_size);
    } else if (input_type == "Random") { // random
        create_random_array(localArray, rank, local_size);
    } else if (input_type == "ReverseSorted") { // reverse sorted
        create_reverse_sorted(localArray, rank, local_size);
    } else if (input_type == "1_perc_perturbed") { // 1% perturbed
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

// Perform sample sort in parallel
int main(int argc, char* argv[]) {
    CALI_CXX_MARK_FUNCTION;

    // Main start
    CALI_MARK_BEGIN("main");

    // Input should be size of array (exponent) and input type
    int exponent, array_size;
    string inputType;
    if (argc == 3) {
        // take in the size of the array x with 2^x elements
        exponent = atoi(argv[1]);  
        // exponentiate 2 to the power of x to get final array size
        array_size = pow(2, exponent);
        // get the input type
        inputType = argv[2];  // Get the input type
    } else {
        printf("\n Please provide the size of the array and the input type\n");
        return 0;
    }

    bool validInputType = inputType == "Sorted" || inputType == "Random" || inputType == "ReverseSorted" || inputType == "1_perc_perturbed";
    // input validation
     if(exponent % 2 != 0 /*|| exponent < 16*/ || exponent > 28) {
        printf("\n Please provide a valid size of the array as an even number between 16 and 28, inclusive\n");
        return 0;
    } else if(!validInputType) {
        printf("\n Please provide a valid input type.\n");
        printf("Valid options: Sorted, Random, ReverseSorted, 1_perc_perturbed \n");
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
        implementation_source = "AI and Online"; // choices: ("online", "ai", "handwritten")
        // online sources: Class Notes, https://www.geeksforgeeks.org/cpp-program-for-quicksort/ and https://en.wikipedia.org/wiki/Samplesort#:~:text=sequential%%2C%%20sorting%%20algorithm.-,Pseudocode,-%5Bedit%5D)
    int 
        group_number = 26,
        size_of_data_type = sizeof(int),
        input_size = array_size,
        num_procs; // number of processors (MPI ranks)

    // convert to choices wanted in adiak
    if(inputType == "Random") {
        input_type = "Random";
    } else if(inputType == "Sorted") {
        input_type = "Sorted";
    } else if (inputType == "ReverseSorted") {
        input_type = "Reverse sorted";
    } else if (inputType == "1_perc_perturbed") {
        input_type = "1%% perturbed";
    } else {
        printf("invalid input somehow");
        return 0;
    }

    // Processor counts: 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024

    // offical start of program
    CALI_MARK_BEGIN("MPI_Init");
    MPI_Init(&argc, &argv);
    CALI_MARK_END("MPI_Init");

    CALI_MARK_BEGIN("MPI_Comm_rank");
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    CALI_MARK_END("MPI_Comm_rank");

    CALI_MARK_BEGIN("MPI_Comm_size");
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    CALI_MARK_END("MPI_Comm_size");

    num_procs = numtasks;

    cali::ConfigManager mgr;
    mgr.start();

    // data_init_runtime start
    CALI_MARK_BEGIN("data_init_runtime");
    int local_size = array_size / numtasks;
    int* localArray = new int[local_size];

    // distribute elements of array to m buckets
    // generate arrays in each process
    create_input(localArray, local_size, input_type, taskid);

    // end data_init_runtime
    CALI_MARK_END("data_init_runtime");

    // sort each local array with quicksort
    CALI_MARK_BEGIN("comp");
    CALI_MARK_BEGIN("comp_small"); 
    // using small computation because only local arrays are involved
    quicksort(localArray, 0, local_size-1);
    CALI_MARK_END("comp_small");

    // draw sample of size s
    CALI_MARK_BEGIN("comp_small");
    // sample 10% of elements in each array
    double sampling_fraction = 0.1; 
    int sample_size = std::max(1, static_cast<int>(std::round(local_size * sampling_fraction)));
    sample_size = std::min(sample_size, local_size); 

    int* local_sample = new int[sample_size];
    for(int i = 0; i < sample_size; i++) {
        int index = ((i+1) * local_size) / (sample_size + 1);

        if(index >= local_size) {
            index = local_size - 1;
        }
        local_sample[i] = localArray[index];
    }
    CALI_MARK_END("comp_small");
    CALI_MARK_END("comp");

    CALI_MARK_BEGIN("comm");
    CALI_MARK_BEGIN("comm_large");
    
    // gather all samples in master 
    int* gathered_sample = NULL;
    if(taskid == MASTER){
        gathered_sample = new int[sample_size * num_procs];
    }

    // create a dummy buffer for non-master processes
    int dummy = 0;

    CALI_MARK_BEGIN("MPI_Gather");
    // Corrected MPI_Gather call with valid recvbuf for all processes
    MPI_Gather(local_sample, sample_size, MPI_INT, 
               (taskid == MASTER) ? gathered_sample : &dummy, 
               sample_size, MPI_INT, MASTER, MPI_COMM_WORLD);
    CALI_MARK_END("MPI_Gather");

    CALI_MARK_END("comm_large");
    CALI_MARK_END("comm");

    CALI_MARK_BEGIN("comp");

    // master sorts sample and selects m-1 elements to be splitters
    int* splitters = new int[num_procs - 1];
    if(taskid == MASTER){ 
        CALI_MARK_BEGIN("comp_small");
        // gathered_sample is of size sample_size * num_procs
        int gathered_sample_size = sample_size * num_procs;
        quicksort(gathered_sample, 0, gathered_sample_size - 1);

        for(int i = 1; i < num_procs; i++) {
            int sample_index = i * sample_size - 1; 
            if(sample_index >= gathered_sample_size) {
                printf("ERROR: sample index out of range on master process.\n");
                sample_index = gathered_sample_size - 1;
            }
            splitters[i-1] = gathered_sample[sample_index];
        }

        CALI_MARK_END("comp_small");
    }

    CALI_MARK_END("comp");

    // globally broadcast splitters to all processes
    CALI_MARK_BEGIN("comm");
    CALI_MARK_BEGIN("comm_large");
    CALI_MARK_BEGIN("MPI_Bcast");
    MPI_Bcast(splitters, num_procs - 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    CALI_MARK_END("MPI_Bcast");
    CALI_MARK_END("comm_large");
    CALI_MARK_END("comm");

    CALI_MARK_BEGIN("comp");
    // partition local arrays into buckets based on splitters
    // allocate send_counts and send_displs
    int* send_counts = new int[num_procs]();
    for(int i = 0; i < local_size; i++) {
        int bucket = 0;
        while(bucket < (num_procs -1) && localArray[i] > splitters[bucket]){
            bucket++;
        }
        send_counts[bucket]++;
    }

    // calculate send_displs
    int* send_displs = new int[num_procs];
    send_displs[0] = 0;
    for(int i = 1; i < num_procs; i++) {
        send_displs[i] = send_displs[i-1] + send_counts[i-1];
    }

    // allocate send_buffer
    int total_send = 0;
    for(int i = 0; i < num_procs; i++) {
        total_send += send_counts[i];
    }
    int* send_buffer = new int[total_send];

    // initialize current positions for each bucket
    int* current_positions = new int[num_procs];
    for(int i = 0; i < num_procs; i++) {
        current_positions[i] = send_displs[i];
    }

    // populate send_buffer
    for(int i = 0; i < local_size; i++) {
        int bucket = 0;
        while(bucket < (num_procs - 1) && localArray[i] > splitters[bucket]){
            bucket++;
        }
        send_buffer[current_positions[bucket]++] = localArray[i];
    }

    CALI_MARK_END("comp");

    CALI_MARK_BEGIN("comm");

    // exchange buckets between processes
    // all-to-all communication to get recv_counts
    int* recv_counts_array = new int[num_procs];
    CALI_MARK_BEGIN("comm_small");
    CALI_MARK_BEGIN("MPI_Alltoall");
    MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts_array, 1, MPI_INT, MPI_COMM_WORLD);
    CALI_MARK_END("MPI_Alltoall");
    CALI_MARK_END("comm_small");
    CALI_MARK_END("comm");

    CALI_MARK_BEGIN("comp");
    CALI_MARK_BEGIN("comp_small");

    // calculate recv_displs and total_recv
    int* recv_displs = new int[num_procs];
    recv_displs[0] = 0;
    for(int i = 1; i < num_procs; i++) {
        recv_displs[i] = recv_displs[i-1] + recv_counts_array[i-1];
    }
    int total_recv = 0;
    for(int i = 0; i < num_procs; i++) {
        total_recv += recv_counts_array[i];
    }

    CALI_MARK_END("comp_small");
    CALI_MARK_END("comp");

    CALI_MARK_BEGIN("comm");
    // allocate recv_buffer
    int* recv_buffer = new int[total_recv];

    // perform all-to-all variable communication  
    CALI_MARK_BEGIN("comm_small");
    CALI_MARK_BEGIN("MPI_Alltoallv");
    MPI_Alltoallv(send_buffer, send_counts, send_displs, MPI_INT,
                 recv_buffer, recv_counts_array, recv_displs, MPI_INT,
                 MPI_COMM_WORLD);
    CALI_MARK_END("MPI_Alltoallv");
    CALI_MARK_END("comm_small");
    CALI_MARK_END("comm");

    // sort received buckets
    CALI_MARK_BEGIN("comp");
    CALI_MARK_BEGIN("comp_large");
    quicksort(recv_buffer, 0, total_recv - 1);
    CALI_MARK_END("comp_large");
    CALI_MARK_END("comp");

    // gather sorted subarrays back to master process
    int* recv_counts_gathered = NULL;
    int* recv_displs_gathered = NULL;
    if(taskid == MASTER){
        recv_counts_gathered = new int[num_procs];
    }

    CALI_MARK_BEGIN("comm");
    CALI_MARK_BEGIN("comm_large");
    CALI_MARK_BEGIN("MPI_Gather");
    MPI_Gather(&total_recv, 1, MPI_INT, recv_counts_gathered, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    CALI_MARK_END("MPI_Gather");
    CALI_MARK_END("comm_large");
    CALI_MARK_END("comm");

    if(taskid == MASTER){
        recv_displs_gathered = new int[num_procs];
        recv_displs_gathered[0] = 0;
        for(int i = 1; i < num_procs; i++) {
            recv_displs_gathered[i] = recv_displs_gathered[i-1] + recv_counts_gathered[i-1];
        }
    }

    int* sortedArray = NULL;
    if(taskid == MASTER) {
        sortedArray = new int[array_size];
    }

    CALI_MARK_BEGIN("comm");
    CALI_MARK_BEGIN("comm_small");
    CALI_MARK_BEGIN("MPI_Gatherv");
    MPI_Gatherv(recv_buffer, total_recv, MPI_INT,
                sortedArray, recv_counts_gathered, recv_displs_gathered, MPI_INT,
                MASTER, MPI_COMM_WORLD);
    CALI_MARK_END("MPI_Gatherv");
    CALI_MARK_END("comm_small");
    CALI_MARK_END("comm");

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
    }

    // Required adiak code
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

    mgr.stop();
    mgr.flush();

    // clean up dynamically allocated memory
    delete[] localArray;
    delete[] local_sample;
    if(gathered_sample){ 
        delete[] gathered_sample;
    }
    delete[] splitters;
    delete[] send_counts;
    delete[] send_displs;
    delete[] send_buffer;
    delete[] current_positions;
    delete[] recv_counts_array;
    delete[] recv_displs;
    delete[] recv_buffer;
    if(taskid == MASTER){
        delete[] recv_counts_gathered;
        delete[] recv_displs_gathered;
        delete[] sortedArray;
    }

    CALI_MARK_BEGIN("MPI_Finalize");
    MPI_Finalize();
    CALI_MARK_END("MPI_Finalize");

    CALI_MARK_END("main");
} // end main

