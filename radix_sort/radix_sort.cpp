#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <caliper/cali.h>
#include <caliper/cali-manager.h>
#include <adiak.hpp>

const int BASE = 10; // Base for decimal numbers

// Counting sort function (sorting portion only)
void counting_sort(int* arr, int size, int digitPlacement, int* buckets) {
    int output[size]; 
    int count[BASE] = {0};

    // Count occurrences of digits based on current exponent
    for (int i = 0; i < size; i++) {
        count[(arr[i] / digitPlacement) % BASE]++;
    }
    
    for (int i = 0; i < 10; i++){
		buckets[i] = count[i];
    }

    // Transform count array to represent positions in the output array
    for (int i = 1; i < BASE; i++) {
        count[i] += count[i - 1];
    }

    // Build the output array
    for (int i = size - 1; i >= 0; i--) {
        output[count[(arr[i] / digitPlacement) % BASE] - 1] = arr[i];
        count[(arr[i] / digitPlacement) % BASE]--;
    }

    // Copy the output array back to arr[]
    for (int i = 0; i < size; i++) {
        arr[i] = output[i];
    }
}

int main(int argc, char* argv[]) {
    // Initialize MPI environment
    MPI_Init(&argc, &argv);
    
    int taskid, numTasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    MPI_Comm_size(MPI_COMM_WORLD, &numTasks);
    
    if(argc < 2) {
        if(taskid == 0) {
            std::cerr << "Usage: " << argv[0] << " <size_of_array>" << std::endl;
        }
        MPI_Finalize();
        return EXIT_FAILURE;
    }
    
    int sizeOfArray = atoi(argv[1]);
    int sizeOfLocalArray = sizeOfArray / numTasks;
    
    // Initialize profiling tools
    cali::ConfigManager mgr;
    mgr.start();
    
    CALI_MARK_BEGIN("main"); 
    
    srand(time(NULL) + taskid);
    CALI_MARK_BEGIN("data-init-runtime");
    std::vector<int> arr(sizeOfLocalArray);
    for(int i = 0; i < sizeOfLocalArray; i++) {
        arr[i] = rand() % 10000;
    }
    CALI_MARK_END("data-init-runtime");
  
    // Find local maximum
    int local_max = (sizeOfLocalArray > 0) ? arr[0] : 0;
    for(int j = 1; j < sizeOfLocalArray; j++) {
        if(arr[j] > local_max) local_max = arr[j];
    }

    // Find global maximum
    int global_max = 0;
    CALI_MARK_BEGIN("comm");
    CALI_MARK_BEGIN("comm-small");
    MPI_Allreduce(&local_max, &global_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    CALI_MARK_END("comm-small");
    CALI_MARK_END("comm");

    // Calculate the number of digits in the global maximum
    CALI_MARK_BEGIN("comp");
    CALI_MARK_BEGIN("comp-small");
    int numDigits = 0;
    int temp = global_max;
    do {
        numDigits++;
        temp /= 10;
    } while(temp > 0);
    CALI_MARK_END("comp-small");
    CALI_MARK_END("comp");
    
    
    // Radix Sort Loop
    
    int* localBuffer = new int[sizeOfLocalArray];
    int* allBuckets = new int[10 * numTasks];
    for(int d = 0; d < numDigits; d++) {

        int buckets[10] = {0};
        int allDigitsTotal[10] = {0};
        int allDigitsPrefixTotal[10] ={0};
        int allDigitsLeftTotal[10] = {0};
        int LeastDigitTracking[10] = {0}; 
        MPI_Request request;
        MPI_Status status;

        int digitPlacement = static_cast<int>(pow(10, d));
        CALI_MARK_BEGIN("comp");
        CALI_MARK_BEGIN("comp-small");
        counting_sort(arr.data(), sizeOfLocalArray, digitPlacement, buckets);
        CALI_MARK_END("comp-small");
        CALI_MARK_END("comp");
    

        CALI_MARK_BEGIN("comm");
        CALI_MARK_BEGIN("comm-small");
        MPI_Allgather(buckets, 10, MPI_INTEGER, allBuckets, 10, MPI_INTEGER, MPI_COMM_WORLD);
        CALI_MARK_END("comm-small");
        CALI_MARK_END("comm");

        CALI_MARK_BEGIN("comp");
        CALI_MARK_BEGIN("comp-large");
        for (int i = 0; i < 10 * numTasks; i++) {
            int leastDigit = i % 10;
            int rank = i / 10; 
            int num = allBuckets[i];

            if (rank < taskid) {
                allDigitsLeftTotal[leastDigit] += num;
            }
            allDigitsTotal[leastDigit] += num;
            allDigitsPrefixTotal[leastDigit] += num;

        }
        for (int i = 1; i < 10; i++) {
          allDigitsPrefixTotal[i] += allDigitsPrefixTotal[i - 1];
        }
        CALI_MARK_END("comp-large");
        CALI_MARK_END("comp");
        
        
        MPI_Barrier(MPI_COMM_WORLD);
        
        int digitIndex[2];
        int digit, leastDigit, destIdx, destProcess, localDestIdx;
        for (int i = 0; i < sizeOfLocalArray; i++)
        {
            CALI_MARK_BEGIN("comp");
            CALI_MARK_BEGIN("comp-small");
            digit = arr[i];
            leastDigit = (arr[i] / digitPlacement) % 10;

            destIdx = allDigitsPrefixTotal[leastDigit] - allDigitsTotal[leastDigit] + allDigitsLeftTotal[leastDigit] + LeastDigitTracking[leastDigit];
            
            
            LeastDigitTracking[leastDigit]++;
            destProcess = destIdx / sizeOfLocalArray; 

            digitIndex[0] = digit;
            digitIndex[1] = destIdx; 
            CALI_MARK_END("comp-small");
            CALI_MARK_END("comp");
            
          
            CALI_MARK_BEGIN("comm");
            CALI_MARK_BEGIN("comm-large");
            MPI_Send(&digitIndex, 2, MPI_INT, destProcess, 0, MPI_COMM_WORLD);
            MPI_Recv(digitIndex, 2, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            CALI_MARK_END("comm-large");
            CALI_MARK_END("comm");
        
            CALI_MARK_BEGIN("comp");
            CALI_MARK_BEGIN("comp-small");
            localDestIdx = digitIndex[1] % sizeOfLocalArray; 
            localBuffer[localDestIdx] = digitIndex[0]; 
            CALI_MARK_END("comp-small");
            CALI_MARK_END("comp");
        
        }
        CALI_MARK_BEGIN("comp");
        CALI_MARK_BEGIN("comp-small");
        for (int i = 0; i < sizeOfLocalArray; i++)
        {
            arr[i] = localBuffer[i];
        }
        CALI_MARK_END("comp-small");
        CALI_MARK_END("comp");
        
        MPI_Barrier(MPI_COMM_WORLD);
    }
    delete[] localBuffer;
    delete[] allBuckets;
    // Verification 
    CALI_MARK_BEGIN("correctness-check");

    bool local_sorted = true;
    if (!arr.empty()) {
        for(size_t i = 1; i < arr.size(); i++) {
            if(arr[i] < arr[i-1]) {
                return false;
            }
        }
    }

    if(!local_sorted) {
        std::cout << "Process " << taskid << ": Local Array not sorted" << std::endl;
    }
    
    int local_max_value = (sizeOfLocalArray > 0) ? arr[sizeOfLocalArray - 1] : INT_MIN;
    int received_max = INT_MIN;
    
    if (taskid != numTasks - 1) {
        MPI_Send(&local_max_value, 1, MPI_INT, taskid + 1, 1, MPI_COMM_WORLD);
    }
    if (taskid != 0) {
        MPI_Recv(&received_max, 1, MPI_INT, taskid - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (received_max > ((sizeOfLocalArray > 0) ? arr[0] : INT_MIN)) {
            std::cout << "Array not sorted between processor " << taskid - 1 << " and " << taskid << std::endl;
            local_sorted = false;
        }
    }
    
    // Gather the global sorted status
    bool global_sorted = false;
    MPI_Allreduce(&local_sorted, &global_sorted, 1, MPI_CXX_BOOL, MPI_LAND, MPI_COMM_WORLD);
    
    if(global_sorted && taskid == 0) {
        std::cout << "Array Sorted Globally!" << std::endl;
    } else {
        if(!local_sorted) {
            std::cout << "Process " << taskid << ": Local array not sorted or global order incorrect." << std::endl;
        }
    }
    
    // Print the sorted array in order
    int rank = 0;
    while (rank < numTasks) {
        if (rank == taskid) {
            if (rank == 0) {
                printf("Finished Arrays\n");
            }
            
            printf("Processor number %d, array detail:", taskid);
            for(auto currentElement : arr) {
                printf(" %d", currentElement);
            }
            printf("\n");  

            fflush(stdout);
        }
        rank++;
        MPI_Barrier(MPI_COMM_WORLD); 
    }
    CALI_MARK_END("correctness-check");

    CALI_MARK_END("main");
    
    adiak::init(NULL);
    adiak::launchdate();    
    adiak::libraries();    
    adiak::cmdline();       
    adiak::clustername();   
    adiak::value("algorithm", "radix"); // The name of the algorithm
    adiak::value("programming_model", "mpi"); // Programming model used
    adiak::value("data_type", "int"); // Data type of input elements
    adiak::value("size_of_data_type", sizeof(int)); // Size of data type in bytes
    adiak::value("input_size", sizeOfArray); // Number of elements in the input dataset
    adiak::value("input_type", "Random"); // Type of input data
    adiak::value("num_procs", numTasks); // Number of MPI processes
    adiak::value("scalability", "strong"); // Scalability type
    adiak::value("group_num", 26); // Group number
    adiak::value("implementation_source", "online & ai"); // Source of implementation
    
    // Finalize MPI environment
    MPI_Finalize();
    return 0;
}
