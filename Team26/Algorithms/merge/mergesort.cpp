#include <iostream>
#include <mpi.h>
#include <vector>
#include <algorithm>
#include <random>
#include <chrono>
#include <cstdlib>

#include <caliper/cali.h>
#include <caliper/cali-manager.h>
#include <adiak.hpp>


//-----Function Definitions-----//

// Function to merge two sorted arrays
std::vector<int> merge(const std::vector<int>& left, const std::vector<int>& right) {
    std::vector<int> result(left.size() + right.size());
    size_t i = 0, j = 0, k = 0;

    while (i < left.size() && j < right.size()) {
        if (left[i] < right[j]) {
            result[k++] = left[i++];
        } else {
            result[k++] = right[j++];
        }
    }

    // Copy remaining elements
    while (i < left.size()) {
        result[k++] = left[i++];
    }
    while (j < right.size()) {
        result[k++] = right[j++];
    }

    return result;
}

// Function generates the local array portion based on input
void generate_array(std::vector<int>& local_array, int array_size, const std::string& input_type, int rank) {
    // randomization methodology
    std::srand(time(0) + rank);
    

    if (input_type == "Sorted") {
        for (int i = 0; i < array_size; ++i) {
            local_array[i] = i;
        }
    } else if (input_type == "Random") {
        for (int i = 0; i < array_size; ++i) {
            local_array[i] = std::rand() % RAND_MAX;
        }
    } else if (input_type == "Reverse") {
        for (int i = 0; i < array_size; ++i) {
            local_array[i] = array_size - i;
        }
    } else if (input_type == "Perturbed") {
        for (int i = 0; i < array_size; ++i) {
            local_array[i] = i;
        }
        // ensure 1% of elements are modified
        // don't need to check for 0 since we know
        // array size already
        for (int i = 0; i < array_size / 100; ++i) {
            int idx1 = rand() % array_size;
            int idx2 = rand() % array_size;
            std::swap(local_array[idx1], local_array[idx2]);
        }
    }
}

// slow correctness check function because i forgot the faster method
bool correctness_check(const std::vector<int>& sorted_array) {
    for (size_t i = 0; i < sorted_array.size() - 1; ++i) {
        if (sorted_array[i] > sorted_array[i + 1]) {
            return false;  // The array is not sorted correctly
        }
    }
    return true;  // The array is sorted correctly
}

//-----Main-----//

int main(int argc, char** argv) {
    // Create caliper ConfigManager object
    cali::ConfigManager mgr;
    mgr.start();
    
    CALI_CXX_MARK_FUNCTION;
    
    // track main time
    CALI_MARK_BEGIN("main");
    
    CALI_MARK_BEGIN("data_init");
    
    MPI_Init(&argc, &argv);

    int rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    // check cl args
    if (argc != 3) {
        if (rank == 0) {
            std::cerr << "Usage: " << argv[0] << " <power of 2 for array size> <input_type (Sorted, Random, Reverse, Perturbed)>\n";
        }
        MPI_Finalize();
        return 1;
    }

    int power = std::stoi(argv[1]);  // this will be raised as a power
    std::string input_type = argv[2];  // Array generation type

    int total_array_size = 1 << power;  // bit shift to raise to this power
    int local_array_size = total_array_size / num_procs;
    
    
    // have each process make its array
    std::vector<int> local_array(local_array_size);
    generate_array(local_array, local_array_size, input_type, rank);
    CALI_MARK_END("data_init");

    
    CALI_MARK_BEGIN("comp");
    CALI_MARK_BEGIN("comp_small");
    // step 1: sort local array using any simple sort,
    // i am just using the cpp standard sorting algorithm
    std::sort(local_array.begin(), local_array.end());
    CALI_MARK_END("comp_small");
    CALI_MARK_END("comp");
    
    
    // step 2: merge sub-arrays with neighboring processes
    int step = 1;
    while (step < num_procs) {
        if (rank % (2 * step) == 0) {
            if (rank + step < num_procs) {
                CALI_MARK_BEGIN("comm");
                CALI_MARK_BEGIN("comm_large");
                // receive sorted array from neighboring process
                int received_size = local_array_size * step;
                std::vector<int> received_array(received_size);
                
                MPI_Recv(received_array.data(), received_size, MPI_INT, rank + step, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                CALI_MARK_END("comm_large");
                CALI_MARK_END("comm");
                
                CALI_MARK_BEGIN("comp");
                CALI_MARK_BEGIN("comp_large");
                // merge local array with the received array
                local_array = merge(local_array, received_array);
                CALI_MARK_END("comp_large");
                CALI_MARK_END("comp");
            }
        } else {
            CALI_MARK_BEGIN("comm");
            CALI_MARK_BEGIN("comm_large");
            // send local array to neighboring process
            MPI_Send(local_array.data(), local_array.size(), MPI_INT, rank - step, 0, MPI_COMM_WORLD);
            CALI_MARK_END("comm_large");
            CALI_MARK_END("comm");
            break;
        }
        
        // double the step size since our local arrays are all 2x larger
        step *= 2;  
    }
    
    
    CALI_MARK_BEGIN("correctness_check");
    // Step 3: ensure only root process checks the correctness
    if (rank == 0) {
        if (correctness_check(local_array)) {
            std::cout << "The array is sorted correctly." << std::endl;
        } else {
            std::cout << "The array is NOT sorted correctly." << std::endl;
        }
    }
    CALI_MARK_END("correctness_check");
    
    
    // required metadata tracking
    adiak::init(NULL);
    adiak::launchdate();    // launch date of the job
    adiak::libraries();     // Libraries used
    adiak::cmdline();       // Command line used to launch the job
    adiak::clustername();   // Name of the cluster
    adiak::value("algorithm", "Merge Sort"); // The name of the algorithm you are using (e.g., "merge", "bitonic")
    adiak::value("programming_model", "mpi"); // e.g. "mpi"
    adiak::value("data_type", "int"); // The datatype of input elements (e.g., double, int, float)
    adiak::value("size_of_data_type", sizeof(int)); // sizeof(datatype) of input elements in bytes (e.g., 1, 2, 4)
    adiak::value("input_size", total_array_size); // The number of elements in input dataset (1000)
    adiak::value("input_type", input_type); // For sorting, this would be choices: ("Sorted", "ReverseSorted", "Random", "1_perc_perturbed")
    adiak::value("num_procs", num_procs); // The number of processors (MPI ranks)
    adiak::value("scalability", "weak (assumption)"); // The scalability of your algorithm. choices: ("strong", "weak")
    adiak::value("group_num", "26"); // The number of your group (integer, e.g., 1, 10)
    adiak::value("implementation_source", "ai & handwritten"); // Where you got the source code of your algorithm. choices: ("online", "ai", "handwritten").
    
    CALI_MARK_END("main");
    // Flush Caliper output before finalizing MPI
    mgr.stop();
    mgr.flush();

    MPI_Finalize();
    
    
    return 0;
}



