/******************************************************************************
* FILE: columnsort.cpp
* DESCRIPTION:
* AUTHOR: Griffin Beaudreau
* LAST REVISED: 10/14/2024
*******************************************************************************/

/*
    Begins with n = r*s values place in an r x s matrix, where r % s = 0 and r >= 2(2-1)^2.
    
    Expected input is the number of values n and the number of processors p.
    r = (n + p - 1) / p
    

    Steps:
        1) Sort the values in each column
        2) transpose the matrix
        3) Sort the values in each column
        4) transpose the matrix
        5) Sort the values in each column
        6) shift the values foward by r/2 positions. Creates an extra column. Fill vacant positions with -infinity
            excess values with +infinity
        7) Sort the values in each column
        8) Unshift the values back by r/2 positions
*/

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <algorithm>
#include <vector>
#include <iostream>
#include <caliper/cali.h>
#include <caliper/cali-manager.h>
#include <adiak.hpp>

#define MASTER 0
const int NEG_INF = INT_MIN;
const int POS_INF = INT_MAX;

// Helpers
void initialize_matrix(std::vector<int>& matrix, int n, int p, int r, int rank) {
    srand(time(NULL) + rank); // Seed that is unique to each process
    for (int i = 0; i < n; ++i) {
        matrix[i] = rand() % 100;
    }
}

void print_matrix(const std::vector<int>& matrix, int r, int p) {
    for (int i = 0; i < r; ++i) {
        for (int j = 0; j < p; ++j) {
            printf("%d ", matrix[i * p + j]);
        }
        printf("\n");
    }
}

void block_transpose(std::vector<int>&local_vector, int r, int p, int rank,  MPI_Comm comm, bool reverse) {
    /*
 
     p = 3
     r = 2

     Transpose Example:
	a0 b0 c0
	a1 b1 c1

     Matrix Index view:
        0 1 2 = rank  
       -------|------
 	0 2 4 | i = 0
	1 3 5 | i = 1

     Result:
	a0 a1 b0
	b1 c0 c1

     Untranspose Example:
     	a0 b0
	a1 b1
	a2 b2

     Matrix Index view:
        0 1 = rank
       -----|------
        0 1 | i = 0       0 3
        2 3 | i = 1  -->  1 4
        4 5 | i = 2       2 5

     Result:
	a0 b1
	b0 a2
	a1 b2
    */

    struct Instruction {
        int value;
	int index;
    };

    std::vector<int> buffer(r);

    // Send instructions
    int matrix_index, target_rank, target_index;
    for (int i = 0; i < r; ++i) {
        if (!reverse) {
            matrix_index = rank * r + i;
	    target_rank = matrix_index % p;
	    //printf("[NORMAL] target_rank: %d = %d / %d", target_rank, matrix_index, p);
	    target_index = matrix_index / p;
	} else {
	    matrix_index = p * i + rank;
	    target_rank = matrix_index / r;
	    //printf("[REVERSE] target_rank: %d = %d / %d", target_rank, matrix_index, p);
	    target_index = matrix_index % r;
	}

	//printf("[SEND] rank:%d, mi:%d, ti:%d, tr:%d,  p:%d, i:%i, r:%d\n", rank, matrix_index, target_index, target_rank, p, i, r);

	if (target_rank == rank) {
	    buffer[target_index] = local_vector[i];
	 //   printf("Rank %d: Self-update, placed %d at buffer[%d]\n", rank, local_vector[i], target_index);
	} else {
	    Instruction instruction = {local_vector[i], target_index};
	   // printf("Rank %d: Sending instruction to rank %d (value=%d, index=%d)\n", rank, target_rank, instruction.value, instruction.index);

	    MPI_Send(&instruction, 2, MPI_INT, target_rank, 0, comm);
	}
    }

    // Retrieve instructions
    MPI_Status status;
    for (int i = 0; i < r; ++i) {
        if (!reverse) {
            target_rank = (rank * r + i) % p;
	} else {
	    target_rank = (p * i + rank) / r;
	}
	//printf("target_rank=%d", target_rank);
	if (target_rank != rank) {
	    Instruction instruction;

	    MPI_Recv(&instruction, 2, MPI_INT, MPI_ANY_SOURCE, 0, comm, &status);
	    buffer[instruction.index] = instruction.value;
	  //  printf("Rank %d: Received instruction from rank %d, placed %d at buffer[%d]\n", rank, target_rank, instruction.value, instruction.index	}
	}
    }

    std::copy(buffer.begin(), buffer.end(), local_vector.begin());
}

void shift_vector(std::vector<int>&local_vector, int r, int p, int rank, int shift, MPI_Comm comm) {
    MPI_Status status;
    std::vector<int> recv_buffer(shift);
    std::vector<int> send_buffer(shift);

    if (p <= 1) return; // No need to shift as all elements are in same column / processor

    // let n = shift =  r/2
    if (rank < p - 1) {
	// Send last n elements to the next processor
	std::copy(local_vector.end() - shift, local_vector.end(), send_buffer.begin());
	MPI_Send(send_buffer.data(), shift, MPI_INT, rank + 1, 0, comm);
    }

    if (rank > 0) {
	MPI_Recv(recv_buffer.data(), shift, MPI_INT, rank - 1, 0, comm, &status);
    } else {
	std::fill(recv_buffer.begin(), recv_buffer.end(), NEG_INF);
    }

    std::vector<int> temp_vector(r + shift);
    std::copy(recv_buffer.begin(), recv_buffer.end(), temp_vector.begin());
    std::copy(local_vector.begin(), local_vector.begin() + (r - shift), temp_vector.begin() + shift);

    if (rank == p - 1) {
    	std::copy(local_vector.end() - shift, local_vector.end(), temp_vector.end() - shift);
    }

    if (rank < p - 1) {
        local_vector.resize(r);
	std::copy(temp_vector.begin(), temp_vector.begin() + r, local_vector.begin());
    } else {
        local_vector = temp_vector;
    }
}

void unshift_vector(std::vector<int>& local_vector, int r, int p, int rank, int shift, MPI_Comm comm) {
    MPI_Status status;
    std::vector<int> send_buffer(shift);
    std::vector<int> recv_buffer(shift);

    if (rank > 0) {
	std::copy(local_vector.begin(), local_vector.begin() + shift, send_buffer.begin());
	MPI_Send(send_buffer.data(), shift, MPI_INT, rank - 1, 0, comm);
    }

    if (rank < p - 1) {
	MPI_Recv(recv_buffer.data(), shift, MPI_INT, rank + 1, 0, comm, &status);
    } else {
	std::copy(local_vector.end() - shift, local_vector.end(), recv_buffer.begin());
    }

    std::vector<int> temp_vector(r);
    std::copy(local_vector.begin() + shift, local_vector.begin() + r, temp_vector.begin()); 
    
    if (rank < p - 1) {
        std::copy(recv_buffer.begin(), recv_buffer.end(), temp_vector.end() - shift);
    } else {
	std::copy(local_vector.end() - shift, local_vector.end(), temp_vector.end() - shift);
    }

    local_vector = temp_vector;
}

void print_vector(const std::vector<int>& vec, int rank) {
    printf("Processor %d's vector:\n", rank);
    for (int val : vec) {
        if (val == NEG_INF) printf("-inf ");
        else if (val == POS_INF) printf("+inf ");
        else printf("%d ", val);
    }
    printf("\n");
}

void distribute_columns(const std::vector<int>& matrix, std::vector<int>& column_data, int r, int p) {
    for (int j = 0; j < p; ++j) {
        for (int i = 0; i < r; ++i) {
            column_data[j * r + i] = matrix[i * p + j];
        }
    }
}

int main(int argc, char *argv[]) {
    CALI_CXX_MARK_FUNCTION;

    // Whole Computation START
    CALI_MARK_BEGIN("main");

    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 3) {
        if (rank == MASTER) {
            printf("Usage: %s <number of elements (n)> <number of processors (p)>\n", argv[0]);
        }
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    // Parse command line arguments
    int n = atoi(argv[1]);
    int p = atoi(argv[2]);
    int r = (n + p - 1) / p;
    int shift = r / 2;

    if (n <= 0 || p <= 0 || p != size) {
        if (rank == MASTER) {
            printf("Error: Invalid input. n and p must be positive integers and p must match the number of MPI processes.\n");
        }
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    if (rank == MASTER) {
        printf("n = %d, p = %d, r = %d, shift = %d\n", n, p, r, shift);
    }

    cali::ConfigManager mgr;
    mgr.start();

    // Create the matrix with filler values (includes padding)
    CALI_MARK_BEGIN("create_matrix");
    std::vector<int> matrix(r * p, POS_INF);
    std::vector<int> local_vector(r);
    std::vector<int> transposed_matrix(r * p, POS_INF);

    if (rank == MASTER) {
        initialize_matrix(matrix, n, p, r, rank);
        printf("Matrix before sorting:\n");
        print_matrix(matrix, r, p);
        printf("\n");

        std::vector<int> column_data(r * p);
        distribute_columns(matrix, column_data, r, p);
        MPI_Scatter(column_data.data(), r, MPI_INT, local_vector.data(), r, MPI_INT, MASTER, MPI_COMM_WORLD);
    } else {
        MPI_Scatter(NULL, r, MPI_INT, local_vector.data(), r, MPI_INT, MASTER, MPI_COMM_WORLD);
    }

    CALI_MARK_END("create_matrix");


    // Step 1 - Sort
    CALI_MARK_BEGIN("step1_sort");
    std::sort(local_vector.begin(), local_vector.end());
    CALI_MARK_END("step1_sort");

    // Step 2 - Inplace Transpose
    CALI_MARK_BEGIN("step2_transpose");
    block_transpose(local_vector, r, p, rank, MPI_COMM_WORLD, false);
    CALI_MARK_END("step2_transpose");

    // Step 3 - Sort
    CALI_MARK_BEGIN("step3_sort");
    std::sort(local_vector.begin(), local_vector.end());
    CALI_MARK_END("step3_sort");

    // Step 4 - Reverse Inplace Tranpose
    CALI_MARK_BEGIN("step4_transpose");
    block_transpose(local_vector, r, p, rank, MPI_COMM_WORLD, true);
    CALI_MARK_END("step4_transpose");

    // Step 5 - Sort
    CALI_MARK_BEGIN("step5_sort");
    std::sort(local_vector.begin(), local_vector.end());
    CALI_MARK_END("step5_sort");

    // Step 6 - Shift
    CALI_MARK_BEGIN("step6_shift");
    shift_vector(local_vector, r, p, rank, shift, MPI_COMM_WORLD);
    CALI_MARK_END("step6_shift");

    // Step 7 - Sort
    CALI_MARK_BEGIN("step7_sort");
    std::sort(local_vector.begin(), local_vector.end());
    CALI_MARK_END("step7_sort");

    // Step 8 - Unshift
    CALI_MARK_BEGIN("step8_unshift");
    unshift_vector(local_vector, r, p, rank, shift, MPI_COMM_WORLD);
    CALI_MARK_END("step8_unshift");

    // Gather column-wise (MPI_Gather defaults to row-wise)
    CALI_MARK_BEGIN("gather_final");
    for(int i = 0; i < r; ++i) {
        MPI_Gather(&local_vector[i], 1, MPI_INT, &matrix[i * p], 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    }
    CALI_MARK_END("gather_final");

    // Whole Computation END
    CALI_MARK_END("main");   

    if (rank == MASTER) {
        printf("Matrix after sorting:\n");
        print_matrix(matrix, r, p);
        printf("\n");
    }
    mgr.stop();
    mgr.flush();
    MPI_Finalize();
    return 0;
}
