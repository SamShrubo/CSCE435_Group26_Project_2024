#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <limits.h>

#include <caliper/cali.h>
#include <caliper/cali-manager.h>
#include <adiak.hpp>

#define MASTER 0 
void swap(int* i, int* j, int direction) {
    CALI_MARK_BEGIN("comp");
    CALI_MARK_BEGIN("comp_large");
    if (direction == *i > *j) {
        int temp = *i;
        *i = *j;
        *j = temp;
    }
    CALI_MARK_END("comp_large");
    CALI_MARK_END("comp");
}

void merge(int start, int length, int array[], int direction) {
    int half = length / 2;
    if(length > 1){
        for(int i = start; i < start + half; i++){
            swap(&array[i], &array[i] + half, direction);
        }
        merge(start, half, array, direction);
        merge(start + half, half, array, direction);
    }
}

void partial_bitonic(int start, int length, int array[], int direction, int min){
    if (length > min){
        int half = length/2;
        partial_bitonic(start, half, array, 1, min);
        partial_bitonic(start + half, half, array, 0, min);
        merge(start, length, array, direction);
    }
}

void bitonic(int start, int length, int array[], int direction){

    if (length > 1){
        int half = length/2;
        bitonic(start, half, array, 1);
        bitonic(start + half, half, array, 0);
        merge(start, length, array, direction);
    }
}

void check_sorted(int start, int length, int array[], int direction){
	CALI_MARK_BEGIN("correctness_check");
    for(int i = start+1; i < start + length; i++) {
		if(array[i-1] != array[i] && direction != array[i-1] < array[i]) {
			std::cout << array[i-1] << ", " << array[i] << " WRONG" << std::endl;
			return;
		}
	}
	std::cout << "CORRECT" << std::endl;
	CALI_MARK_END("correctness_check");
}

int main(int argc, char *argv[]) {
    CALI_CXX_MARK_FUNCTION;
    MPI_Init(&argc, &argv);
   
    int rank;
    int num_tasks;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);

    CALI_MARK_BEGIN("Whole Computation");

    if(rank == MASTER){ // Master process
        int array_size = std::stoi(argv[1]);
        int *array = new int[array_size];
        int portion =  array_size / (num_tasks);

        CALI_MARK_BEGIN("comm");
        for(int i = 1; i < num_tasks; i++){
            CALI_MARK_BEGIN("comm_small");
            MPI_Send(&portion, 1, MPI_INT, i, 0, MPI_COMM_WORLD); // Tell worker size of their portion
            //MPI_SEND(buf: array.sub(i * portion, (i+1) * portion), dest: i) // Send worker their portion
            CALI_MARK_END("comm_small");
        }

        for(int i = 0; i < portion; i++) {
		array[i] = rand() % portion;
	}
	
	bitonic(0, portion, array, 0);
	check_sorted(0, portion, array, 0);

        for(int i = 1; i < num_tasks; i++){
            CALI_MARK_BEGIN("comm_large");
			std::cout << "receiving" << i  << " at " << (i-1)*portion << " for " << portion << std::endl;
			MPI_Recv(array + (i)*(portion), portion, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);	
            CALI_MARK_END("comm_large");
        }

        CALI_MARK_END("comm");
        partial_bitonic(0, array_size, array, 1, portion);    // Combine worker's array with a partial bitonic  
	check_sorted(0, array_size, array, 1);
	delete[] array;

	return 0;
    } else { // Worker process
        int size;
        int *array;
        MPI_Recv(&size, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // Get size of portion
		
        CALI_MARK_BEGIN("data_init_runtime");
        array = new int[size];  
        for(int i = 0; i < size; i++) {
		array[i] = rand() % size;
	}
        CALI_MARK_END("data_init_runtime");
        
		bitonic(0, size, array, rank % 2); // Sort incre/desc based on worker index
		
		check_sorted(0, size, array, rank % 2);
        MPI_Send(array, size, MPI_INT, MASTER, 0, MPI_COMM_WORLD); // Send sorted sub array to master
		delete[] array;
    	CALI_MARK_END("Whole Computation");
    }

    return 0;
}
