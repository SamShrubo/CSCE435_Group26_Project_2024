# CSCE 435 Group project

## 0. Group number: 26

## 1. Group members:
1. Caroline Jia
   - Username: cjia2003
   - Algorithm: Radix Sort
2. Griffin Beaudreau
   - Username: CyberGriffin
   - Algorithm: Column Sort
3. Kaitlyn Griffin
   - Username: kaitlyngrif
   - Algorithm: Sample Sort
4. Samuel Bush
   - Username: SamShrubo
   - Algorithm: Merge Sort
5. Zhongyou Wu
   - Username: ZhongyouWuTAMU
   - Algorithm: Bitonic Sort

### 1a. Team Communication:
We will be using Discord for our team communications.

## 2. Project topic: Parallel Sorting Algorithms

### 2a. Brief project description (what algorithms will you be comparing and on what architectures)

- Bitonic Sort: A sorting algorithm that scales well with parallel computing architecture. The algorithm requires the array to be size of 2^n where n is a positive integer value. It process the array recursively as subarrays and sort them into bitonic arrays. It will be implemented on the Grace Cluster using MPI.
- Sample Sort: Implemented using MPI on the Grace cluster. Appropriately sized arrays will be generated on each process and sorted using quicksort. Samples are taken from each process's local array, and a sample is then taken from the set of samples. Using the final sample, array elements are sorted into buckets, rearranged, and sorted locally. All sub-arrays are conglomerated back together, resulting in a sorted global array.   
- Merge Sort: Implement using MPI on the Grace cluster, Multiple sub-arrays that will make up the final sorted array will be generated within each process within MPI, from there they are sorted locally using std::sort, then the process of merging each process with its neighboring processes in a loop continues until the final array is produced and all sub-arrays are merged.
- Radix Sort: Implemented using MPI on Grace cluster. The inital array will be split into multiple smaller arrays across the nodes and processors, and sorted locally using a counting sort on the current digit. Afterwards a prefix sum will be caluclated to help put the integers back into a paritally sorted array until all the digit places have een sorted resulting in a sorted global array. 
- Column Sort: A parallel sorting algorithm that is well suited for sorting data arranged in a 2D grid. The matrix is sorted column-wise, transposed, and sorted again row-wise. This process is repeated until the matrix is sorted. This algorithm will be implemented using MPI on the Grace cluster.

### 2b. Pseudocode for each parallel algorithm
- For MPI programs, include MPI calls you will use to coordinate between processes

**---Merge Sort Pseudocode---**

```c++
# Full sorting algorithm
def MergeSort(Array, arraySize) {
   MPI_Init(arguments to set up mpi)
   MPI_Comm_rank(MPI_COMM_WORLD, rank)        # set rank to this process rank
   MPI_Comm_size(MPI_COMM_WORLD, numProcs)    # get number of processes

   # divide array into local sub-arrays for each process to sort individually
   if arraySize % numProcs != 0:
      # round up for array size and any unfilled space in the sub-array can be accounted for as null   
      localArraySize = ceiling(arraySize / numProcs) 
   else:
      localArraySize = arraySize / numProcs
      
   # make local array
   arrayOffset = rank * localArraySize
   localArray = array[localArraySize]    # create buffer for receiving scattered data

   # scatter to all processes
   MPI_Scatter(Array, localArraySize, MPI_INT, localArray, localArraySize, MPI_INT, root=0, MPI_COMM_WORLD)

   # each process sorts locally using quicksort
   localQuickSort(localArray, localArraySize)

   # begin merging with neighbor processes
   step = 1
   while step < num_procs:
      # combine even and odd processes in the even process
      if (rank % (2 * step) == 0):
         if (rank + step < num_procs):
            # get sorted array from neighbor process
            receivedSize = localArraySize * step
            receivedArray = new array[receivedSize]

            MPI_Recv(receivedArray, receivedSize, MPI_INT, rank + step, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE)
   
            # Merge local array and received array using helper function defined below
            localArray = Merge(localArray, localArraySize, receivedArray, receivedSize)
            # Update local size after merge
            localArraySize = localArraySize + receivedSize  
   
      elif (rank % step == 0):
         # Send local array to neighboring process
         MPI_Send(localArray, localArraySize, MPI_INT, rank - step, 0, MPI_COMM_WORLD)

      # double step size each iteration
      step = step * 2 

   # only main process gets final sorted array
   if rank == 0 
      sortedArray = array[arraySize]
   else
      sortedArray = null

   # from all scattered processes gather into final sorted array
   MPI_Gather(localArray, localArraySize, MPI_INT, sortedArray, localArraySize, MPI_INT, root=0, MPI_COMM_WORLD)

   # only main can display / output the sorted array
   if rank == 0
      Display(sortedArray)

   # call in every process
   MPI_Finalize()
}

# helper function merges 2 already sorted arrays into 1
def Merge(Array1, array1Size, Array2, array2Size) {
   mergedArray = array[array1Size + array2Size]
   i = 0
   j = 0
   k = 0

   # compare each element until completed 1 array
   while (i < array1Size) and (j < array2Size)
      if Array1[i] < Array2[j]
         mergedArray[k] = Array1[i]
         i++
      else
         mergedArray[k] = Array2[j]
         j++
      k++

   # if unread elements in Array1 copy them to the end of mergedArray
   while i < array1Size
      mergedArray[k] = Array1[i]
      i++
      k++

   # if unread elements in Array2 copy them to the end of mergedArray
   while j < array2Size
      mergedArray[k] = Array2[j]
      j++
      k++

   return mergedArray
}
```

**---Sample Sort Pseudocode---**
```
// perform sample sort in parallel
int main(int argc, char* argv[]) {
    CALI_CXX_MARK_FUNCTION;

    double total_time, time_per_process = 0.0;
    int numtasks, taskid, N;
    int* globalArray = NULL;

    if (argc == 2) {
        N = atoi(argv[1]);  // Get the size of the array
    } else {
        printf("\n Please provide the size of the array\n");
        return 0;
    }

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

    int localArraySize = N / numtasks;
    // should also handle when numtasks does not cleanly divide N

    int *localArray = (int*) malloc(localArraySize * sizeof(int));

    if(taskid = MASTER) { // only generate array in master process
        globalArray = generate_input_array(N);
    }

    // START OF PARALLEL SECTION
    double total_time_start = MPI_Wtime();
    double process_time_start = MPI_Wtime();
    
    // distribute elements of array to m buckets
    MPI_Scatter(globalArray, localArraySize, MPI_INT, localArray, localArraySize, MPI_INT, MASTER, MPI_COMM_WORLD);

    // sort each local array with quicksort
    quicksort(localArray, 0, localArraySize - 1);

    // draw sample of size s
    int sampleSize = numtasks - 1;
    int* localSamples = select_samples(localArray, sampleSize);

    // gather samples 
    int sampleSize = numtasks - 1;
    int* localSamples = select_samples(localArray, sampleSize);

    // sort sample and select pivots

    // globally broadcast pivots to all processes
    MPI_Bcast(pivots, numtasks - 1, MPI_INT, MASTER, MPI_COMM_WORLD);

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

    MPI_Finalize();
}
```

**---Radix Sort Pseudocode---**
```
MPI_Init()
MPI_Comm_rank(comm, rank);
MPI_Comm_size(comm, size);

// Scatter the array across processes
if (rank == 0) {
    int *arr = generate_input_array(N);
    // Scatter the array to all processes
    MPI_Scatter(variables);
} else {
    // Other processes
    MPI_Scatter(variables);
}

// Perform radix sort on the local portion of the array

while i < max_digits {
    // Perform counting sort at current digit
    int local_count = counting_sort_by_digit(local_arr, local_size, digit_pos, base);

    // Gather global counts from other processes
    MPI_Allgather(local_count, base, MPI_INT, global_count, base, MPI_INT, comm);
    
    // Compute prefix sums on global counts to determine offsets
    int prefix_sum = compute_prefix_sums(global_count, base);
    
    // Redistribute elements based on the computed prefix sums
    int sorted_local_arr = redistribute_elements(local_arr, local_size, digit_pos, prefix_sum, base);

    // Replace local array with the newly sorted portion
    local_arr = sorted_local_arr;
}

// Gather the locally sorted arrays back into the root process
if (rank == 0) {
    MPI_Gather(local_arr, local_size, MPI_INT, sorted_arr, local_size, MPI_INT, 0, MPI_COMM_WORLD);
} else {
    MPI_Gather(local_arr, local_size, MPI_INT, NULL, local_size, MPI_INT, 0, MPI_COMM_WORLD);
}
MPI_Finalize();

```

**---Column Sort Pseudocode---**
```
Steps:
 1: Arrange data in a matrix with r rows and c columns, where r is the number of processors and c is the numner of items per process.
 2: Sort each column using a sequential sorting algorithm (may change algorithm depending on input size).
 3: Transpose the matrix.
 4: Sort each row independently
 5: Sort each column again.
 6: Tranpose the matrix.
 7: Final column sort.
```
```c++
/*
* Include MPI header
* Include Caliper header
* Include any additional headers
* Define Constants (MASTER)
*/

int main(int argc, char *argv[]) {
   CALI_CXX_MARK_FUNCTION;

   if args invalid return 0;

   // Set up MPI environment (needs to include arguments)
   MPI_INIT();
   MPI_Comm_rank()
   MPI_Comm_size()

   // Initialize variables
   int N = atoi(argv[1])
   int P = num processors;
   int rows, N_padded, padding_size;

   // To handle cases where the total number of data elements isn't divisible by the number of processors
   // Example: N = 10 and P = 4:
   // rows = (10 + 4 - 1) / 4 = 3
   // N_padded = 3 * 4 = 12 (what it will be when padded)
   // padding_size = 12 - 10 = 2 (num of padding elements)
   rows = (N + P - 1) / P;
   N_padded = rows * P;
   padding_size = N_padded - N;

   // Local data for each processor
   int *local_array = (int*)malloc(rows * sizeof(int));

   if (taskid == MASTER) {
      int *global_data = (int*)malloc(N_padded * sizeof(int));
      
      /* Add actual data to global_data */

      /* Add padded data */
      for (int i = N; i < N_padded; ++i) {
         global_data[i] = INT_MAX;
      }

      for (int p = 0; p < P; ++p) {
         if (p == MASTER) {
            // copy local data to global data
         } else {
            // MPI_Send
         }
      }
      free(global_data)
   } else {
      // MPI_Recv
   }

   /* local column sort, marked with caliper */

   /* transpose, marked with caliper */

   /* lcaol row sort, marked with caliper */

   /* transpose, marked with caliper */

   /* local column sort, marked with caliper */

   if (taskid == MASTER) {
      int *sorted_data = (int*)malloc(N * sizeof(int));
      // copy local_data into sorted_data

      foreach p = 1; p < column, MPI_Recv(sorted_data[p*rows],rows,MPI_INIT,p,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

      // Verify that array is sorted

      free(sorted_data);
   } else {
      MPI_SNED(local_data, rows, MPI_INT, MASTER, 0, MPI_COMM_WORLD);
   }

   free(local_data);

   mgr.flush();
   MPI_Finalize();

   return 0;
}
```

**---Bitonic Sort Pseudocode---**
```c++
swap(int* i, int* j, int direction):
   if direction == *i > *j:
      temp = *i
      i = *j
      *j = temp

merge(int start, int length, int array[], int direction):
   half = length / 2
   if(length > 1):
      for(int i = start, half):
         swap(array[i], array[i] + half, direction)
      merge(start, length, array, direction)

master_bitonic(int start, int length, int array[], int direction, int min):
   if length > min:
      int half = length/2
      master_bitonic(start, half, array, 1, min)
      master_bitonic(start + half, half, array, 0, min)
      merge(start, length, array, direction)

bitonic(int start, int length, int array[], int direction):
   if length > 1:
      int half = length/2
      bitonic(start, half, array, 1)
      bitonic(start + half, half, array, 0)
      merge(start, length, array, direction)

main():
   MPI_init()
   rank = MPI_Comm_rank()
   size = MPI_Comm_size()
   

   if(rank == MASTER): // Master process
      array = input()
      portion =  array.size / size

      for(int i = 0, size):
         MPI_SEND(buf: portion, dest: i) // Tell worker size of their portion
         MPI_SEND(buf: array.sub(i * portion, (i+1) * portion), dest: i) // Send worker their portion
      
      for(int i = 0, size):
         MPI_RECV(buf: array + i*(portion), count: portion, source: i)

      master_bitonic(0, size, array, 1, portion)    // Combine worker's array with a partial bitonic  

   else : // Worker process
      size
      array
      MPI_RECV(buf: size, count: 1, source: MASTER) // Get size of portion
      MPI_RECV(buf: array, count: size, source: MASTER) // Get portion
      bitonic(0, size, array, i % 2) // Sort incre/desc based on worker index

      MPI_SEND(buf: array, dest: MASTER) // Send sorted sub array to master
      


```
### 2c. Evaluation plan - what and how will you measure and compare
- Evaluating with multiple process counts, the total process count should always be a power of 2 (2^n processors):
  - Processor count: 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024
- Using Caliper + thicket to calculate the total time taken, time per process, etc, for each algorithm with the same inputs and processor count across each
  - This method can allow us to determine which algorithms are fastest in what input context
- Adjust the following in each evaluation case to test each algorithm:
  - Input sizes, Input types:
    - Input sizes: 2^16, 2^18, 2^20, 2^22, 2^24, 2^26, 2^28
    - Input types: Sorted, Random, Reverse sorted, 1% perturbed
  - Strong scaling (same problem size, increase the number of processors/nodes)
  - Weak scaling (increase problem size, increase the number of processors)

### 3a. Caliper instrumentation

Radix Sort Call Tree
```
0.030 main
├─ 0.000 data-init-runtime
├─ 0.017 comm
│  ├─ 0.002 comm-small
│  │  ├─ 0.001 MPI_Allreduce
│  │  └─ 0.000 MPI_Allgather
│  └─ 0.013 comm-large
│     ├─ 0.001 MPI_Send
│     └─ 0.010 MPI_Recv
├─ 0.005 comp
│  ├─ 0.002 comp-small
│  └─ 0.000 comp-large
├─ 0.006 MPI_Barrier
└─ 0.000 correctness-check
   ├─ 0.000 MPI_Send
   ├─ 0.000 MPI_Recv
   ├─ 0.000 MPI_Allreduce
   └─ 0.000 MPI_Barrier
0.000 MPI_Finalize
0.000 MPI_Initialized
0.000 MPI_Finalized
0.001 MPI_Comm_dup
```

Sample Sort Call Tree
```
1.099 main
└─ 1.099 main
   ├─ 0.718 MPI_Init
   ├─ 0.000 MPI_Comm_rank
   ├─ 0.000 MPI_Comm_size
   ├─ 0.001 data_init_runtime
   ├─ 0.112 comp
   │  ├─ 0.019 comp_small
   │  └─ 0.092 comp_large
   ├─ 0.263 comm
   │  ├─ 0.257 comm_large
   │  │  ├─ 0.016 MPI_Gather
   │  │  └─ 0.241 MPI_Bcast
   │  └─ 0.006 comm_small
   │     ├─ 0.000 MPI_Alltoall
   │     ├─ 0.002 MPI_Alltoallv
   │     └─ 0.003 MPI_Gatherv
   └─ 0.001 correctness_check
```

Merge Sort Call Tree
```
content here
```

Column Sort Call Tree
```
content here
```

Bitonic Sort Call Tree
```
content here
```

### 3b. Collect Metadata

Radix Sort Metadata

```
           cali.caliper.version  mpi.world.size  \
profile                                           
3298923211               2.11.0               4   
                                                 spot.metrics  \
profile                                                         
3298923211  min#inclusive#sum#time.duration,max#inclusive#...   
          spot.timeseries.metrics  spot.format.version  \
profile                                                   
3298923211                                            2   

                                                 spot.options  spot.channels  \
profile                                                                        
3298923211  time.variance,profile.mpi,node.order,region.co...  regionprofile   

           cali.channel spot:node.order spot:output spot:profile.mpi  \
profile                                                                
3298923211         spot            true   p4-a.cali             true   

           spot:region.count spot:time.exclusive spot:time.variance  \
profile                                                               
3298923211              true                true               true   

           		 launchdate                     libraries  \
profile                                                                     
3298923211  1729089290  [/scratch/group/csce435-f24/Caliper/caliper/li...   

                        cmdline cluster algorithm programming_model data_type  \
profile                                                                         
3298923211  [./radix_sort, 256]       c     radix               mpi       int   

            size_of_data_type  input_size input_type  num_procs scalability  \
profile                                                                       
3298923211                  4         256     Random          4      strong   

            group_num implementation_source  
profile                                      
3298923211         26           online & ai  
```

Sample Sort Metadata
```
         cali.caliper.version  mpi.world.size  \
profile                                          
743475682               2.11.0              16   

                                                spot.metrics  \
profile                                                        
743475682  min#inclusive#sum#time.duration,max#inclusive#...   

          spot.timeseries.metrics  spot.format.version  \
profile                                                  
743475682                                            2   

                                     spot.options  spot.channels cali.channel  \
profile                                                                         
743475682  node.order,region.count,time.exclusive  regionprofile         spot   

          spot:node.order                      spot:output spot:region.count  \
profile                                                                        
743475682            true  cali-samp-262144-p16-type1.cali              true   

          spot:time.exclusive  launchdate  \
profile                                     
743475682                true  1729115470   

                                                   libraries  \
profile                                                        
743475682  [/scratch/group/csce435-f24/Caliper/caliper/li...   

                         cmdline cluster   algorithm programming_model  \
profile                                                                  
743475682  [./samplesort, 18, 1]       c  samplesort               mpi   

          data_type  size_of_data_type  input_size input_type  num_procs  \
profile                                                                    
743475682       int                  4      262144     Random         16   

          scalability  group_num  \
profile                            
743475682      strong         26   

                                       implementation_source  
profile                                                       
743475682  AI (ChatGPT) and Online (Class Notes, https://www.geeksforgeeks.org/cpp-program-for-quicksort/, and https://en.wikipedia.org/wiki/Samplesort#:~:text=sequential%2C%20sorting%20algorithm.-,Pseudocode,-%5Bedit%5D)  
```

Merge Sort Metadata
```
content here
```

Column Sort Metadata
```
content here
```

Bitonic Sort Metadata
```
content here
```
