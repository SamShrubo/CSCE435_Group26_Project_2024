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
0.889 main
├─ 0.867 main
│  ├─ 0.864 data_init
│  │  └─ 0.000 MPI_Init
│  ├─ 0.002 comp
│  │  ├─ 0.001 comp_small
│  │  └─ 0.001 comp_large
│  ├─ 0.002 comm
│  │  └─ 0.002 comm_large
│  │     ├─ 0.003 MPI_Recv
│  │     └─ 0.000 MPI_Send
│  └─ 0.000 correctness_check
├─ 0.000 MPI_Finalize
├─ 0.000 MPI_Initialized
├─ 0.000 MPI_Finalized
└─ 0.022 MPI_Comm_dup
```

Column Sort Call Tree
```
0.614 main
└─ 0.614 main
   ├─ 0.002 create_matrix
   ├─ 0.000 gather_final
   ├─ 0.000 step1_sort
   ├─ 0.001 step2_transpose
   ├─ 0.000 step3_sort
   ├─ 0.002 step4_transpose
   ├─ 0.000 step5_sort
   ├─ 0.000 step6_shift
   ├─ 0.000 step7_sort
   └─ 0.000 step8_unshift
```

Bitonic Sort Call Tree
```
332.070 main
├─ 0.000 MPI_Init
├─ 331.693 Whole Computation
│  ├─ 281.564 comm
│  │  ├─ 0.000 comm_small
│  │  │  └─ 0.000 MPI_Send
│  │  ├─ 201.302 comp
│  │  │  └─ 71.301 comp_large
│  │  ├─ 0.002 correctness_check
│  │  └─ 21.190 comm_large
│  │     └─ 21.190 MPI_Recv
│  ├─ 0.000 MPI_Recv
│  ├─ 0.006 data_init_runtime
│  ├─ 201.979 comp
│  │  └─ 71.532 comp_large
│  ├─ 0.004 correctness_check
│  └─ 0.001 MPI_Send
├─ 0.000 MPI_Initialized
├─ 0.000 MPI_Finalized
└─ 0.000 MPI_Comm_dup
0.000 MPI_Initialized
0.000 MPI_Finalized
249.995 MPI_Comm_dup
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
           cali.caliper.version  mpi.world.size  \
profile                                           
1981606483               2.11.0              16   

                                                 spot.metrics  \
profile                                                         
1981606483  min#inclusive#sum#time.duration,max#inclusive#...   

           spot.timeseries.metrics  spot.format.version  \
profile                                                   
1981606483                                            2   

                                                 spot.options  spot.channels  \
profile                                                                        
1981606483  time.variance,profile.mpi,node.order,region.co...  regionprofile   

           cali.channel spot:node.order   spot:output spot:profile.mpi  \
profile                                                                  
1981606483         spot            true  p16-a16.cali             true   

           spot:region.count spot:time.exclusive spot:time.variance  \
profile                                                               
1981606483              true                true               true   

            launchdate                                          libraries  \
profile                                                                     
1981606483  1729133312  [/scratch/group/csce435-f24/Caliper/caliper/li...   

                              cmdline cluster   algorithm programming_model  \
profile                                                                       
1981606483  [./mergesort, 16, Random]       c  Merge Sort               mpi   

           data_type  size_of_data_type  input_size input_type  num_procs  \
profile                                                                     
1981606483       int                  4       65536     Random         16   

                  scalability  group_num implementation_source  
profile                                                         
1981606483  weak (assumption)         26      ai & handwritten  
```

Column Sort Metadata
```
cali.caliper.version  mpi.world.size  \
profile                                           
202565079                2.11.0              10   
977574015                2.11.0               3   
3024513808               2.11.0              10   
3635606938               2.11.0               2   
4276463802               2.11.0              10   

                                                 spot.metrics  \
profile                                                         
202565079   min#inclusive#sum#time.duration,max#inclusive#...   
977574015   min#inclusive#sum#time.duration,max#inclusive#...   
3024513808  min#inclusive#sum#time.duration,max#inclusive#...   
3635606938  min#inclusive#sum#time.duration,max#inclusive#...   
4276463802  min#inclusive#sum#time.duration,max#inclusive#...   

           spot.timeseries.metrics  spot.format.version  \
profile                                                   
202565079                                             2   
977574015                                             2   
3024513808                                            2   
3635606938                                            2   
4276463802                                            2   

                                      spot.options  spot.channels  \
profile                                                             
202565079   node.order,region.count,time.exclusive  regionprofile   
977574015   node.order,region.count,time.exclusive  regionprofile   
3024513808  node.order,region.count,time.exclusive  regionprofile   
3635606938  node.order,region.count,time.exclusive  regionprofile   
4276463802  node.order,region.count,time.exclusive  regionprofile   

           cali.channel spot:node.order       spot:output spot:region.count  \
profile                                                                       
202565079          spot            true    p10-a1000.cali              true   
977574015          spot            true       p3-a11.cali              true   
3024513808         spot            true   p10-a10000.cali              true   
3635606938         spot            true        p2-a4.cali              true   
4276463802         spot            true  p10-a100000.cali              true   

           spot:time.exclusive  
profile                         
202565079                 true  
977574015                 true  
3024513808                true  
3635606938                true  
4276463802                true  
```

Bitonic Sort Metadata
```
 	cali.caliper.version 	mpi.world.size 	spot.metrics 	spot.timeseries.metrics 	spot.format.version 	spot.options 	spot.channels 	cali.channel 	spot:node.order 	spot:output 	spot:region.count 	spot:time.exclusive 	user 	launchdate 	libraries 	cmdline 	cluster 	num_procs 	matrix_size 	program_name 	matrix_datatype_size 	MPI_Reduce-whole_computation_time 	MPI_Reduce-master_initialization_time 	MPI_Reduce-master_send_receive_time 	MPI_Reduce-worker_receive_time_max 	MPI_Reduce-worker_receive_time_min 	MPI_Reduce-worker_receive_time_average 	MPI_Reduce-worker_calculation_time_max 	MPI_Reduce-worker_calculation_time_min 	MPI_Reduce-worker_calculation_time_average 	MPI_Reduce-worker_send_time_max 	MPI_Reduce-worker_send_time_min 	MPI_Reduce-worker_send_time_average
profile 																																	
170336558 	2.11.0 	8 	min#inclusive#sum#time.duration,max#inclusive#... 		2 	node.order,region.count,time.exclusive 	regionprofile 	spot 	true 	128-8.cali 	true 	true 	zhongyouwu 	1727466897 	[/scratch/group/csce435-f24/Caliper/caliper/li... 	[./mpi_mm, 128] 	c 	8 	128 	master_worker_matrix_multiplication 	8 	0.024461 	0.000229 	0.024208 	0.055580 	0.023114 	0.041921 	0.011459 	0.001422 	0.003019 	0.000048 	0.000031 	0.000040
458074348 	2.11.0 	16 	min#inclusive#sum#time.duration,max#inclusive#... 		2 	node.order,region.count,time.exclusive 	regionprofile 	spot 	true 	128-16.cali 	true 	true 	zhongyouwu 	1727466897 	[/scratch/group/csce435-f24/Caliper/caliper/li... 	[./mpi_mm, 128] 	c 	16 	128 	master_worker_matrix_multiplication 	8 	0.077896 	0.000243 	0.077624 	0.123626 	0.000429 	0.074863 	0.000741 	0.000642 	0.000696 	0.000048 	0.000035 	0.000041
802824523 	2.11.0 	16 	min#inclusive#sum#time.duration,max#inclusive#... 		2 	node.order,region.count,time.exclusive 	regionprofile 	spot 	true 	1024-16.cali 	true 	true 	zhongyouwu 	1727466954 	[/scratch/group/csce435-f24/Caliper/caliper/li... 	[./mpi_mm, 1024] 	c 	16 	1024 	master_worker_matrix_multiplication 	8 	5.962221 	0.112963 	5.849179 	0.525303 	0.206468 	0.364798 	5.638342 	5.240062 	5.426046 	0.030461 	0.000195 	0.003998
1325292579 	2.11.0 	32 	min#inclusive#sum#time.duration,max#inclusive#... 		2 	node.order,region.count,time.exclusive 	regionprofile 	spot 	true 	1024-32.cali 	true 	true 	zhongyouwu 	1727466954 	[/scratch/group/csce435-f24/Caliper/caliper/li... 	[./mpi_mm, 1024] 	c 	32 	1024 	master_worker_matrix_multiplication 	8 	8.474886 	0.220471 	8.254303 	1.513168 	0.375166 	1.002473 	7.094080 	5.674065 	6.319820 	0.066136 	0.000126 	0.007246
2875060761 	2.11.0 	64 	min#inclusive#sum#time.duration,max#inclusive#... 		2 	node.order,region.count,time.exclusive 	regionprofile 	spot 	true 	128-64.cali 	true 	true 	zhongyouwu 	1727466898 	[/scratch/group/csce435-f24/Caliper/caliper/li... 	[./mpi_mm, 128] 	c 	64 	128 	master_worker_matrix_multiplication 	8 	0.193155 	0.000243 	0.192825 	0.250448 	0.000171 	0.085454 	0.014162 	0.000168 	0.000631 	0.000063 	0.000028 	0.000033
2896985760 	2.11.0 	2 	min#inclusive#sum#time.duration,max#inclusive#... 		2 	node.order,region.count,time.exclusive 	regionprofile 	spot 	true 	1024-2.cali 	true 	true 	zhongyouwu 	1727466897 	[/scratch/group/csce435-f24/Caliper/caliper/li... 	[./mpi_mm, 1024] 	c 	2 	1024 	master_worker_matrix_multiplication 	8 	10.373385 	0.010415 	10.362917 	0.021158 	0.021158 	0.021158 	10.351604 	10.351604 	10.351604 	0.003412 	0.003412 	0.003412
3146970658 	2.11.0 	2 	min#inclusive#sum#time.duration,max#inclusive#... 		2 	node.order,region.count,time.exclusive 	regionprofile 	spot 	true 	128-2.cali 	true 	true 	zhongyouwu 	1727465561 	[/scratch/group/csce435-f24/Caliper/caliper/li... 	[./mpi_mm, 128] 	c 	2 	128 	master_worker_matrix_multiplication 	8 	0.019298 	0.000279 	0.018984 	0.010334 	0.010334 	0.010334 	0.010818 	0.010818 	0.010818 	0.000068 	0.000068 	0.000068
3258936128 	2.11.0 	32 	min#inclusive#sum#time.duration,max#inclusive#... 		2 	node.order,region.count,time.exclusive 	regionprofile 	spot 	true 	128-32.cali 	true 	true 	zhongyouwu 	1727466897 	[/scratch/group/csce435-f24/Caliper/caliper/li... 	[./mpi_mm, 128] 	c 	32 	128 	master_worker_matrix_multiplication 	8 	0.217835 	0.000286 	0.217514 	0.358161 	0.146064 	0.253260 	0.000667 	0.000332 	0.000360 	0.001665 	0.000029 	0.000090
3263988617 	2.11.0 	64 	min#inclusive#sum#time.duration,max#inclusive#... 		2 	node.order,region.count,time.exclusive 	regionprofile 	spot 	true 	1024-64.cali 	true 	true 	zhongyouwu 	1727466955 	[/scratch/group/csce435-f24/Caliper/caliper/li... 	[./mpi_mm, 1024] 	c 	64 	1024 	master_worker_matrix_multiplication 	8 	12.363227 	0.318166 	12.044952 	4.062414 	0.473372 	2.416818 	10.249945 	6.819659 	8.160158 	0.117756 	0.000096 	0.005823
3434837044 	2.11.0 	4 	min#inclusive#sum#time.duration,max#inclusive#... 		2 	node.order,region.count,time.exclusive 	regionprofile 	spot 	true 	128-4.cali 	true 	true 	zhongyouwu 	1727466846 	[/scratch/group/csce435-f24/Caliper/caliper/li... 	[./mpi_mm, 128] 	c 	4 	128 	master_worker_matrix_multiplication 	8 	0.026028 	0.000254 	0.025743 	0.033771 	0.015406 	0.027525 	0.003430 	0.003333 	0.003391 	0.000077 	0.000051 	0.000063
4202622048 	2.11.0 	8 	min#inclusive#sum#time.duration,max#inclusive#... 		2 	node.order,region.count,time.exclusive 	regionprofile 	spot 	true 	1024-8.cali 	true 	true 	zhongyouwu 	1727466954 	[/scratch/group/csce435-f24/Caliper/caliper/li... 	[./mpi_mm, 1024] 	c 	8 	1024 	master_worker_matrix_multiplication 	8 	6.465010 	0.040723 	6.424182 	0.192274 	0.076241 	0.152816 	6.308338 	6.006143 	6.124993 	0.024461 	0.000866 	0.006158
4288893745 	2.11.0 	4 	min#inclusive#sum#time.duration,max#inclusive#... 		2 	node.order,region.count,time.exclusive 	regionprofile 	spot 	true 	1024-4.cali 	true 	true 	zhongyouwu 	1727466897 	[/scratch/group/csce435-f24/Caliper/caliper/li... 	[./mpi_mm, 1024] 	c 	4 	1024 	master_worker_matrix_multiplication 	8 	7.564211 	0.031770 	7.532343 	0.093590 	0.067238 	0.076938 	7.483138 	7.354114 	7.410794 	0.012188 	0.002114 	0.008659
```
### 4a. Vary the following parameters
For input_size's:
- 2^16, 2^18, 2^20, 2^22, 2^24, 2^26, 2^28

For input_type's:
- Sorted, Random, Reverse sorted, 1%perturbed

MPI: num_procs:
- 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024

This should result in 4x7x10=280 Caliper files for your MPI experiments.

### 4b. Hints for performance analysis

To automate running a set of experiments, parameterize your program.

- input_type: "Sorted" could generate a sorted input to pass into your algorithms
- algorithm: You can have a switch statement that calls the different algorithms and sets the Adiak variables accordingly
- num_procs: How many MPI ranks you are using

When your program works with these parameters, you can write a shell script 
that will run a for loop over the parameters above (e.g., on 64 processors, 
perform runs that invoke algorithm2 for Sorted, ReverseSorted, and Random data).  

### 4c. You should measure the following performance metrics
- `Time`
    - Min time/rank
    - Max time/rank
    - Avg time/rank
    - Total time
    - Variance time/rank
 
#### Radix Sort Performance Evaluation:
Hydra errors in the grace cluster made it borderline impossible to consistently run 512 and 1024 processes which is why I only have one 512 run. The code becomes to slow for 2^26 elements at lower processor numbers when it takes about 2 hours for one run at 64 processors. The reason for the delay is probably that there is a lot of sending happening that could probably be optimized but I don’t know how to optimize it yet so I can’t run the largest problem size. 
##### Graphs 
![alt text](Graphs/main-65536.png)
![alt text](Graphs/main-262144.png)
![alt text](Graphs/main-1048576.png)
![alt text](Graphs/main-4194304.png)
![alt text](Graphs/main-16777216.png)
![alt text](Graphs/main-67108864.png)
There is a general downward trend for all the graph and almost all of them have a spike when reaching 32 processors. The decreasing trend can be attributed to the local array shriking and the proccors having the communicate less when the number of processors goes up since the problem size stays the same. 
![alt text](Graphs/comm-65536.png)
![alt text](Graphs/comm-262144.png)
![alt text](Graphs/comm-1048576.png)
![alt text](Graphs/comm-4194304.png)
![alt text](Graphs/comm-16777216.png)
![alt text](Graphs/comm-67108864.png)
There is a general downward trend for all the graph but other than that there isn't much similarity across all the graphs. The decreasing trend can be attributed to the local array shriking and the proccors having the communicate less when the number of processors goes up since the problem size stays the same. 
![alt text](Graphs/comp-large-65536.png)
![alt text](Graphs/comp-large-262144.png)
![alt text](Graphs/comp-large-1048576.png)
![alt text](Graphs/comp-large-4194304.png)
![alt text](Graphs/comp-large-16777216.png)
![alt text](Graphs/comp-large-67108864.png)
There is very little change in the speed up for computation as the number of processros increase and this is most likely due to the way the caliper barriers were placed during implemenation as there were other place it should have been placed but I didn't notice it until much later. I plan on rerunning my code to get more accurate timings on my algorithim. 

#### Sample Sort Performance Evaluation
I had issues running jobs with 512 and 1024 processors on Grace. Hydra consistently would error, and the runs would not complete. Additionally, Grace's exceptionally long queue times on 10/21 and 10/22 made it impossible for me to run more jobs. This is why the graphs for 2^24, 2^26, and 2^28 are sparse. 


##### 2^16 array element graphs
![Alt text](Graphs/samplesort_main_65536.png)
![Alt text](Graphs/samplesort_comp_65536.png)
![Alt text](Graphs/samplesort_comm_65536.png)
- Computation time decreases as the number of processors increases while the computation time increases. Overall, the general pattern of runtimes is that the runtime decreases as the number of processors increases, but there is a point of dimininishing returns at smaller input sizes. 

##### 2^18 array element graphs
![Alt text](Graphs/samplesort_main_262144.png)
![Alt text](Graphs/samplesort_comp_262144.png)
![Alt text](Graphs/samplesort_comm_262144.png)
- Computation time decreases as the number of processors increases while the computation time increases. Communication had strange spikes for the 1% perturbed and reverse sorted inputs. The sorted input had a consistently low communication time due to array objects not having to be sent between buckets. 

##### 2^20 array element graphs
![Alt text](Graphs/samplesort_main_1048576.png)
![Alt text](Graphs/samplesort_comp_1048576.png)
![Alt text](Graphs/samplesort_comm_1048576.png)
- Computation and communication times generally decrease as the number of processors increase. Overall, the general pattern of runtimes is that the runtime decreases as the number of processors increases. All data in these graphs looks as expected. 

##### 2^22 array element graphs
![Alt text](Graphs/samplesort_main_4194304.png)
![Alt text](Graphs/samplesort_comp_4194304.png)
![Alt text](Graphs/samplesort_comm_4194304.png)
- There are holes in this graph due to the excessive time that the 2 and 4 processor jobs take with 2^22 array elements as inputs. In general, computation and communication times generally decrease as the number of processors increase.

##### 2^24 array element graphs
![Alt text](Graphs/samplesort_main_16777216.png)
![Alt text](Graphs/samplesort_comp_16777216.png)
![Alt text](Graphs/samplesort_comm_16777216.png)
- Performance time appears to decrease as the number of processors increases. More time is spent communicating than performing computations. 

##### 2^26 array element graphs
![Alt text](Graphs/samplesort_main_67108864.png)
![Alt text](Graphs/samplesort_comp_67108864.png)
![Alt text](Graphs/samplesort_comm_67108864.png)
- The only job that has been run with a 2^26 array input size was the 1024 process job. There is no other data to compare this point against. Comparatively, the most time was spent communicating, not performing computations. 

##### 2^28 array element graphs
- At this point in time, no jobs have been run with 2^28 array elements in them. This is due to the current implementation of samplesort having poor scaling and timing out within a reasonable amount of time (four hours).

##### Comments
In general, the current implementation of sample sort appears to not scale well due to the amount of time spent communicating. Performance time also almost always appears to decrease as the number of processors increases. I intend on reconfiguring my code to reduce communication times and rerunning all of my jobs again in the next week to have more comprehensive data.  

#### Merge Sort Performance Evaluation:
There were many issues at the larger processor counts of 512 and 1024, consistent errors resulting in a lack of data for many points within 512 processors and all in 1024, this is most likely due to the excessive communication loop required to make merge sort work and it lacks scalability to higher processor counts.

##### 2^16 array element graphs
![Alt text](Graphs/ms_main_16.png)
![Alt text](Graphs/ms_comm_16.png)
![Alt text](Graphs/ms_comp_16.png)
- Here computation time has a negligible influence on the total run time, meaning the algorithm is far less efficient with more processors at a higher run time due to higher communication overhead being the driving factor.
##### 2^18 array element graphs
![Alt text](Graphs/ms_main_18.png)
![Alt text](Graphs/ms_comm_18.png)
![Alt text](Graphs/ms_comp_18.png)
- The same trend is observed at this array size as in 2^16
##### 2^20 array element graphs
![Alt text](Graphs/ms_main_20.png)
![Alt text](Graphs/ms_comm_20.png)
![Alt text](Graphs/ms_comp_20.png)
- The same trend is observed at this array size as in 2^16
##### 2^22 array element graphs
![Alt text](Graphs/ms_main_22.png)
![Alt text](Graphs/ms_comm_22.png)
![Alt text](Graphs/ms_comp_22.png)
- The same trend is observed at this array size as in 2^16, though a slight performance gain can be found at around 4-8 processors
##### 2^24 array element graphs
![Alt text](Graphs/ms_main_24.png)
![Alt text](Graphs/ms_comm_24.png)
![Alt text](Graphs/ms_comp_24.png)
- An exponential decay is more prevalent here and optimal time is around 16-32 procs
##### 2^26 array element graphs
![Alt text](Graphs/ms_main_26.png)
![Alt text](Graphs/ms_comm_26.png)
![Alt text](Graphs/ms_comp_26.png)
- A strong exponential decay can be observed here showing that computation is the driving factor and significantly beneifts from the smaller array chunks with more processors
##### 2^28 array element graphs
![Alt text](Graphs/ms_main_28.png)
![Alt text](Graphs/ms_comm_28.png)
![Alt text](Graphs/ms_comp_28.png)
- Same trend as 2^26

##### Comments:
It is clear that at the larger array sizes the driving factor in run time of the program is mainly influenced by computation time, as computation time makes up a far higher proportion at the larger array sizes, but for smaller arrays the time benefit is negated by communication overhead and other smaller factors.

For the communication there is a strange spike at 32 processors for some of the larger arrays, I would presume this has something to do with the communication at that point in the node or possibly a limit on memory usage, however the general communication numbers are so small that the ultimate impact this spike has on total run time of the program is frankly negligible.

For the computation time it would appear that the sorted array takes consistently the longest time, i presume this has to do with the local sort and that my code does not have a catch to auto-pass an array if it is sorted, meaning it most likely goes through the entire process to sort the array resulting in a slower average time.

#### Column Sort Performance Evaluation:
```
Content to be added later
```

#### Bitonic Sort Performance Evaluation:
```
Content to be added later
```
