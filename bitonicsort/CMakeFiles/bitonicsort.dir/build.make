# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.12

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /sw/eb/sw/CMake/3.12.1/bin/cmake

# The command to remove a file.
RM = /sw/eb/sw/CMake/3.12.1/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /scratch/user/zhongyouwu/Project/MPI_Builds

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /scratch/user/zhongyouwu/Project/MPI_Builds

# Include any dependencies generated for this target.
include CMakeFiles/bitonicsort.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/bitonicsort.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/bitonicsort.dir/flags.make

CMakeFiles/bitonicsort.dir/Bitonic.cpp.o: CMakeFiles/bitonicsort.dir/flags.make
CMakeFiles/bitonicsort.dir/Bitonic.cpp.o: Bitonic.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/scratch/user/zhongyouwu/Project/MPI_Builds/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/bitonicsort.dir/Bitonic.cpp.o"
	/sw/eb/sw/GCCcore/8.3.0/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bitonicsort.dir/Bitonic.cpp.o -c /scratch/user/zhongyouwu/Project/MPI_Builds/Bitonic.cpp

CMakeFiles/bitonicsort.dir/Bitonic.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bitonicsort.dir/Bitonic.cpp.i"
	/sw/eb/sw/GCCcore/8.3.0/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /scratch/user/zhongyouwu/Project/MPI_Builds/Bitonic.cpp > CMakeFiles/bitonicsort.dir/Bitonic.cpp.i

CMakeFiles/bitonicsort.dir/Bitonic.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bitonicsort.dir/Bitonic.cpp.s"
	/sw/eb/sw/GCCcore/8.3.0/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /scratch/user/zhongyouwu/Project/MPI_Builds/Bitonic.cpp -o CMakeFiles/bitonicsort.dir/Bitonic.cpp.s

# Object files for target bitonicsort
bitonicsort_OBJECTS = \
"CMakeFiles/bitonicsort.dir/Bitonic.cpp.o"

# External object files for target bitonicsort
bitonicsort_EXTERNAL_OBJECTS =

bitonicsort: CMakeFiles/bitonicsort.dir/Bitonic.cpp.o
bitonicsort: CMakeFiles/bitonicsort.dir/build.make
bitonicsort: /scratch/group/csce435-f24/Caliper/caliper/lib64/libcaliper.so.2.11.0
bitonicsort: /sw/eb/sw/impi/2019.9.304-iccifort-2020.4.304/intel64/lib/libmpicxx.so
bitonicsort: /sw/eb/sw/impi/2019.9.304-iccifort-2020.4.304/intel64/lib/release/libmpi.so
bitonicsort: /lib64/librt.so
bitonicsort: /lib64/libpthread.so
bitonicsort: /lib64/libdl.so
bitonicsort: CMakeFiles/bitonicsort.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/scratch/user/zhongyouwu/Project/MPI_Builds/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable bitonicsort"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/bitonicsort.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/bitonicsort.dir/build: bitonicsort

.PHONY : CMakeFiles/bitonicsort.dir/build

CMakeFiles/bitonicsort.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/bitonicsort.dir/cmake_clean.cmake
.PHONY : CMakeFiles/bitonicsort.dir/clean

CMakeFiles/bitonicsort.dir/depend:
	cd /scratch/user/zhongyouwu/Project/MPI_Builds && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /scratch/user/zhongyouwu/Project/MPI_Builds /scratch/user/zhongyouwu/Project/MPI_Builds /scratch/user/zhongyouwu/Project/MPI_Builds /scratch/user/zhongyouwu/Project/MPI_Builds /scratch/user/zhongyouwu/Project/MPI_Builds/CMakeFiles/bitonicsort.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/bitonicsort.dir/depend
