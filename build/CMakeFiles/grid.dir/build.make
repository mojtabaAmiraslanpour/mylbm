# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.18

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /media/mojtaba/data/MyCode/mylbm

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /media/mojtaba/data/MyCode/mylbm/build

# Include any dependencies generated for this target.
include CMakeFiles/grid.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/grid.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/grid.dir/flags.make

CMakeFiles/grid.dir/grid.C.o: CMakeFiles/grid.dir/flags.make
CMakeFiles/grid.dir/grid.C.o: ../grid.C
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/media/mojtaba/data/MyCode/mylbm/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/grid.dir/grid.C.o"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/grid.dir/grid.C.o -c /media/mojtaba/data/MyCode/mylbm/grid.C

CMakeFiles/grid.dir/grid.C.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/grid.dir/grid.C.i"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /media/mojtaba/data/MyCode/mylbm/grid.C > CMakeFiles/grid.dir/grid.C.i

CMakeFiles/grid.dir/grid.C.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/grid.dir/grid.C.s"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /media/mojtaba/data/MyCode/mylbm/grid.C -o CMakeFiles/grid.dir/grid.C.s

# Object files for target grid
grid_OBJECTS = \
"CMakeFiles/grid.dir/grid.C.o"

# External object files for target grid
grid_EXTERNAL_OBJECTS =

libgrid.a: CMakeFiles/grid.dir/grid.C.o
libgrid.a: CMakeFiles/grid.dir/build.make
libgrid.a: CMakeFiles/grid.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/media/mojtaba/data/MyCode/mylbm/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libgrid.a"
	$(CMAKE_COMMAND) -P CMakeFiles/grid.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/grid.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/grid.dir/build: libgrid.a

.PHONY : CMakeFiles/grid.dir/build

CMakeFiles/grid.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/grid.dir/cmake_clean.cmake
.PHONY : CMakeFiles/grid.dir/clean

CMakeFiles/grid.dir/depend:
	cd /media/mojtaba/data/MyCode/mylbm/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /media/mojtaba/data/MyCode/mylbm /media/mojtaba/data/MyCode/mylbm /media/mojtaba/data/MyCode/mylbm/build /media/mojtaba/data/MyCode/mylbm/build /media/mojtaba/data/MyCode/mylbm/build/CMakeFiles/grid.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/grid.dir/depend
