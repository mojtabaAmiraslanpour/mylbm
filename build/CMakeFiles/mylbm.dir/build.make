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
include CMakeFiles/mylbm.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/mylbm.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/mylbm.dir/flags.make

CMakeFiles/mylbm.dir/main.cpp.o: CMakeFiles/mylbm.dir/flags.make
CMakeFiles/mylbm.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/media/mojtaba/data/MyCode/mylbm/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/mylbm.dir/main.cpp.o"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mylbm.dir/main.cpp.o -c /media/mojtaba/data/MyCode/mylbm/main.cpp

CMakeFiles/mylbm.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mylbm.dir/main.cpp.i"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /media/mojtaba/data/MyCode/mylbm/main.cpp > CMakeFiles/mylbm.dir/main.cpp.i

CMakeFiles/mylbm.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mylbm.dir/main.cpp.s"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /media/mojtaba/data/MyCode/mylbm/main.cpp -o CMakeFiles/mylbm.dir/main.cpp.s

# Object files for target mylbm
mylbm_OBJECTS = \
"CMakeFiles/mylbm.dir/main.cpp.o"

# External object files for target mylbm
mylbm_EXTERNAL_OBJECTS =

mylbm: CMakeFiles/mylbm.dir/main.cpp.o
mylbm: CMakeFiles/mylbm.dir/build.make
mylbm: CMakeFiles/mylbm.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/media/mojtaba/data/MyCode/mylbm/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable mylbm"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mylbm.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/mylbm.dir/build: mylbm

.PHONY : CMakeFiles/mylbm.dir/build

CMakeFiles/mylbm.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/mylbm.dir/cmake_clean.cmake
.PHONY : CMakeFiles/mylbm.dir/clean

CMakeFiles/mylbm.dir/depend:
	cd /media/mojtaba/data/MyCode/mylbm/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /media/mojtaba/data/MyCode/mylbm /media/mojtaba/data/MyCode/mylbm /media/mojtaba/data/MyCode/mylbm/build /media/mojtaba/data/MyCode/mylbm/build /media/mojtaba/data/MyCode/mylbm/build/CMakeFiles/mylbm.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/mylbm.dir/depend

