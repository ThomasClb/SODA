cmake_minimum_required (VERSION 3.1.1)

project(UNIT_TEST CXX)

# Gathering source files
file(GLOB_RECURSE include_files ../include/*)
file(GLOB_RECURSE source_files ../source/*)
file(GLOB_RECURSE test_files ../test/*)

# Find gtest
find_package(GTest CONFIG REQUIRED)

# Include

# DACE
include_directories(../../lib/dace/interfaces/cxx/include)
include_directories(../../lib/dace/build/core/include)
include_directories(../../lib/dace/core/include)

# Local files
include_directories(../test)
include_directories(../include)

# Linking librairies
link_directories(../../lib/dace/build/interfaces/cxx)

# Main executable
add_executable(
	unit_test
	${source_files}
	${include_files}
	${test_files}
	./run_test.cpp)

# Properties
set_target_properties(unit_test PROPERTIES
    CXX_STANDARD 14)
	
# Link
target_link_libraries(
    unit_test
    dace
    gtest
    gtest_main)
