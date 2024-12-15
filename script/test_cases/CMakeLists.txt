cmake_minimum_required (VERSION 3.1.1)

project(TEST_CASES CXX)

# Gathering source files
file(GLOB_RECURSE source_files ../../source/*)
file(GLOB_RECURSE test_files ../../script/test_cases/*)

# Include

# DACE
include_directories(../../../lib/dace/interfaces/cxx/include)
include_directories(../../../lib/dace/build/core/include)
include_directories(../../../lib/dace/core/include)

# Local files
include_directories(../../include)
include_directories(../../source)
include_directories(../../script/test_cases)
include_directories(../../script)

# Linking librairies (DACE)
link_directories(../../../lib/dace/build/interfaces/cxx)

# Main executable
add_executable(
	test_cases
	${source_files} 
	${test_files}
	../../script/test_cases/run_test_cases.cpp)

# Properties
set_target_properties(test_cases PROPERTIES
    CXX_STANDARD 14)
	
# Link
target_link_libraries(
    test_cases
    dace)
