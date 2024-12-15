/**
	test_IO.cpp

	Purpose: Test of the implementation of the inputs and outputs of data.
    methods.
    
	@author Thomas Caleb

	@version 1.0 13/12/2024
*/

#include "gtest/gtest.h"
#include "../include/soda.h"

using namespace DACE;
using namespace std;

TEST(TestIO, Split) {
	// Init
	string test="BLA,BLE,BLI,BLO,BLU,BLY";
	string delimiter=",";
	vector<string> list_split=split(test, delimiter);

	// Test
	EXPECT_EQ(list_split.size(), 6);
}
