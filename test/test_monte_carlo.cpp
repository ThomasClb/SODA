/**
	test_monte_carlo.cpp

	Purpose: Test of the Monte-Carlo validation methods.
    
	@author Thomas Caleb

	@version 1.0 13/12/2024
*/

#include "gtest/gtest.h"
#include "../include/soda.h"
#include "../include/monte_carlo.h"

using namespace DACE;
using namespace std;

TEST(TestMC, GenerateNormalSample) {

	size_t size_sample=100;
	vectordb mean{-1, -5, -3};
	matrixdb Sigma(mean.size(), mean.size(), 0.0);
	for (size_t i=0; i<mean.size(); i++) {
		Sigma.at(i, i) = i*1.0 + 1.0;
	}
	matrixdb sample = generate_normal_sample(size_sample, mean, Sigma);

	// Tests
	EXPECT_EQ(sample.ncols(), size_sample);
	EXPECT_EQ(sample.nrows(), mean.size());
}
