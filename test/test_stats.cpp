/**
	test_stats.cpp

	Purpose: Test of the implementation of the statistical functions.
    
	@author Thomas Caleb

	@version 1.0 13/12/2024
*/

#include "gtest/gtest.h"
#include "../include/soda.h"

using namespace DACE;
using namespace std;

TEST(TestStats, CompareDA) {

	// Init DACE
	size_t nb_variables = 3;
	DA::init(2, nb_variables);

	// Init
	DA a = 1.0 + DA(1);
	DA b = 2.0 + DA(2);
	DA c = -1.5 + DA(3);

	// Test
	EXPECT_TRUE(compare_DA(a, b));
	EXPECT_TRUE(compare_DA(c, a));
}
TEST(TestStats, SortVector) {

	// Init DACE
	size_t nb_variables = 3;
	DA::init(2, nb_variables);

	// Init
	DA a = 1.0 + DA(1);
	DA b = 2.0 + DA(2);
	DA c = -1.5 + DA(3);
	vectorDA vect_DA{a, b, c};
	vectordb vect = vect_DA.cons();
	vect_DA = sort_vector(vect_DA);
	vect = sort_vector(vect);

	// Test
	for (size_t i=0; i<vect.size() - 1; i++) {
		EXPECT_TRUE(compare_DA(vect_DA[i], vect_DA[i+1]));
		EXPECT_TRUE(vect[i] < vect[i+1]);
	}
}
TEST(TestStats, IncBeta) {
	// Init
	double d = 5.0;
	double a = (d - 1)/2;

	// Test
	EXPECT_EQ(inc_beta(a, 0.0), 0.0);
	EXPECT_EQ(inc_beta(a, 1.0), 1.0);
}
TEST(TestStats, Chi2CDF) {
	// Init
	double d = 5.0;

	// Test
	EXPECT_EQ(chi_2_cdf(d, 0.0), 0.0);
	EXPECT_EQ(chi_2_cdf(d, 1e13), 1.0);
}
TEST(TestStats, InvChi2CDF) {
	// Init
	double d = 5.0;
	double r = 1.0;

	// Test
	EXPECT_NEAR(inv_chi_2_cdf(d, 0.0), 0.0, 1e-8);
	EXPECT_NEAR(inv_chi_2_cdf(d, chi_2_cdf(d, r)), r, 1e-8);
}
TEST(TestStats, DeltaPhi) {
	// Init
	double d = 5.0;
	double r_1 = 1.0;
	double r_2 = 0.5;

	// Test
	EXPECT_EQ(delta_phi(d, r_1, r_1), 0.0);
	EXPECT_NEAR(delta_phi(d, 1/EPS, EPS), 1.0, EPS);
	EXPECT_TRUE(delta_phi(d, r_1, r_2) > 0.0);
}
TEST(TestStats, Conservatism) {
	// Init
	double beta = 1e-2;
	double beta_r_1 = 1e-1;
	double beta_r_2 = 1e-2;
	double beta_r_3 = 1e-3;

	// Test
	EXPECT_EQ(conservatism(beta, beta_r_2), 1.0);
	EXPECT_TRUE(conservatism(beta, beta_r_1) < 1);
	EXPECT_TRUE(conservatism(beta, beta_r_3) > 1);
}
TEST(TestStats, GetListDistance) {
	// Init
	vectordb mean{-1.0, -1.5};
	matrixdb Sigma(2, 2, 0.0);
	Sigma.at(0, 0) = 0.1;
	Sigma.at(1, 1) = 0.2;
	vectordb list_distance = get_list_distance(mean, get_diag_vector_(Sigma));

	// Test
	EXPECT_EQ(list_distance.size(), mean.size());
	for (size_t i=0; i<list_distance.size() - 1; i++) {
		EXPECT_TRUE(list_distance[i] < list_distance[i+1]);
	}
}
TEST(TestStats, FirstOrderRisk) {
	// Init
	vectordb mean{-1.0, -1.5};
	matrixdb Sigma(2, 2, 0.0);
	Sigma.at(0, 0) = 0.1;
	Sigma.at(1, 1) = 0.2;
	double beta = first_order_risk_estimation(mean, get_diag_vector_(Sigma));

	// Test
	EXPECT_TRUE(beta >= 0);
	EXPECT_TRUE(beta <= 1);
}
TEST(TestStats, DthOrderRisk) {
	// Init
	vectordb mean{-1.0, -1.5};
	matrixdb Sigma(2, 2, 0.0);
	Sigma.at(0, 0) = 0.1;
	Sigma.at(1, 1) = 0.2;
	double beta_1 = first_order_risk_estimation(mean, get_diag_vector_(Sigma));
	double beta_d = dth_order_risk_estimation(mean, get_diag_vector_(Sigma));

	// Test
	EXPECT_TRUE(beta_d <= beta_1);
	EXPECT_TRUE(beta_d >= 0);
	EXPECT_TRUE(beta_d <= 1);
}