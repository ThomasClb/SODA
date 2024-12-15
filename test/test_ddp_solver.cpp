/**
	test_ddp_solver.cpp

	Purpose: Test of the implementation of the DDPSolver class.

	@author Thomas Caleb

	@version 1.0 05/09/2024
*/

#include "gtest/gtest.h"
#include "../include/soda.h"

using namespace DACE;
using namespace std;


TEST(TestDDPSolver, Setters) {
	// Init DACE
	DA::init(2, 3);

	// Init
	vectordb x0(SIZE_VECTOR + 2, 1.0);
	vectordb xg(SIZE_VECTOR + 2, 1.0);
	vectorDA u(3, 1.0);
	Dynamics dynamics = get_tbp_SUN_lt_dynamics();
	SpacecraftParameters spacecraft_p(dynamics.constants());
	SolverParameters solver_p;
	DDPSolver ddp_solver(solver_p, spacecraft_p, dynamics);
	vector<vectordb> list_lambda(solver_p.N(), vectordb(3));
	vector<vectordb> list_mu(solver_p.N(), vectordb(3));
	double ToF = 10.0;
	double homotopy_coefficient = 0.5;
	double huber_loss_coefficient = 1e-2;
	double path_quantile = 5;
	double terminal_quantile = 8;
	matrixdb navigation_error_covariance(x0.size(), x0.size(), 0.0);

	// Set
	ddp_solver.set_list_lambda(list_lambda);
	ddp_solver.set_list_mu(list_mu);
	ddp_solver.set_homotopy_coefficient(homotopy_coefficient);
	ddp_solver.set_ToF(ToF);
	ddp_solver.set_huber_loss_coefficient(huber_loss_coefficient);
	ddp_solver.set_path_quantile(path_quantile);
	ddp_solver.set_terminal_quantile(terminal_quantile);
	ddp_solver.set_navigation_error_covariance(navigation_error_covariance);

	// Tests
	EXPECT_EQ(ddp_solver.solver_parameters().homotopy_coefficient(), homotopy_coefficient);
	EXPECT_EQ(ddp_solver.solver_parameters().ToF(), ToF);
	EXPECT_EQ(ddp_solver.solver_parameters().huber_loss_coefficient(), huber_loss_coefficient);
	EXPECT_EQ(ddp_solver.solver_parameters().path_quantile(), path_quantile);
	EXPECT_EQ(ddp_solver.solver_parameters().terminal_quantile(), terminal_quantile);
	for (size_t i = 0; i < solver_p.N(); i++) {
		vectordb lambda_i = ddp_solver.solver_parameters().list_lambda()[i];
		vectordb lambda_target_i = list_lambda[i];
		vectordb mu_i = ddp_solver.solver_parameters().list_mu()[i];
		vectordb mu_target_i = list_mu[i];
		for (size_t j = 0; j < mu_target_i.size(); j++) {
			EXPECT_EQ(lambda_i[j], lambda_target_i[j]);
			EXPECT_EQ(mu_i[j], mu_target_i[j]);
		}
	}
	for (size_t i=0; i<navigation_error_covariance.nrows(); i++) {
		for (size_t j=0; j<navigation_error_covariance.nrows(); j++) {
			EXPECT_EQ(
				ddp_solver.solver_parameters().navigation_error_covariance().at(i, j),
				navigation_error_covariance.at(i, j));
		}
	}
}
