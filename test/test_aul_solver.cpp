/**
	test_aul_solver.cpp

	Purpose: Test of the implementation of the AULSolver class.

	@author Thomas Caleb

	@version 2.0 14/12/2024
*/

#include "gtest/gtest.h"
#include "../include/soda.h"

using namespace DACE;
using namespace std;

TEST(TestAULSolver, Setters) {
	// Init DACE
	DA::init(2, 3);

	// Init
	vectordb x0(SIZE_VECTOR + 2, 1.0);
	vectordb xg(SIZE_VECTOR + 2, 1.0);
	vectorDA u(3, 1.0);
	Dynamics dynamics = get_tbp_SUN_lt_dynamics();
	SpacecraftParameters spacecraft_p(dynamics.constants());
	SolverParameters solver_p;
	AULSolver aul_solver(solver_p, spacecraft_p, dynamics);
	vector<vectordb> list_lambda(solver_p.N(), vectordb(3));
	vector<vectordb> list_mu(solver_p.N(), vectordb(3));
	double ToF = 10.0;
	double homotopy_coefficient = 0.5;
	double huber_loss_coefficient = 1e-2;
	double path_quantile = 5;
	double terminal_quantile = 8;
	matrixdb navigation_error_covariance(x0.size(), x0.size(), 0.0);

	// Set
	aul_solver.set_homotopy_coefficient(homotopy_coefficient);
	aul_solver.set_ToF(ToF);
	aul_solver.set_huber_loss_coefficient(huber_loss_coefficient);
	aul_solver.set_path_quantile(path_quantile);
	aul_solver.set_terminal_quantile(terminal_quantile);
	aul_solver.set_navigation_error_covariance(navigation_error_covariance);
	DDPSolver ddp_solver = aul_solver.DDPsolver();
	SolverParameters solver_parameters = ddp_solver.solver_parameters();

	// Tests
	EXPECT_EQ(solver_parameters.homotopy_coefficient(), homotopy_coefficient);
	EXPECT_EQ(solver_parameters.ToF(), ToF);
	EXPECT_EQ(solver_parameters.huber_loss_coefficient(), huber_loss_coefficient);
	EXPECT_EQ(ddp_solver.solver_parameters().path_quantile(), path_quantile);
	EXPECT_EQ(ddp_solver.solver_parameters().terminal_quantile(), terminal_quantile);
	for (size_t i=0; i<navigation_error_covariance.nrows(); i++) {
		for (size_t j=0; j<navigation_error_covariance.nrows(); j++) {
			EXPECT_EQ(
				ddp_solver.solver_parameters().navigation_error_covariance().at(i, j),
				navigation_error_covariance.at(i, j));
		}
	}
}