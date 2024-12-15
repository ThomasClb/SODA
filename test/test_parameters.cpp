/**
	test_parameters.cpp

	Purpose: Test of the implementation of the SpacecraftParameters
	and SolverParameters classes.

	@author Thomas Caleb

	@version 2.0 12/12/2024
*/

#include "gtest/gtest.h"
#include "../include/soda.h"

using namespace DACE;
using namespace std;

/*

	SPACECRAFT PARAMETERS

*/

// Constructors
TEST(TestSpacecraftParameters, EmptyConstructor) {
	// Init
	Constants constants;
	SpacecraftParameters spacecraft_parameters;

	// Default values
	double initial_mass = 1000 / constants.massu(); // [MASSU]
	double dry_mass = 500 / constants.massu(); // [MASSU]
	double thrust = 0.5 / constants.thrustu(); // [THRUSTU]
	double Isp = 2000 / constants.tu(); // [TU]

	// Additionnal values
	double initial_wet_mass = initial_mass - dry_mass; // [MASSU]
	double ejection_velocity = G_0/(1000 * constants.vu()/ constants.tu()) * Isp; // [VU]
	double mass_flow = thrust / ejection_velocity; // [MASSU/TU]

	// Tests
	EXPECT_EQ(spacecraft_parameters.initial_mass(), initial_mass);
	EXPECT_EQ(spacecraft_parameters.thrust(), thrust);
	EXPECT_EQ(spacecraft_parameters.Isp(), Isp);
	EXPECT_EQ(spacecraft_parameters.dry_mass(), dry_mass);
	EXPECT_EQ(spacecraft_parameters.initial_wet_mass(), initial_wet_mass);
	EXPECT_EQ(spacecraft_parameters.ejection_velocity(), ejection_velocity);
	EXPECT_EQ(spacecraft_parameters.mass_flow(), mass_flow);
}
TEST(TestSpacecraftParameters, FilledConstructors) {
	// Init
	Constants constants;
	double initial_mass = 1000 / constants.massu(); // [MASSU]
	double dry_mass = 500 / constants.massu(); // [MASSU]
	double thrust = 0.5 / constants.thrustu(); // [THRUSTU]
	double Isp = 2000 / constants.tu(); // [TU]
	SpacecraftParameters spacecraft_parameters(
		constants,
		initial_mass, dry_mass, thrust, Isp);

	// Additionnal values
	double initial_wet_mass = initial_mass - dry_mass; // [MASSU]
	double ejection_velocity = G_0 / (1000 * constants.vu() / constants.tu()) * Isp; // [VU]
	double mass_flow = thrust / ejection_velocity; // [MASSU/TU]

	// Tests
	EXPECT_EQ(spacecraft_parameters.initial_mass(), initial_mass);
	EXPECT_EQ(spacecraft_parameters.thrust(), thrust);
	EXPECT_EQ(spacecraft_parameters.Isp(), Isp);
	EXPECT_EQ(spacecraft_parameters.dry_mass(), dry_mass);
	EXPECT_EQ(spacecraft_parameters.initial_wet_mass(), initial_wet_mass);
	EXPECT_EQ(spacecraft_parameters.ejection_velocity(), ejection_velocity);
	EXPECT_EQ(spacecraft_parameters.mass_flow(), mass_flow);

	// Second constructor
	spacecraft_parameters = SpacecraftParameters(
		constants);

	// Tests
	EXPECT_EQ(spacecraft_parameters.initial_mass(), initial_mass);
	EXPECT_EQ(spacecraft_parameters.thrust(), thrust);
	EXPECT_EQ(spacecraft_parameters.Isp(), Isp);
	EXPECT_EQ(spacecraft_parameters.dry_mass(), dry_mass);
	EXPECT_EQ(spacecraft_parameters.initial_wet_mass(), initial_wet_mass);
	EXPECT_EQ(spacecraft_parameters.ejection_velocity(), ejection_velocity);
	EXPECT_EQ(spacecraft_parameters.mass_flow(), mass_flow);
}
TEST(TestSpacecraftParameters, CopyConstructor) {
	// Init
	Constants constants;
	double initial_mass = 1000 / constants.massu(); // [MASSU]
	double dry_mass = 500 / constants.massu(); // [MASSU]
	double thrust = 0.5 / constants.thrustu(); // [THRUSTU]
	double Isp = 2000 / constants.tu(); // [TU]
	SpacecraftParameters spacecraft_parameters(
		constants,
		initial_mass, dry_mass, thrust, Isp);
	SpacecraftParameters spacecraft_parameters_copy = spacecraft_parameters;

	// Additionnal values
	double initial_wet_mass = initial_mass - dry_mass; // [MASSU]
	double ejection_velocity = G_0 / (1000 * constants.vu() / constants.tu()) * Isp; // [VU]
	double mass_flow = thrust / ejection_velocity; // [MASSU/TU]

	// Tests
	EXPECT_EQ(spacecraft_parameters_copy.initial_mass(), initial_mass);
	EXPECT_EQ(spacecraft_parameters_copy.thrust(), thrust);
	EXPECT_EQ(spacecraft_parameters_copy.Isp(), Isp);
	EXPECT_EQ(spacecraft_parameters_copy.dry_mass(), dry_mass);
	EXPECT_EQ(spacecraft_parameters_copy.initial_wet_mass(), initial_wet_mass);
	EXPECT_EQ(spacecraft_parameters_copy.ejection_velocity(), ejection_velocity);
	EXPECT_EQ(spacecraft_parameters_copy.mass_flow(), mass_flow);
}
TEST(TestSpacecraftParameters, IOFunctions) {
	// Init
	double mu = MU_MOON / (MU_EARTH + MU_MOON); // [-]
	double lu = EARTH_MOON_DISTANCE; // [km]
	double wu = sqrt((MU_EARTH + MU_MOON) / pow(EARTH_MOON_DISTANCE, 3)); // [s^-1]
	double tu = 1 / wu; // [s]
	double vu = lu * wu; // [m.s^-1]
	double massu = 1000; // [kg]
	double thrustu = 1000 * vu * massu * wu; // [N]
	Constants constants(mu, lu, wu, massu);
	double initial_mass = 1000 / constants.massu(); // [MASSU]
	double dry_mass = 500 / constants.massu(); // [MASSU]
	double thrust = 0.5 / constants.thrustu(); // [THRUSTU]
	double Isp = 2000 / constants.tu(); // [TU]
	double initial_wet_mass = initial_mass - dry_mass; // [MASSU]
	double ejection_velocity = G_0 / (1000 * constants.vu() / constants.tu()) * Isp; // [VU]
	double mass_flow = thrust / ejection_velocity; // [MASSU/TU]
	SpacecraftParameters spacecraft_parameters(
		constants,
		initial_mass, dry_mass, thrust, Isp);
	SpacecraftParameters spacecraft_parameters_copy;

	// Open file
	string file_name_("../data/spacecraft_parameters/test_spacecraft_parameters.dat");
	ofstream ofs(file_name_);

	// Store the object to file
	ofs << spacecraft_parameters;

	ofs.close();

	// Open file
	ifstream ifs(file_name_);

	// Load data
	ifs >> spacecraft_parameters_copy;

	ifs.close();

	// Tests
	EXPECT_EQ(spacecraft_parameters_copy.initial_mass(), initial_mass);
	EXPECT_EQ(spacecraft_parameters_copy.thrust(), thrust);
	EXPECT_EQ(spacecraft_parameters_copy.Isp(), Isp);
	EXPECT_EQ(spacecraft_parameters_copy.dry_mass(), dry_mass);
	EXPECT_EQ(spacecraft_parameters_copy.initial_wet_mass(), initial_wet_mass);
	EXPECT_EQ(spacecraft_parameters_copy.ejection_velocity(), ejection_velocity);
	EXPECT_EQ(spacecraft_parameters_copy.mass_flow(), mass_flow);
}
TEST(TestSpacecraftParameters, LoaderConstructor) {
	// Init
	double mu = MU_MOON / (MU_EARTH + MU_MOON); // [-]
	double lu = EARTH_MOON_DISTANCE; // [km]
	double wu = sqrt((MU_EARTH + MU_MOON) / pow(EARTH_MOON_DISTANCE, 3)); // [s^-1]
	double tu = 1 / wu; // [s]
	double vu = lu * wu; // [m.s^-1]
	double massu = 1000; // [kg]
	double thrustu = 1000 * vu * massu * wu; // [N]
	Constants constants(mu, lu, wu, massu);
	double initial_mass = 1000 / constants.massu(); // [MASSU]
	double dry_mass = 500 / constants.massu(); // [MASSU]
	double thrust = 0.5 / constants.thrustu(); // [THRUSTU]
	double Isp = 2000 / constants.tu(); // [TU]
	double initial_wet_mass = initial_mass - dry_mass; // [MASSU]
	double ejection_velocity = G_0 / (1000 * constants.vu() / constants.tu()) * Isp; // [VU]
	double mass_flow = thrust / ejection_velocity; // [MASSU/TU]
	SpacecraftParameters spacecraft_parameters(
		constants,
		initial_mass, dry_mass, thrust, Isp);
	string file_name_("../data/spacecraft_parameters/test_spacecraft_parameters.dat");
	spacecraft_parameters.save(file_name_);
	SpacecraftParameters spacecraft_parameters_copy(file_name_);

	// Tests
	EXPECT_EQ(spacecraft_parameters_copy.initial_mass(), initial_mass);
	EXPECT_EQ(spacecraft_parameters_copy.thrust(), thrust);
	EXPECT_EQ(spacecraft_parameters_copy.Isp(), Isp);
	EXPECT_EQ(spacecraft_parameters_copy.dry_mass(), dry_mass);
	EXPECT_EQ(spacecraft_parameters_copy.initial_wet_mass(), initial_wet_mass);
	EXPECT_EQ(spacecraft_parameters_copy.ejection_velocity(), ejection_velocity);
	EXPECT_EQ(spacecraft_parameters_copy.mass_flow(), mass_flow);
}

/*

	SOLVER PARAMETERS

*/

// Constructors
TEST(TestSolverParameters, EmptyConstructor) {
	// Init
	SolverParameters solver_parameters;

	// Default values
	unsigned int N = 40;
	unsigned int Nx = (SIZE_VECTOR + 1) + 1;
	unsigned int Nu = SIZE_VECTOR / 2;
	unsigned int Nineq = 2;
	unsigned int Ntineq = 0;
	bool with_J2 = false;
	double stage_cost_gain = 1e-2;
	double terminal_cost_gain = 1e4;
	matrixdb terminal_cost_inv_covariance = inv(make_diag_matrix_(sqr(vectordb(Nx, 1e-4))));
	matrixdb navigation_error_covariance = make_diag_matrix_(sqr(vectordb(Nx, 1e-4)));
	double transcription_beta(5e-2);
	double terminal_quantile = sqrt(inv_chi_2_cdf(Ntineq + 1, 1 - transcription_beta)); 
	double path_quantile = sqrt(inv_chi_2_cdf(Nineq + 1, 1 - transcription_beta));
	double mass_leak = 1e-8;
	double homotopy_coefficient = 0.0;
	double huber_loss_coefficient = 5e-3;
	vectordb homotopy_coefficient_sequence(1, 0);
	vectordb huber_loss_coefficient_sequence(1, 1e-2);
	unsigned int DDP_type = 0;
	double DDP_tol = 1e-4;
	double AUL_tol = 1e-4;
	double PN_tol = 1e-12;
	double PN_active_constraint_tol = 1e-13;
	unsigned int max_iter = 10000;
	unsigned int DDP_max_iter = 100;
	unsigned int AUL_max_iter = max_iter / DDP_max_iter;
	unsigned int PN_max_iter = 60;
	vectordb lambda_parameters{ 0.0, 1e8 };
	vectordb mu_parameters{ 1, 1e8, 10 };
	vectordb line_search_parameters{ 1e-8, 10.0, 0.5, 20 };
	bool backward_sweep_regulation = true;
	vectordb backward_sweep_regulation_parameters{ 0, 1e-8, 1e8, 1.6 };
	double PN_regularisation(1e-8);
	double PN_cv_rate_threshold(1.1);
	double PN_alpha(1.0); double PN_gamma(0.5);
	double ToF = 0.0;
	unsigned int verbosity = 0;
	unsigned int saving_iterations = 0;

	// Tests
	EXPECT_EQ(solver_parameters.N(), N);
	EXPECT_EQ(solver_parameters.Nx(), Nx);
	EXPECT_EQ(solver_parameters.Nu(), Nu);
	EXPECT_EQ(solver_parameters.Nineq(), Nineq);
	EXPECT_EQ(solver_parameters.Ntineq(), Ntineq);
	EXPECT_EQ(solver_parameters.ToF(), ToF);
	EXPECT_EQ(solver_parameters.with_J2(), with_J2);
	EXPECT_EQ(solver_parameters.stage_cost_gain(), stage_cost_gain);
	EXPECT_EQ(solver_parameters.terminal_cost_gain(), terminal_cost_gain);
	EXPECT_EQ(solver_parameters.transcription_beta(), transcription_beta);
	EXPECT_EQ(solver_parameters.terminal_quantile(), terminal_quantile);
	EXPECT_EQ(solver_parameters.path_quantile(), path_quantile);
	for (size_t i=0; i<solver_parameters.navigation_error_covariance().nrows(); i++) {
		for (size_t j=0; j<solver_parameters.navigation_error_covariance().ncols(); j++) {
			EXPECT_EQ(solver_parameters.navigation_error_covariance().at(i, j), navigation_error_covariance.at(i, j));
			EXPECT_EQ(solver_parameters.terminal_cost_inv_covariance().at(i, j), terminal_cost_inv_covariance.at(i, j));
		}
	}
	EXPECT_EQ(solver_parameters.mass_leak(), mass_leak);
	EXPECT_EQ(solver_parameters.homotopy_coefficient(), homotopy_coefficient);
	EXPECT_EQ(solver_parameters.huber_loss_coefficient(), huber_loss_coefficient);
	EXPECT_EQ(solver_parameters.homotopy_coefficient_sequence(), homotopy_coefficient_sequence);
	EXPECT_EQ(solver_parameters.huber_loss_coefficient_sequence(), huber_loss_coefficient_sequence);
	EXPECT_EQ(solver_parameters.DDP_type(), DDP_type);
	EXPECT_EQ(solver_parameters.DDP_tol(), DDP_tol);
	EXPECT_EQ(solver_parameters.AUL_tol(), AUL_tol);
	EXPECT_EQ(solver_parameters.PN_tol(), PN_tol);
	EXPECT_EQ(solver_parameters.DDP_max_iter(), DDP_max_iter);
	EXPECT_EQ(solver_parameters.AUL_max_iter(), AUL_max_iter);
	EXPECT_EQ(solver_parameters.PN_max_iter(), PN_max_iter);
	EXPECT_EQ(solver_parameters.line_search_parameters(), line_search_parameters);
	EXPECT_EQ(solver_parameters.backward_sweep_regulation(), backward_sweep_regulation);
	EXPECT_EQ(solver_parameters.backward_sweep_regulation_parameters(), backward_sweep_regulation_parameters);
	EXPECT_EQ(solver_parameters.lambda_parameters(), lambda_parameters);
	EXPECT_EQ(solver_parameters.mu_parameters(), mu_parameters);
	EXPECT_EQ(solver_parameters.PN_regularisation(), PN_regularisation);
	EXPECT_EQ(solver_parameters.PN_active_constraint_tol(), PN_active_constraint_tol);
	EXPECT_EQ(solver_parameters.PN_cv_rate_threshold(), PN_cv_rate_threshold);
	EXPECT_EQ(solver_parameters.PN_alpha(), PN_alpha);
	EXPECT_EQ(solver_parameters.PN_gamma(), PN_gamma);
	EXPECT_EQ(solver_parameters.verbosity(), verbosity);
	EXPECT_EQ(solver_parameters.saving_iterations(), saving_iterations);
	for (size_t i = 0; i < N; i++) {
		vectordb lambda_i = solver_parameters.list_lambda()[i];
		vectordb mu_i = solver_parameters.list_mu()[i];
		for (size_t j = 0; j < Nineq; j++) {
			EXPECT_EQ(lambda_i[j], lambda_parameters[0]);
			EXPECT_EQ(mu_i[j], mu_parameters[0]);
		}
	}
	vectordb lambda_N = solver_parameters.list_lambda()[N];
	vectordb mu_N = solver_parameters.list_mu()[N];
	for (size_t j = 0; j < Ntineq; j++) {
		EXPECT_EQ(lambda_N[j], lambda_parameters[0]);
		EXPECT_EQ(mu_N[j], mu_parameters[0]);
	}
}
TEST(TestSolverParameters, FilledConstructor) {
	// Init
	unsigned int N = 40;
	unsigned int Nx = (SIZE_VECTOR + 1) + 1;
	unsigned int Nu = SIZE_VECTOR / 2;
	unsigned int Nineq = 2;
	unsigned int Ntineq = 0;
	bool with_J2 = true;
	double stage_cost_gain = 1e-2;
	double terminal_cost_gain = 1e4;
	matrixdb terminal_cost_inv_covariance = inv(make_diag_matrix_(sqr(vectordb(Nx, 2e-4))));
	matrixdb navigation_error_covariance = make_diag_matrix_(sqr(vectordb(Nx, 2e-4)));
	double transcription_beta(1e-2);
	double terminal_quantile = sqrt(inv_chi_2_cdf(Ntineq + 1, 1 - transcription_beta)); 
	double path_quantile = sqrt(inv_chi_2_cdf(Nineq + 1, 1 - transcription_beta));
	double mass_leak = 1e-10;
	double homotopy_coefficient = 0.0;
	double huber_loss_coefficient = 5e-3;
	vectordb homotopy_coefficient_sequence(2, 0);
	vectordb huber_loss_coefficient_sequence(2, 1e-2);
	unsigned int DDP_type = 0;
	double DDP_tol = 1e-4;
	double AUL_tol = 1e-4;
	double PN_tol = 1e-12;
	double PN_active_constraint_tol = 1e-13;
	unsigned int max_iter = 10000;
	unsigned int DDP_max_iter = 100;
	unsigned int AUL_max_iter = max_iter / DDP_max_iter;
	unsigned int PN_max_iter = 60;
	vectordb lambda_parameters{ 0.0, 1e8 };
	vectordb mu_parameters{ 1, 1e8, 10 };
	vectordb line_search_parameters{ 1e-8, 10.0, 0.5, 20 };
	bool backward_sweep_regulation = true;
	vectordb backward_sweep_regulation_parameters{ 0, 1e-8, 1e8, 1.6 };
	double PN_regularisation(1e-8);
	double PN_cv_rate_threshold(1.1);
	double PN_alpha(1.0); double PN_gamma(0.5);
	double ToF = 0.0;
	unsigned int verbosity = 1;
	unsigned int saving_iterations = 1;
	SolverParameters solver_parameters(
		N, Nx, Nu,
		Nineq, Ntineq, with_J2,
		stage_cost_gain, terminal_cost_gain,
		terminal_cost_inv_covariance,
		navigation_error_covariance,
		transcription_beta,mass_leak,
		homotopy_coefficient, huber_loss_coefficient,
		homotopy_coefficient_sequence,
		huber_loss_coefficient_sequence,
		DDP_type,
		DDP_tol, AUL_tol, PN_tol,
		DDP_max_iter, AUL_max_iter, PN_max_iter,
		line_search_parameters,
		backward_sweep_regulation,
		backward_sweep_regulation_parameters,
		lambda_parameters, mu_parameters,
		PN_regularisation, PN_active_constraint_tol,
		PN_cv_rate_threshold, PN_alpha, PN_gamma, verbosity,
		saving_iterations);

	// Tests
	EXPECT_EQ(solver_parameters.N(), N);
	EXPECT_EQ(solver_parameters.Nx(), Nx);
	EXPECT_EQ(solver_parameters.Nu(), Nu);
	EXPECT_EQ(solver_parameters.Nineq(), Nineq);
	EXPECT_EQ(solver_parameters.Ntineq(), Ntineq);
	EXPECT_EQ(solver_parameters.ToF(), ToF);
	EXPECT_EQ(solver_parameters.with_J2(), with_J2);
	EXPECT_EQ(solver_parameters.stage_cost_gain(), stage_cost_gain);
	EXPECT_EQ(solver_parameters.terminal_cost_gain(), terminal_cost_gain);
	EXPECT_EQ(solver_parameters.transcription_beta(), transcription_beta);
	EXPECT_EQ(solver_parameters.terminal_quantile(), terminal_quantile);
	EXPECT_EQ(solver_parameters.path_quantile(), path_quantile);
	for (size_t i=0; i<solver_parameters.navigation_error_covariance().nrows(); i++) {
		for (size_t j=0; j<solver_parameters.navigation_error_covariance().ncols(); j++) {
			EXPECT_EQ(solver_parameters.navigation_error_covariance().at(i, j), navigation_error_covariance.at(i, j));
			EXPECT_EQ(solver_parameters.terminal_cost_inv_covariance().at(i, j), terminal_cost_inv_covariance.at(i, j));
		}
	}
	EXPECT_EQ(solver_parameters.mass_leak(), mass_leak);
	EXPECT_EQ(solver_parameters.homotopy_coefficient(), homotopy_coefficient);
	EXPECT_EQ(solver_parameters.huber_loss_coefficient(), huber_loss_coefficient);
	EXPECT_EQ(solver_parameters.homotopy_coefficient_sequence(), homotopy_coefficient_sequence);
	EXPECT_EQ(solver_parameters.huber_loss_coefficient_sequence(), huber_loss_coefficient_sequence);
	EXPECT_EQ(solver_parameters.DDP_type(), DDP_type);
	EXPECT_EQ(solver_parameters.DDP_tol(), DDP_tol);
	EXPECT_EQ(solver_parameters.AUL_tol(), AUL_tol);
	EXPECT_EQ(solver_parameters.PN_tol(), PN_tol);
	EXPECT_EQ(solver_parameters.DDP_max_iter(), DDP_max_iter);
	EXPECT_EQ(solver_parameters.AUL_max_iter(), AUL_max_iter);
	EXPECT_EQ(solver_parameters.PN_max_iter(), PN_max_iter);
	EXPECT_EQ(solver_parameters.line_search_parameters(), line_search_parameters);
	EXPECT_EQ(solver_parameters.backward_sweep_regulation(), backward_sweep_regulation);
	EXPECT_EQ(solver_parameters.backward_sweep_regulation_parameters(), backward_sweep_regulation_parameters);
	EXPECT_EQ(solver_parameters.lambda_parameters(), lambda_parameters);
	EXPECT_EQ(solver_parameters.mu_parameters(), mu_parameters);
	EXPECT_EQ(solver_parameters.PN_regularisation(), PN_regularisation);
	EXPECT_EQ(solver_parameters.PN_active_constraint_tol(), PN_active_constraint_tol);
	EXPECT_EQ(solver_parameters.PN_cv_rate_threshold(), PN_cv_rate_threshold);
	EXPECT_EQ(solver_parameters.PN_alpha(), PN_alpha);
	EXPECT_EQ(solver_parameters.PN_gamma(), PN_gamma);
	EXPECT_EQ(solver_parameters.verbosity(), verbosity);
	EXPECT_EQ(solver_parameters.saving_iterations(), saving_iterations);
	for (size_t i = 0; i < N; i++) {
		vectordb lambda_i = solver_parameters.list_lambda()[i];
		vectordb mu_i = solver_parameters.list_mu()[i];
		for (size_t j = 0; j < Nineq; j++) {
			EXPECT_EQ(lambda_i[j], lambda_parameters[0]);
			EXPECT_EQ(mu_i[j], mu_parameters[0]);
		}
	}
	vectordb lambda_N = solver_parameters.list_lambda()[N];
	vectordb mu_N = solver_parameters.list_mu()[N];
	for (size_t j = 0; j < Ntineq; j++) {
		EXPECT_EQ(lambda_N[j], lambda_parameters[0]);
		EXPECT_EQ(mu_N[j], mu_parameters[0]);
	}
}
TEST(TestSolverParameters, CopyConstructor) {
	// Init
	unsigned int N = 40;
	unsigned int Nx = (SIZE_VECTOR + 1) + 1;
	unsigned int Nu = SIZE_VECTOR / 2;
	unsigned int Nineq = 2;
	unsigned int Ntineq = 0;
	bool with_J2 = true;
	double stage_cost_gain = 1e-2;
	double terminal_cost_gain = 1e4;
	matrixdb terminal_cost_inv_covariance = inv(make_diag_matrix_(sqr(vectordb(Nx, 2e-4))));
	matrixdb navigation_error_covariance = make_diag_matrix_(sqr(vectordb(Nx, 2e-4)));
	double transcription_beta(1e-2);
	double terminal_quantile = sqrt(inv_chi_2_cdf(Ntineq + 1, 1 - transcription_beta)); 
	double path_quantile = sqrt(inv_chi_2_cdf(Nineq + 1, 1 - transcription_beta));
	double mass_leak = 1e-10;
	double homotopy_coefficient = 0.0;
	double huber_loss_coefficient = 5e-3;
	vectordb homotopy_coefficient_sequence(2, 0);
	vectordb huber_loss_coefficient_sequence(2, 1e-2);
	unsigned int DDP_type = 0;
	double DDP_tol = 1e-4;
	double AUL_tol = 1e-4;
	double PN_tol = 1e-12;
	double PN_active_constraint_tol = 1e-13;
	unsigned int max_iter = 10000;
	unsigned int DDP_max_iter = 100;
	unsigned int AUL_max_iter = max_iter / DDP_max_iter;
	unsigned int PN_max_iter = 60;
	vectordb lambda_parameters{ 0.0, 1e8 };
	vectordb mu_parameters{ 1, 1e8, 10 };
	vectordb line_search_parameters{ 1e-8, 10.0, 0.5, 20 };
	bool backward_sweep_regulation = true;
	vectordb backward_sweep_regulation_parameters{ 0, 1e-8, 1e8, 1.6 };
	double PN_regularisation(1e-8);
	double PN_cv_rate_threshold(1.1);
	double PN_alpha(1.0); double PN_gamma(0.5);
	double ToF = 0.0;
	unsigned int verbosity = 1;
	unsigned int saving_iterations = 1;
	SolverParameters solver_parameters(
		N, Nx, Nu,
		Nineq, Ntineq, with_J2,
		stage_cost_gain, terminal_cost_gain,
		terminal_cost_inv_covariance,
		navigation_error_covariance,
		transcription_beta,mass_leak,
		homotopy_coefficient, huber_loss_coefficient,
		homotopy_coefficient_sequence,
		huber_loss_coefficient_sequence,
		DDP_type,
		DDP_tol, AUL_tol, PN_tol,
		DDP_max_iter, AUL_max_iter, PN_max_iter,
		line_search_parameters,
		backward_sweep_regulation,
		backward_sweep_regulation_parameters,
		lambda_parameters, mu_parameters,
		PN_regularisation, PN_active_constraint_tol,
		PN_cv_rate_threshold, PN_alpha, PN_gamma, verbosity,
		saving_iterations);
	SolverParameters solver_parameters_copy = solver_parameters;

	// Tests
	EXPECT_EQ(solver_parameters_copy.N(), N);
	EXPECT_EQ(solver_parameters_copy.Nx(), Nx);
	EXPECT_EQ(solver_parameters_copy.Nu(), Nu);
	EXPECT_EQ(solver_parameters_copy.Nineq(), Nineq);
	EXPECT_EQ(solver_parameters_copy.Ntineq(), Ntineq);
	EXPECT_EQ(solver_parameters_copy.with_J2(), with_J2);
	EXPECT_EQ(solver_parameters_copy.ToF(), ToF);
	EXPECT_EQ(solver_parameters_copy.stage_cost_gain(), stage_cost_gain);
	EXPECT_EQ(solver_parameters_copy.terminal_cost_gain(), terminal_cost_gain);
	EXPECT_EQ(solver_parameters_copy.transcription_beta(), transcription_beta);
	EXPECT_EQ(solver_parameters_copy.terminal_quantile(), terminal_quantile);
	EXPECT_EQ(solver_parameters_copy.path_quantile(), path_quantile);
	for (size_t i=0; i<solver_parameters.navigation_error_covariance().nrows(); i++) {
		for (size_t j=0; j<solver_parameters.navigation_error_covariance().ncols(); j++) {
			EXPECT_EQ(solver_parameters_copy.navigation_error_covariance().at(i, j), navigation_error_covariance.at(i, j));
			EXPECT_EQ(solver_parameters_copy.terminal_cost_inv_covariance().at(i, j), terminal_cost_inv_covariance.at(i, j));
		}
	}
	EXPECT_EQ(solver_parameters_copy.mass_leak(), mass_leak);
	EXPECT_EQ(solver_parameters_copy.homotopy_coefficient(), homotopy_coefficient);
	EXPECT_EQ(solver_parameters_copy.huber_loss_coefficient(), huber_loss_coefficient);
	EXPECT_EQ(solver_parameters_copy.homotopy_coefficient_sequence(), homotopy_coefficient_sequence);
	EXPECT_EQ(solver_parameters_copy.huber_loss_coefficient_sequence(), huber_loss_coefficient_sequence);
	EXPECT_EQ(solver_parameters_copy.DDP_type(), DDP_type);
	EXPECT_EQ(solver_parameters_copy.DDP_tol(), DDP_tol);
	EXPECT_EQ(solver_parameters_copy.AUL_tol(), AUL_tol);
	EXPECT_EQ(solver_parameters_copy.PN_tol(), PN_tol);
	EXPECT_EQ(solver_parameters_copy.DDP_max_iter(), DDP_max_iter);
	EXPECT_EQ(solver_parameters_copy.AUL_max_iter(), AUL_max_iter);
	EXPECT_EQ(solver_parameters_copy.PN_max_iter(), PN_max_iter);
	EXPECT_EQ(solver_parameters_copy.line_search_parameters(), line_search_parameters);
	EXPECT_EQ(solver_parameters_copy.backward_sweep_regulation(), backward_sweep_regulation);
	EXPECT_EQ(solver_parameters_copy.backward_sweep_regulation_parameters(), backward_sweep_regulation_parameters);
	EXPECT_EQ(solver_parameters_copy.lambda_parameters(), lambda_parameters);
	EXPECT_EQ(solver_parameters_copy.mu_parameters(), mu_parameters);
	EXPECT_EQ(solver_parameters_copy.PN_regularisation(), PN_regularisation);
	EXPECT_EQ(solver_parameters_copy.PN_active_constraint_tol(), PN_active_constraint_tol);
	EXPECT_EQ(solver_parameters_copy.PN_cv_rate_threshold(), PN_cv_rate_threshold);
	EXPECT_EQ(solver_parameters_copy.PN_alpha(), PN_alpha);
	EXPECT_EQ(solver_parameters_copy.PN_gamma(), PN_gamma);
	EXPECT_EQ(solver_parameters_copy.verbosity(), verbosity);
	EXPECT_EQ(solver_parameters_copy.saving_iterations(), saving_iterations);
	for (size_t i = 0; i < N; i++) {
		vectordb lambda_i = solver_parameters_copy.list_lambda()[i];
		vectordb mu_i = solver_parameters_copy.list_mu()[i];
		for (size_t j = 0; j < Nineq; j++) {
			EXPECT_EQ(lambda_i[j], lambda_parameters[0]);
			EXPECT_EQ(mu_i[j], mu_parameters[0]);
		}
	}
	vectordb lambda_N = solver_parameters_copy.list_lambda()[N];
	vectordb mu_N = solver_parameters_copy.list_mu()[N];
	for (size_t j = 0; j < Ntineq; j++) {
		EXPECT_EQ(lambda_N[j], lambda_parameters[0]);
		EXPECT_EQ(mu_N[j], mu_parameters[0]);
	}
}

// Setters
TEST(TestSolverParameters, Setters) {
	// Init
	unsigned int N = 40;
	unsigned int Nx = (SIZE_VECTOR + 1) + 1;
	unsigned int Nu = SIZE_VECTOR / 2;
	unsigned int Nineq = 2;
	unsigned int Ntineq = 0;
	double ToF = 10.0;
	double homotopy_coefficient = 0.0;
	double huber_loss_coefficient = 10.0;
	vectordb lambda_parameters{ 0.0, 1e8 };
	vectordb mu_parameters{ 1, 1e8, 10 };
	double transcription_beta = 0.05;
	double terminal_quantile = 3.0; 
	double path_quantile = 5.0;
	matrixdb navigation_error_covariance = make_diag_matrix_(sqr(vectordb(Nx, 2e-4)));
	SolverParameters solver_parameters;

	// Changes
	vector<vectordb> list_lambda(N, vectordb(N, 0.0)), list_mu(N, vectordb(N, 0.0));
	homotopy_coefficient = 42.0;

	// Set
	solver_parameters.set_list_lambda(list_lambda);
	solver_parameters.set_list_mu(list_mu);
	solver_parameters.set_homotopy_coefficient(homotopy_coefficient);
	solver_parameters.set_ToF(ToF);
	solver_parameters.set_huber_loss_coefficient(huber_loss_coefficient);
	solver_parameters.set_terminal_quantile(terminal_quantile);
	solver_parameters.set_path_quantile(path_quantile);
	solver_parameters.set_navigation_error_covariance(navigation_error_covariance);

	// Tests
	EXPECT_EQ(solver_parameters.homotopy_coefficient(), homotopy_coefficient);
	EXPECT_EQ(solver_parameters.ToF(), ToF);
	EXPECT_EQ(solver_parameters.huber_loss_coefficient(), huber_loss_coefficient);
	EXPECT_EQ(solver_parameters.terminal_quantile(), terminal_quantile);
	EXPECT_EQ(solver_parameters.path_quantile(), path_quantile);
	for (size_t i=0; i<solver_parameters.navigation_error_covariance().nrows(); i++) {
		for (size_t j=0; j<solver_parameters.navigation_error_covariance().ncols(); j++) {
			EXPECT_EQ(solver_parameters.navigation_error_covariance().at(i, j), navigation_error_covariance.at(i, j));
		}
	}
	for (size_t i = 0; i < N; i++) {
		vectordb lambda_i = solver_parameters.list_lambda()[i];
		vectordb lambda_target_i = list_lambda[i];
		vectordb mu_i = solver_parameters.list_mu()[i];
		vectordb mu_target_i = list_mu[i];
		for (size_t j = 0; j < N; j++) {
			EXPECT_EQ(lambda_i[j], lambda_target_i[j]);
			EXPECT_EQ(mu_i[j], mu_target_i[j]);
		}
	}
}