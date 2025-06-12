/**
	test_transcription.cpp

	Purpose: Test of the implementation of the the transcription methods.
    
	@author Thomas Caleb

	@version 1.0 13/12/2024
*/

#include "gtest/gtest.h"
#include "../include/soda.h"

using namespace DACE;
using namespace std;



TEST(TestTranscription, FODDPPathConstraints) {
	// Init DACE
	SolverParameters solver_parameters;
	size_t nb_variables = solver_parameters.Nu() + solver_parameters.Nx();
	DA::init(3, nb_variables);

	// Init
	SpacecraftParameters spacecraft_parameters;
	Constants constants;
	vectorDA constraints_eval{-1.0 + DA(1)*DA(3) + DA(1) + DA(2) + DA(3), -1.5 + DA(1) + DA(2) + DA(3) + DA(2) + DA(3)*DA(2)};
	stateDA x_DA(vectorDA(solver_parameters.Nx(), 1.0 + DA(1) + DA(1)*DA(2)*DA(3)));
	controlDA u_DA(vectorDA(solver_parameters.Nu(), 0.01 + DA(3) + 0.01*DA(1)*DA(2)*DA(3)));
	matrixdb Sigma(solver_parameters.Nx(), solver_parameters.Nx(), 0.0);
	matrixdb der_dynamics(solver_parameters.Nx(), solver_parameters.Nx()+solver_parameters.Nu(), 1.0);
	x_DA.set_der_dynamics(der_dynamics);
	for (size_t i=0; i<solver_parameters.Nx(); i++) {
		Sigma.at(i, i) = 0.0001*(1.0+i);
	}
	x_DA.set_Sigma(Sigma);
	matrixdb feedback_gain(solver_parameters.Nu(), solver_parameters.Nx(), 1.0);
	u_DA.set_feedback_gain(feedback_gain);

	vectorDA constraints_eval_trans = first_order_path_transcription(
		constraints_eval, x_DA, u_DA,
		spacecraft_parameters, constants, solver_parameters);

	// Test
	EXPECT_EQ(constraints_eval_trans.size(), constraints_eval.size());
	for (size_t i=0; i<constraints_eval_trans.size(); i++) {
		EXPECT_TRUE(constraints_eval_trans[i].cons() > constraints_eval[i].cons());
	}
}
TEST(TestTranscription, SRDDPPathConstraints) {
	// Init DACE
	SolverParameters solver_parameters;
	size_t nb_variables = solver_parameters.Nu() + solver_parameters.Nx();
	DA::init(3, nb_variables);

	// Init
	SpacecraftParameters spacecraft_parameters;
	Constants constants;
	vectorDA constraints_eval{-1.0 + DA(1)*DA(3) + DA(1) + DA(2) + DA(3), -1.5 + DA(1) + DA(2) + DA(3) + DA(2) + DA(3)*DA(2)};
	stateDA x_DA(vectorDA(solver_parameters.Nx(), 1.0 + DA(1) + DA(1)*DA(2)*DA(3)));
	controlDA u_DA(vectorDA(solver_parameters.Nu(), 0.01 + DA(3) + 0.01*DA(1)*DA(2)*DA(3)));
	matrixdb Sigma(solver_parameters.Nx(), solver_parameters.Nx(), 0.0);
	matrixdb der_dynamics(solver_parameters.Nx(), solver_parameters.Nx()+solver_parameters.Nu(), 1.0);
	x_DA.set_der_dynamics(der_dynamics);
	for (size_t i=0; i<solver_parameters.Nx(); i++) {
		Sigma.at(i, i) = 0.0001*(1.0+i);
	}
	x_DA.set_Sigma(Sigma);
	matrixdb feedback_gain(solver_parameters.Nu(), solver_parameters.Nx(), 1.0);
	u_DA.set_feedback_gain(feedback_gain);

	vectorDA constraints_eval_trans = spectral_radius_path_transcription(
		constraints_eval, x_DA, u_DA,
		spacecraft_parameters, constants, solver_parameters);

	// Test
	EXPECT_EQ(constraints_eval_trans.size(), constraints_eval.size());
	for (size_t i=0; i<constraints_eval_trans.size(); i++) {
		EXPECT_TRUE(constraints_eval_trans[i].cons() > constraints_eval[i].cons());
	}
}
TEST(TestTranscription, FODDPTerminalConstraints) {
	// Init DACE
	SolverParameters solver_parameters;
	size_t nb_variables = solver_parameters.Nu() + solver_parameters.Nx();
	DA::init(3, nb_variables);

	// Init
	SpacecraftParameters spacecraft_parameters;
	Constants constants;
	vectorDA constraints_eval{-1.0 + DA(1)*DA(3) + DA(1) + DA(2) + DA(3), -1.5 + DA(1) + DA(2) + DA(3) + DA(2) + DA(3)*DA(2)};
	stateDA x_DA(vectorDA(solver_parameters.Nx(), 1.0 + DA(1) + DA(1)*DA(2)*DA(3)));
	matrixdb Sigma(solver_parameters.Nx(), solver_parameters.Nx(), 0.0);
	matrixdb der_dynamics(solver_parameters.Nx(), solver_parameters.Nx()+solver_parameters.Nu(), 1.0);
	x_DA.set_der_dynamics(der_dynamics);
	for (size_t i=0; i<solver_parameters.Nx(); i++) {
		Sigma.at(i, i) = 0.0001*(1.0+i);
	}
	x_DA.set_Sigma(Sigma);
	statedb x_goal(x_DA.cons());

	vectorDA constraints_eval_trans = first_order_terminal_transcription(
		constraints_eval, x_DA, x_goal,
		spacecraft_parameters, constants, solver_parameters);

	// Test
	EXPECT_EQ(constraints_eval_trans.size(), constraints_eval.size());
	for (size_t i=0; i<constraints_eval_trans.size(); i++) {
		EXPECT_TRUE(constraints_eval_trans[i].cons() > constraints_eval[i].cons());
	}
}
TEST(TestTranscription, SRDDPTerminalConstraints) {
	// Init DACE
	SolverParameters solver_parameters;
	size_t nb_variables = solver_parameters.Nu() + solver_parameters.Nx();
	DA::init(3, nb_variables);

	// Init
	SpacecraftParameters spacecraft_parameters;
	Constants constants;
	vectorDA constraints_eval{-1.0 + DA(1)*DA(3) + DA(1) + DA(2) + DA(3), -1.5 + DA(1) + DA(2) + DA(3) + DA(2) + DA(3)*DA(2)};
	stateDA x_DA(vectorDA(solver_parameters.Nx(), 1.0 + DA(1) + DA(1)*DA(2)*DA(3)));
	statedb x_goal(x_DA.cons());
	matrixdb Sigma(solver_parameters.Nx(), solver_parameters.Nx(), 0.0);
	matrixdb der_dynamics(solver_parameters.Nx(), solver_parameters.Nx()+solver_parameters.Nu(), 1.0);
	x_DA.set_der_dynamics(der_dynamics);
	for (size_t i=0; i<solver_parameters.Nx(); i++) {
		Sigma.at(i, i) = 0.0001*(1.0+i);
	}
	x_DA.set_Sigma(Sigma);

	vectorDA constraints_eval_trans = spectral_radius_terminal_transcription(
		constraints_eval, x_DA, x_goal,
		spacecraft_parameters, constants, solver_parameters);

	// Test
	EXPECT_EQ(constraints_eval_trans.size(), constraints_eval.size());
	for (size_t i=0; i<constraints_eval_trans.size(); i++) {
		EXPECT_TRUE(constraints_eval_trans[i].cons() > constraints_eval[i].cons());
	}
}

TEST(TestTranscription, FOPN) {
	// Init DACE
	SolverParameters solver_parameters;
	size_t nb_variables = solver_parameters.Nu() + solver_parameters.Nx();
	DA::init(3, nb_variables);

	// Init
	SpacecraftParameters spacecraft_parameters;
	Constants constants;
	solver_parameters.set_path_quantile(10);
	double eta = 0.75;
	size_t N = solver_parameters.N();
	size_t Nx = solver_parameters.Nx();
	size_t Nu = solver_parameters.Nu();
	size_t Nineq = solver_parameters.Nineq();
	size_t Ntineq = solver_parameters.Ntineq();
	vectorDA constraints_eval(N*Nineq + Ntineq);
	for (size_t i=0; i<constraints_eval.size(); i++) {
		if (i%3 == 0)
			constraints_eval[i] = -1.0*i + DA(1)*DA(3) + DA(1) + DA(5) + DA(7);
		else if (i%3 == 1)
			constraints_eval[i] = (0.33*i + DA(5) + DA(2) + DA(8) + DA(2) + DA(3)*DA(2))*-1.0 + DA(1)*DA(3) + DA(5) + DA(2) + DA(3);
		else if (i%3 == 2)
			constraints_eval[i] = -0.27*i + DA(1) + DA(2) + DA(3) + DA(10) + DA(3)*DA(3);
	}

	matrixdb Sigma(solver_parameters.Nx(), solver_parameters.Nx(), 0.0);
	matrixdb feedback_gain(solver_parameters.Nu(), solver_parameters.Nx(), 0.01);
	for (size_t i=0; i<solver_parameters.Nx(); i++) {
		Sigma.at(i, i) = 0.0001*(1.0+i);
	}
	vector<matrixdb> list_Sigma(N+1, Sigma);
	vector<matrixdb> list_feedback_gain(N, feedback_gain);
	vectorDA constraints_eval_trans = first_order_transcription(
		constraints_eval, list_Sigma, list_feedback_gain,
		eta,
		spacecraft_parameters, constants, solver_parameters);
 
	// Test
	EXPECT_EQ(constraints_eval_trans.size(), constraints_eval.size());
	for (size_t i=0; i<constraints_eval_trans.size(); i++) {
		EXPECT_TRUE(constraints_eval_trans[i].cons() > constraints_eval[i].cons());
	}
}