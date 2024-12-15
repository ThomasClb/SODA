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
	size_t nb_variables = 3;
	DA::init(3, nb_variables);

	// Init
	SolverParameters solver_parameters;
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

	vectorDA constraints_eval_trans = first_order_path_inequality_transcription(
		constraints_eval, x_DA, u_DA,
		spacecraft_parameters, constants, solver_parameters);

	// Test
	EXPECT_EQ(constraints_eval_trans.size(), constraints_eval.size());
	for (size_t i=0; i<constraints_eval_trans.size(); i++) {
		EXPECT_TRUE(constraints_eval_trans[i].cons() > constraints_eval[i].cons());
	}
}
TEST(TestTranscription, FOPnPathConstraints) {
	// Init DACE
	size_t nb_variables = 3;
	DA::init(3, nb_variables);

	// Init
	SolverParameters solver_parameters;
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

	vectorDA constraints_eval_trans = first_order_path_inequality_transcription( // PN version
		constraints_eval, Sigma, feedback_gain,
		spacecraft_parameters, constants, solver_parameters);

	// Test
	EXPECT_EQ(constraints_eval_trans.size(), constraints_eval.size());
	for (size_t i=0; i<constraints_eval_trans.size(); i++) {
		EXPECT_TRUE(constraints_eval_trans[i].cons() > constraints_eval[i].cons());
	}
}
TEST(TestTranscription, SRDDPPathConstraints) {
	// Init DACE
	size_t nb_variables = 3;
	DA::init(3, nb_variables);

	// Init
	SolverParameters solver_parameters;
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

	vectorDA constraints_eval_trans = spectral_radius_path_inequality_transcription(
		constraints_eval, x_DA, u_DA,
		spacecraft_parameters, constants, solver_parameters);

	// Test
	EXPECT_EQ(constraints_eval_trans.size(), constraints_eval.size());
	for (size_t i=0; i<constraints_eval_trans.size(); i++) {
		EXPECT_TRUE(constraints_eval_trans[i].cons() > constraints_eval[i].cons());
	}
}
TEST(TestTranscription, SRPnPathConstraints) {
	// Init DACE
	size_t nb_variables = 3;
	DA::init(3, nb_variables);

	// Init
	SolverParameters solver_parameters;
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

	vectorDA constraints_eval_trans = spectral_radius_path_inequality_transcription( // PN version
		constraints_eval, Sigma, feedback_gain,
		spacecraft_parameters, constants, solver_parameters);

	// Test
	EXPECT_EQ(constraints_eval_trans.size(), constraints_eval.size());
	for (size_t i=0; i<constraints_eval_trans.size(); i++) {
		EXPECT_TRUE(constraints_eval_trans[i].cons() > constraints_eval[i].cons());
	}
}
TEST(TestTranscription, FODDPTerminalConstraints) {
	// Init DACE
	size_t nb_variables = 3;
	DA::init(3, nb_variables);

	// Init
	SolverParameters solver_parameters;
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

	vectorDA constraints_eval_trans = first_order_terminal_inequality_transcription(
		constraints_eval, x_DA, x_goal,
		spacecraft_parameters, constants, solver_parameters);

	// Test
	EXPECT_EQ(constraints_eval_trans.size(), constraints_eval.size());
	for (size_t i=0; i<constraints_eval_trans.size(); i++) {
		EXPECT_TRUE(constraints_eval_trans[i].cons() > constraints_eval[i].cons());
	}
}
TEST(TestTranscription, FOPnTerminalConstraints) {
	// Init DACE
	size_t nb_variables = 3;
	DA::init(3, nb_variables);

	// Init
	SolverParameters solver_parameters;
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
	matrixdb feedback_gain(solver_parameters.Nu(), solver_parameters.Nx(), 1.0);

	vectorDA constraints_eval_trans = first_order_terminal_inequality_transcription( // PN version
		constraints_eval, Sigma,
		spacecraft_parameters, constants, solver_parameters);

	// Test
	EXPECT_EQ(constraints_eval_trans.size(), constraints_eval.size());
	for (size_t i=0; i<constraints_eval_trans.size(); i++) {
		EXPECT_TRUE(constraints_eval_trans[i].cons() > constraints_eval[i].cons());
	}
}
TEST(TestTranscription, SRDDPTerminalConstraints) {
	// Init DACE
	size_t nb_variables = 3;
	DA::init(3, nb_variables);

	// Init
	SolverParameters solver_parameters;
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

	vectorDA constraints_eval_trans = spectral_radius_terminal_inequality_transcription(
		constraints_eval, x_DA, x_goal,
		spacecraft_parameters, constants, solver_parameters);

	// Test
	EXPECT_EQ(constraints_eval_trans.size(), constraints_eval.size());
	for (size_t i=0; i<constraints_eval_trans.size(); i++) {
		EXPECT_TRUE(constraints_eval_trans[i].cons() > constraints_eval[i].cons());
	}
}
TEST(TestTranscription, SRPnTerminalConstraints) {
	// Init DACE
	size_t nb_variables = 3;
	DA::init(3, nb_variables);

	// Init
	SolverParameters solver_parameters;
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

	vectorDA constraints_eval_trans = spectral_radius_terminal_inequality_transcription( // PN version
		constraints_eval, Sigma,
		spacecraft_parameters, constants, solver_parameters);

	// Test
	EXPECT_EQ(constraints_eval_trans.size(), constraints_eval.size());
	for (size_t i=0; i<constraints_eval_trans.size(); i++) {
		EXPECT_TRUE(constraints_eval_trans[i].cons() > constraints_eval[i].cons());
	}
}
