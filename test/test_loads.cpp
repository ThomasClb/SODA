/**
	test_loads.cpp

	Purpose: Test of the implementation of the Low-Order Automatic Domain Splitting
    methods.
    
	@author Thomas Caleb

    @version 1.0 23/10/2024
*/

#include "gtest/gtest.h"
#include "../include/soda.h"

using namespace DACE;
using namespace std;

TEST(TestLoads, Eigenvalues) {
	// Init
	size_t size_state(6);
	size_t size_control(3);
	double val_state(0.1);
	double val_control(0.15);
	matrixdb Sigma_x(size_state, size_state, val_control);
	matrixdb feedback_gain(size_control, size_state, val_state);
	for (size_t i=0; i<size_state; i++) {
		Sigma_x.at(i, i) += 1.0 + i;
		for (size_t j=0; j<size_control; j++) {
			feedback_gain.at(i,j) *= i*j*1.0;
		}
	}
	matrixdb Sigma_u(feedback_gain*Sigma_x*feedback_gain.transpose());	
	pair<vectordb, matrixdb> eigen = get_eigenvalues(Sigma_x, feedback_gain);
	matrixdb D = make_diag_matrix_(eigen.first);
	matrixdb P = eigen.second;
	matrixdb Sigma = P.transpose()*D*P;
	
	// Tests
	EXPECT_NEAR(frobenius_norm_(
		Sigma.submat(size_state - 1, size_state - 1) - Sigma_x), 0, EPS);
	EXPECT_NEAR(frobenius_norm_(
		Sigma.submat(
			size_state, size_state,
			size_control + size_state - 1, size_control + size_state - 1) - Sigma_u), 0, EPS);
}
TEST(TestLoads, Scale) {
	// Init
	size_t size_state(6);
	size_t size_control(3);

	// Init DACE
	size_t nb_variables = size_state + size_control;
	DA::init(2, nb_variables);

	// Init
	double transcription_beta = 5e-2;
	double quantile = inv_chi_2_cdf(size_state, 1 - transcription_beta);
	matrixdb Sigma_x(size_state, size_state, 0.0);
	matrixdb feedback_gain(size_control, size_state, 0.0);
	for (size_t i=0; i<size_state; i++) {
		Sigma_x.at(i, i) += 1.0 + i;
	}
	for (size_t i=0; i<size_control; i++) {
		feedback_gain.at(i, i) = 1.0;
	}
	matrixdb Sigma_u(feedback_gain*Sigma_x*feedback_gain.transpose());	
	vectorDA y(nb_variables, 1.0);
	for (int i=0; i<nb_variables; i++) {
		y[i] += DA(i + 1);
	}
	pair<vectorDA, vectordb> scaling = scale(y, Sigma_x, feedback_gain, transcription_beta);
	vectordb y_scaled_eval = scaling.first.eval(vectordb(nb_variables, 1.0));

	// Tests
	for (size_t i=0; i<size_state; i++) {
		EXPECT_EQ(y_scaled_eval[i], 1.0 + sqrt(quantile*Sigma_x.at(i, i)));
	}
	for (size_t i=0; i<size_control; i++) {
		EXPECT_EQ(y_scaled_eval[i + size_state], 1.0 + sqrt(quantile*Sigma_u.at(i, i)));
	}
	for (size_t i=0; i<size_state; i++) {
		EXPECT_EQ(scaling.second[i], quantile*Sigma_x.at(i, i));
	}
	for (size_t i=0; i<size_control; i++) {
		EXPECT_EQ(scaling.second[i + size_state], quantile*Sigma_u.at(i, i));
	}
}
TEST(TestLoads, NLIndex) {
	// Init DACE
	size_t nb_variables = 1;
	DA::init(2, nb_variables);

	// Tests
	vectorDA y(nb_variables);
	y[0] = 1 + DA(1) + sqr(DA(1));
	EXPECT_EQ(nl_index(y, vectordb(nb_variables, 1.0)), 2);
	y[0] = 1 + DA(1) + 0.5*sqr(DA(1));
	EXPECT_EQ(nl_index(y, vectordb(nb_variables, 1.0)), 1);
	y[0] = 1 + DA(1) + 0*sqr(DA(1));
	EXPECT_EQ(nl_index(y, vectordb(nb_variables, 1.0)), 0);
}
TEST(TestLoads, NLIndexScaled) {
	// Init
	size_t size_state(6);
	size_t size_control(3);

	// Init DACE
	size_t nb_variables = size_state + size_control;
	DA::init(2, nb_variables);

	// Init
	double transcription_beta = 5e-2;
	double quantile = inv_chi_2_cdf(size_state, 1 - transcription_beta);
	matrixdb Sigma_x(size_state, size_state, 0.0);
	matrixdb feedback_gain(size_control, size_state, 0.0);
	for (size_t i=0; i<size_state; i++) {
		Sigma_x.at(i, i) += 1.0 + i;
	}
	for (size_t i=0; i<size_control; i++) {
		feedback_gain.at(i, i) = 1.0;
	}
	matrixdb Sigma_u(feedback_gain*Sigma_x*feedback_gain.transpose());	
	vectorDA y(nb_variables, 1.0);
	for (int i=0; i<nb_variables; i++) {
		y[i] += DA(i + 1);
	}

	double nli_0 = nl_index(y, Sigma_x, feedback_gain, transcription_beta);
	y[0] += DA(1)*DA(1);
	double nli_1 = nl_index(y, Sigma_x, feedback_gain, transcription_beta);

	// Tests
	EXPECT_EQ(nli_0, 0);
	EXPECT_TRUE(nli_1 != 0);
}
TEST(TestLoads, NLIndexTrajectory) {
	// Init
	int size_state(6);
	int size_control(3);

	// Init DACE
	size_t nb_variables = size_state + size_control;
	DA::init(2, nb_variables);

	// Init
	int N = 2;
	vector<statedb> list_x(N + 1);
	vector<controldb> list_u(N);
	vector<vectorDA> list_dynamics_eval(N);
	double transcription_beta = 5e-2;
	for (size_t k=0; k<N+1; k++) {
		matrixdb Sigma_x(size_state, size_state, 1.0);
		for (size_t i=0; i<size_state; i++) {Sigma_x.at(i, i) += 1.0 + 1.0*i;}
		vectorDA y(size_state, 1.0);
		for (int i=0; i<size_state; i++) {
			y[i] += i*1.0 + 0.1*DA(i + 1);
		}
		for (int i=0; i<size_control; i++) {
			y[i] += 0.5*DA(size_state + i + 1);
		}
		if (k < N)
			list_dynamics_eval[k] = y;
		statedb state(y.cons());
		state.set_Sigma(Sigma_x);
		state.set_der_dynamics(matrixdb(size_state, nb_variables, 0.0));
		list_x[k] = state;
	}
	for (size_t k = 0; k < N; k++)	{
		controldb control(vectordb(size_control, 1.0));
		matrixdb feedback_gain(size_control, size_state, 0.1);
		for (size_t i=0; i<size_control; i++) {feedback_gain.at(i, i) = 1.0;}
		control.set_feedback_gain(feedback_gain);
		list_u[k] = control;
	}

	vectordb list_nli = nl_index(
		list_dynamics_eval, list_x, list_u, transcription_beta);

	// Tests
	EXPECT_EQ(list_nli.size(), list_dynamics_eval.size());
	for (size_t i=0; i<list_nli.size(); i++) {
		EXPECT_EQ(list_nli[0], 0.0);
	}
	
}