/**
	test_state.cpp

	Purpose: Test of the implementation of the State class.

	@author Thomas Caleb

	@version 1.0 09/09/2024
*/

#include "gtest/gtest.h"
#include "../include/soda.h"

using namespace DACE;
using namespace std;

TEST(TestState, ConstructorGetSet) {
	// Init DACE
	size_t nb_variables = 3;
	DA::init(2, nb_variables);

	// Init

	// Empty
	State<DA> state_empty;

	// Size
	size_t size_state(6);
	State<DA> state_size(size_state);

	// Size + value
	DA val_state(1.0+DA(1));
	State<DA> state_val(size_state, val_state);

	// State vector
	State<DA> state_vect(state_val.nominal_state());

	// From state
	size_t size_control(3);
	double val_gain(2*val_state.cons());
	matrixdb der_dynamics(size_state, size_state + size_control, val_gain);
	matrixdb Sigma(size_state, size_state, val_gain);
	for (size_t i=0; i<size_state; i++) {
		Sigma.at(i, i) += 1.0;
	}
	State<DA> state_to_copy(state_val.nominal_state());
	state_to_copy.set_Sigma(Sigma);
	state_to_copy.set_der_dynamics(der_dynamics);
	State<DA> state_copy(state_to_copy);

	// Tests

	// Empty
	EXPECT_EQ(state_empty.nominal_state().size(), 0);
	EXPECT_EQ(state_empty.Sigma().ncols(), 0);
	EXPECT_EQ(state_empty.Sigma().nrows(), 0);
	EXPECT_EQ(state_empty.der_dynamics().ncols(), 0);
	EXPECT_EQ(state_empty.der_dynamics().nrows(), 0);

	// Size
	EXPECT_EQ(state_size.nominal_state().size(), size_state);
	EXPECT_EQ(state_size.Sigma().ncols(), 0);
	EXPECT_EQ(state_size.Sigma().nrows(), 0);
	EXPECT_EQ(state_size.der_dynamics().ncols(), 0);
	EXPECT_EQ(state_size.der_dynamics().nrows(), 0);

	// Size + value
	EXPECT_EQ(state_val.nominal_state().size(), size_state);
	for (size_t i=0; i<size_state;i++) {
		EXPECT_EQ((state_val.nominal_state()[i] - val_state).cons(), 0.0);
	}
	EXPECT_EQ(state_val.Sigma().ncols(), 0);
	EXPECT_EQ(state_val.Sigma().nrows(), 0);
	EXPECT_EQ(state_val.der_dynamics().ncols(), 0);
	EXPECT_EQ(state_val.der_dynamics().nrows(), 0);

	// Vector
	EXPECT_EQ(state_vect.nominal_state().size(), size_state);
	for (size_t i=0; i<size_state;i++) {
		EXPECT_EQ((state_vect.nominal_state()[i] - val_state).cons(), 0.0);
	}
	EXPECT_EQ(state_vect.Sigma().ncols(), 0);
	EXPECT_EQ(state_vect.Sigma().nrows(), 0);
	EXPECT_EQ(state_vect.der_dynamics().ncols(), 0);
	EXPECT_EQ(state_vect.der_dynamics().nrows(), 0);

	// Copy
	EXPECT_EQ(state_copy.nominal_state().size(), size_state);
	for (size_t i=0; i<size_state;i++) {
		EXPECT_EQ((state_copy.nominal_state()[i] - state_to_copy.nominal_state()[i]).cons(), 0.0);
	}
	EXPECT_EQ(state_copy.Sigma().ncols(), state_to_copy.Sigma().ncols());
	EXPECT_EQ(state_copy.Sigma().nrows(), state_to_copy.Sigma().nrows());
	for (size_t i=0; i<size_state;i++) {
		for (size_t j=0; j<size_state;j++) {
			EXPECT_EQ((state_copy.Sigma().at(i,j) - state_to_copy.Sigma().at(i,j)), 0.0);
		}
	}
	EXPECT_EQ(state_copy.der_dynamics().ncols(), state_to_copy.der_dynamics().ncols());
	EXPECT_EQ(state_copy.der_dynamics().nrows(), state_to_copy.der_dynamics().nrows());
	for (size_t i=0; i<size_state;i++) {
		for (size_t j=0; j<size_state+size_control;j++) {
			EXPECT_EQ((state_copy.der_dynamics().at(i,j) - state_to_copy.der_dynamics().at(i,j)), 0.0);
		}
	}
}
TEST(TestState, Operators) {
	// Init DACE
	size_t nb_variables = 3;
	DA::init(2, nb_variables);

	// Init
	size_t size_state(6), size_control(3);
	DA val_state(1.0+DA(1));
	State<DA> state_val(size_state);
	for (size_t i=0; i<size_state; i++) {
		state_val[i] = val_state;
	}
	double val_gain(2*val_state.cons());
	matrixdb der_dynamics(size_state, size_state + size_control, val_gain);
	matrixdb Sigma(size_state, size_state, val_gain);
	for (size_t i=0; i<size_state; i++) {
		Sigma.at(i, i) += 1.0;
	}
	State<DA> state_to_copy(state_val.nominal_state());
	state_to_copy.set_Sigma(Sigma);
	state_to_copy.set_der_dynamics(der_dynamics);
	State<DA> state_copy(state_to_copy);
	statedb state_cons(state_to_copy.cons());

	// Tets
	EXPECT_EQ(state_copy.nominal_state().size(), size_state);
	for (size_t i=0; i<size_state;i++) {
		EXPECT_EQ((state_copy[i] - state_to_copy.nominal_state()[i]).cons(), 0.0);
	}
	EXPECT_EQ(state_copy.Sigma().ncols(), state_to_copy.Sigma().ncols());
	EXPECT_EQ(state_copy.Sigma().nrows(), state_to_copy.Sigma().nrows());
	for (size_t i=0; i<size_state;i++) {
		for (size_t j=0; j<size_state;j++) {
			EXPECT_EQ((state_copy.Sigma().at(i,j) - state_to_copy.Sigma().at(i,j)), 0.0);
		}
	}
	EXPECT_EQ(state_copy.der_dynamics().ncols(), state_to_copy.der_dynamics().ncols());
	EXPECT_EQ(state_copy.der_dynamics().nrows(), state_to_copy.der_dynamics().nrows());
	for (size_t i=0; i<size_state;i++) {
		for (size_t j=0; j<size_state+size_control;j++) {
			EXPECT_EQ((state_copy.der_dynamics().at(i,j) - state_to_copy.der_dynamics().at(i,j)), 0.0);
		}
	}

	EXPECT_EQ(state_cons.nominal_state().size(), size_state);
	for (size_t i=0; i<size_state;i++) {
		EXPECT_EQ((state_cons[i] - state_to_copy.nominal_state()[i].cons()), 0.0);
	}
	EXPECT_EQ(state_cons.Sigma().ncols(), state_to_copy.Sigma().ncols());
	EXPECT_EQ(state_cons.Sigma().nrows(), state_to_copy.Sigma().nrows());
	for (size_t i=0; i<size_state;i++) {
		for (size_t j=0; j<size_state;j++) {
			EXPECT_EQ((state_cons.Sigma().at(i,j) - state_to_copy.Sigma().at(i,j)), 0.0);
		}
	}
	EXPECT_EQ(state_cons.der_dynamics().ncols(), state_to_copy.der_dynamics().ncols());
	EXPECT_EQ(state_cons.der_dynamics().nrows(), state_to_copy.der_dynamics().nrows());
	for (size_t i=0; i<size_state;i++) {
		for (size_t j=0; j<size_state+size_control;j++) {
			EXPECT_EQ((state_cons.der_dynamics().at(i,j) - state_to_copy.der_dynamics().at(i,j)), 0.0);
		}
	}
}
TEST(TestState, IO) {
	// Init
	size_t size_state(6);
	double val_state(1.0);
	size_t size_control(3);
	double val_gain(2*val_state);
	matrixdb der_dynamics(size_state, size_state + size_control, val_gain);
	matrixdb Sigma(size_state, size_state, val_gain);
	for (size_t i=0; i<size_state; i++) {
		Sigma.at(i, i) += 1.0;
	}
	State<double> state(vectordb(size_state, val_state));
	state.set_Sigma(Sigma); state.set_der_dynamics(der_dynamics); 
	State<double> state_copy;

	// Open file
	string file_name_("../data/robust_trajectory/test_state.dat");
	ofstream ofs(file_name_);

	// Store the object to file
	ofs << state;

	ofs.close();

	// Open file
	ifstream ifs(file_name_);

	// Load data
	ifs >> state_copy;

	ifs.close();

	// Tests
	EXPECT_EQ(state_copy.nominal_state().size(), size_state);
	for (size_t i=0; i<size_state;i++) {
		EXPECT_EQ((state_copy.nominal_state()[i] - state.nominal_state()[i]), 0.0);
	}
	EXPECT_EQ(state_copy.Sigma().ncols(), state.Sigma().ncols());
	EXPECT_EQ(state_copy.Sigma().nrows(), state.Sigma().nrows());
	for (size_t i=0; i<size_state;i++) {
		for (size_t j=0; j<size_state;j++) {
			EXPECT_EQ((state_copy.Sigma().at(i,j) - state.Sigma().at(i,j)), 0.0);
		}
	}
	EXPECT_EQ(state_copy.der_dynamics().ncols(), state.der_dynamics().ncols());
	EXPECT_EQ(state_copy.der_dynamics().nrows(), state.der_dynamics().nrows());
	for (size_t i=0; i<size_state;i++) {
		for (size_t j=0; j<size_state+size_control;j++) {
			EXPECT_EQ((state_copy.der_dynamics().at(i,j) - state.der_dynamics().at(i,j)), 0.0);
		}
	}
}
TEST(TestState, MakeInitialState) {
	// Init
	size_t size_state(6);
	double val_state(1.0);
	size_t size_control(3);
	matrixdb der_dynamics(size_state, size_state + size_control, 0.0);
	matrixdb Sigma(size_state, size_state, 0.0);
	for (size_t i=0; i<size_state; i++) {
		Sigma.at(i, i) += 1.0;
	}
	State<double> state(vectordb(size_state, val_state));
	state.set_Sigma(Sigma); state.set_der_dynamics(der_dynamics); 
	State<double> state_copy = make_initial_state(
		state.nominal_state(), vectordb(size_state, 1.0), size_control);

	// Tests
	EXPECT_EQ(state_copy.nominal_state().size(), size_state);
	for (size_t i=0; i<size_state;i++) {
		EXPECT_EQ((state_copy.nominal_state()[i] - state.nominal_state()[i]), 0.0);
	}
	EXPECT_EQ(state_copy.Sigma().ncols(), state.Sigma().ncols());
	EXPECT_EQ(state_copy.Sigma().nrows(), state.Sigma().nrows());
	for (size_t i=0; i<size_state;i++) {
		for (size_t j=0; j<size_state;j++) {
			EXPECT_EQ((state_copy.Sigma().at(i,j) - state.Sigma().at(i,j)), 0.0);
		}
	}
	EXPECT_EQ(state_copy.der_dynamics().ncols(), state.der_dynamics().ncols());
	EXPECT_EQ(state_copy.der_dynamics().nrows(), state.der_dynamics().nrows());
	for (size_t i=0; i<size_state;i++) {
		for (size_t j=0; j<size_state+size_control;j++) {
			EXPECT_EQ((state_copy.der_dynamics().at(i,j) - state.der_dynamics().at(i,j)), 0.0);
		}
	}
}
