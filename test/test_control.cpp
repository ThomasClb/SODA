/**
	test_control.cpp

	Purpose: Test of the implementation of the Control class.

	@author Thomas Caleb

	@version 1.0 05/09/2024
*/

#include "gtest/gtest.h"
#include "../include/soda.h"

using namespace DACE;
using namespace std;

TEST(TestControl, ConstructorGetSet) {
	// Init DACE
	size_t nb_variables = 3;
	DA::init(2, nb_variables);

	// Init

	// Empty
	Control<DA> control_empty;

	// Size
	size_t size_control(3);
	Control<DA> control_size(size_control);

	// Size + value
	DA val_control(1.0+DA(1));
	Control<DA> control_val(size_control, val_control);

	// From vector
	Control<DA> control_vect(control_val.nominal_control());

	// From control
	size_t size_state(6);
	double val_gain(2*val_control.cons());
	matrixdb feedback_gain(size_control, size_state, val_gain);
	Control<DA> control_to_copy(control_val.nominal_control());
	control_to_copy.set_feedback_gain(feedback_gain);
	Control<DA> control_copy(control_to_copy);

	// Tests

	// Empty
	EXPECT_EQ(control_empty.nominal_control().size(), 0);
	EXPECT_EQ(control_empty.feedback_gain().ncols(), 0);
	EXPECT_EQ(control_empty.feedback_gain().nrows(), 0);

	// Size
	EXPECT_EQ(control_size.nominal_control().size(), size_control);
	EXPECT_EQ(control_size.feedback_gain().ncols(), 0);
	EXPECT_EQ(control_size.feedback_gain().nrows(), 0);

	// Size + value
	EXPECT_EQ(control_val.nominal_control().size(), size_control);
	for (size_t i=0; i<size_control;i++) {
		EXPECT_EQ((control_val.nominal_control()[i] - val_control).cons(), 0.0);
	}
	EXPECT_EQ(control_val.feedback_gain().ncols(), 0);
	EXPECT_EQ(control_val.feedback_gain().nrows(), 0);

	// Vector
	EXPECT_EQ(control_vect.nominal_control().size(), size_control);
	for (size_t i=0; i<size_control;i++) {
		EXPECT_EQ((control_vect.nominal_control()[i] - val_control).cons(), 0.0);
	}
	EXPECT_EQ(control_vect.feedback_gain().ncols(), 0);
	EXPECT_EQ(control_vect.feedback_gain().nrows(), 0);

	// Copy
	EXPECT_EQ(control_copy.nominal_control().size(), size_control);
	for (size_t i=0; i<size_control;i++) {
		EXPECT_EQ((control_copy.nominal_control()[i] - val_control).cons(), 0.0);
	}
	EXPECT_EQ(control_copy.feedback_gain().ncols(), size_state);
	EXPECT_EQ(control_copy.feedback_gain().nrows(), size_control);
	for (size_t i=0; i<size_control;i++) {
		for (size_t j=0; j<size_state;j++) {
			EXPECT_EQ((control_copy.feedback_gain() - feedback_gain).at(i,j), 0.0);
		}
	}
}
TEST(TestControl, Operators) {
	// Init DACE
	size_t nb_variables = 3;
	DA::init(2, nb_variables);

	// Init
	size_t size_control(3), size_state(6);
	DA val_control(1.0+DA(1));
	double val_gain(2*val_control.cons());
	matrixdb feedback_gain(size_control, size_state, val_gain);
	Control<DA> control(size_control);
	for (size_t i=0; i<size_control; i++) {
		control[i] = val_control;
	}
	control.set_feedback_gain(feedback_gain);
	Control<DA> control_copy;
	control_copy = control;
	controldb control_cons = control.cons();

	// Tests
	EXPECT_EQ(control_copy.nominal_control().size(), control.nominal_control().size());
	for (size_t i=0; i<size_control;i++) {
		EXPECT_EQ((control_copy[i] - control.nominal_control()[i]).cons(), 0.0);
	}
	EXPECT_EQ(control_copy.feedback_gain().ncols(), control.feedback_gain().ncols());
	EXPECT_EQ(control_copy.feedback_gain().nrows(), control.feedback_gain().nrows());
	for (size_t i=0; i<size_control;i++) {
		for (size_t j=0; j<size_state;j++) {
			EXPECT_EQ((control_copy.feedback_gain() - control.feedback_gain()).at(i,j), 0.0);
		}
	}
	EXPECT_EQ(control_cons.nominal_control().size(), control.nominal_control().size());
	for (size_t i=0; i<size_control;i++) {
		EXPECT_EQ((control_cons[i] - control.nominal_control()[i].cons()), 0.0);
	}
	EXPECT_EQ(control_cons.feedback_gain().ncols(), control.feedback_gain().ncols());
	EXPECT_EQ(control_cons.feedback_gain().nrows(), control.feedback_gain().nrows());
	for (size_t i=0; i<size_control;i++) {
		for (size_t j=0; j<size_state;j++) {
			EXPECT_EQ((control_cons.feedback_gain() - control.feedback_gain()).at(i,j), 0.0);
		}
	}
}
TEST(TestControl, IO) {
	// Init
	size_t size_control(3), size_state(6);
	double val_control(1.0);
	double val_gain(2*val_control);
	matrixdb feedback_gain(size_control, size_state, val_gain);
	Control<double> control(size_control);
	for (size_t i=0; i<size_control; i++) {
		control[i] = val_control;
	}
	control.set_feedback_gain(feedback_gain);
	Control<double> control_copy;

	// Open file
	string file_name_("../data/robust_trajectory/test_control.dat");
	ofstream ofs(file_name_);

	// Store the object to file
	ofs << control;

	ofs.close();

	// Open file
	ifstream ifs(file_name_);

	// Load data
	ifs >> control_copy;

	ifs.close();

	// Tests
	EXPECT_EQ(control_copy.nominal_control().size(), control.nominal_control().size());
	for (size_t i=0; i<size_control;i++) {
		EXPECT_EQ((control_copy[i] - control.nominal_control()[i]), 0.0);
	}
	EXPECT_EQ(control_copy.feedback_gain().ncols(), control.feedback_gain().ncols());
	EXPECT_EQ(control_copy.feedback_gain().nrows(), control.feedback_gain().nrows());
	for (size_t i=0; i<size_control;i++) {
		for (size_t j=0; j<size_state;j++) {
			EXPECT_EQ((control_copy.feedback_gain() - control.feedback_gain()).at(i,j), 0.0);
		}
	}
}
TEST(TestControl, FirstGuess) {
	// Init
	size_t size_control(3), size_state(6);
	double val_control(1.0);
	double val_gain(2*val_control);
	matrixdb feedback_gain(size_control, size_state, val_gain);
	Control<double> control(size_control);
	for (size_t i=0; i<size_control; i++) {
		control[i] = val_control;
	}
	control.set_feedback_gain(feedback_gain);
	Control<double> control_copy = make_first_guess(control.nominal_control(), val_gain, size_state);

	// Tests
	EXPECT_EQ(control_copy.nominal_control().size(), control.nominal_control().size());
	for (size_t i=0; i<size_control;i++) {
		EXPECT_EQ((control_copy[i] - control.nominal_control()[i]), 0.0);
	}
	EXPECT_EQ(control_copy.feedback_gain().ncols(), control.feedback_gain().ncols());
	EXPECT_EQ(control_copy.feedback_gain().nrows(), control.feedback_gain().nrows());
	for (size_t i=0; i<size_control;i++) {
		for (size_t j=0; j<size_state;j++) {
			EXPECT_EQ((control_copy.feedback_gain() - control.feedback_gain()).at(i,j), 0.0);
		}
	}
}
