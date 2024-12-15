/**Tu peux me joindre au 07 50 87 93 41 aujourd'hui à partir de 14h30.
	test_dynamics.cpp

	Purpose: Test of the implementation of the dynamics
	methods.

	@author Thomas Caleb

	@version 1.0 13/12/2024
*/

#include "gtest/gtest.h"
#include "../include/soda.h"

using namespace DACE;
using namespace std;

TEST(TestDynamics, AccTBPSUNLT) {
	// Init DACE
	size_t nb_variables = 3;
	DA::init(2, nb_variables);

	// Init
	Constants constants;
	vectorDA x0(SIZE_VECTOR + 2, 1.0 + 2.0 * DA(2));
	vectorDA u(3, 1.0/ constants.thrustu() + DA(1));
	double t = 0; double dt = 1.0;
	SpacecraftParameters sp(constants);
	SolverParameters solver_p;
	vectorDA xf = acceleration_tbp_SUN_lt(x0, u, t, sp, constants, solver_p);

	// Tests
	EXPECT_EQ(xf.size(), x0.size());

	// For doubles
	vectordb xf_db = acceleration_tbp_SUN_lt(
		x0.cons(), u.cons(), t, sp, constants, solver_p);

	// Tests
	EXPECT_EQ(xf_db.size(), x0.cons().size());
	EXPECT_EQ(xf_db, xf.cons());
}
TEST(TestDynamics, AccTBPEARTHLT) {
	// Init DACE
	size_t nb_variables = 3;
	DA::init(2, nb_variables);

	// Init
	Constants constants;
	vectorDA x0(SIZE_VECTOR + 2, 0.5 + 2.0 * DA(2));
	vectorDA u(3, 1.0 / constants.thrustu() + DA(1));
	double t = 0; double dt = 1.0;
	SpacecraftParameters sp(constants);
	SolverParameters solver_p;
	vectorDA xf = acceleration_tbp_EARTH_lt(x0, u, t, sp, constants, solver_p);

	// Tests
	EXPECT_EQ(xf.size(), x0.size());

	// For doubles
	vectordb xf_db = acceleration_tbp_EARTH_lt(
		x0.cons(), u.cons(), t, sp, constants, solver_p);

	// Tests
	EXPECT_EQ(xf_db.size(), x0.cons().size());
	EXPECT_EQ(xf_db, xf.cons());
}
TEST(TestDynamics, AccCR3BPLT) {
	// Init DACE
	size_t nb_variables = 3;
	DA::init(2, nb_variables);

	// Init
	Constants constants;
	vectorDA x0(SIZE_VECTOR + 2, 0.5 + 2.0 * DA(2));
	vectorDA u(3, 1.0 / constants.thrustu() + DA(1));
	double t = 0; double dt = 1.0;
	SpacecraftParameters sp(constants);
	SolverParameters solver_p;
	vectorDA xf = acceleration_cr3bp_lt(x0, u, t, sp, constants, solver_p);

	// Tests
	EXPECT_EQ(xf.size(), x0.size());

	// For doubles
	vectordb xf_db = acceleration_cr3bp_lt(
		vectordb(x0.cons()), vectordb(u.cons()), t, sp, constants, solver_p);

	// Tests
	EXPECT_EQ(xf_db.size(), x0.cons().size());
	EXPECT_EQ(xf_db, xf.cons());
}
TEST(TestDynamics, kep2cart) {
	// Init DACE

	// Init
	double mu = 1;
	vectordb kep_state_vector(SIZE_VECTOR + 2, 0.5);

	// Evaluate
	vectordb cart_state_vector = kep_2_cart(kep_state_vector, mu);
	
	// Tests
	EXPECT_EQ(cart_state_vector.size(), kep_state_vector.size());
	for (size_t i = SIZE_VECTOR; i < cart_state_vector.size(); i++) {
		EXPECT_EQ(cart_state_vector[i], kep_state_vector[i]); // The non-cart values do not change
	}
}
TEST(TestDynamics, equi2kep) {
	// Init DACE

	// Init
	vectordb equi_state_vector(SIZE_VECTOR + 2, 0.5);

	// Evaluate
	vectordb kep_state_vector = equi_2_kep(equi_state_vector);

	// Tests
	EXPECT_EQ(equi_state_vector.size(), kep_state_vector.size());
	for (size_t i = SIZE_VECTOR; i < equi_state_vector.size(); i++) {
		EXPECT_EQ(equi_state_vector[i], kep_state_vector[i]); // The non-cart values do not change
	}
}
TEST(TestDynamics, kep2equi) {
	// Init DACE

	// Init
	vectordb kep_state_vector(SIZE_VECTOR + 2, 0.5);

	// Evaluate
	vectordb equi_state_vector = kep_2_equi(kep_state_vector);

	// Tests
	EXPECT_EQ(equi_state_vector.size(), kep_state_vector.size());
	for (size_t i = SIZE_VECTOR; i < equi_state_vector.size(); i++) {
		EXPECT_EQ(equi_state_vector[i], kep_state_vector[i]); // The non-cart values do not change
	}
}
TEST(TestDynamics, RTN2cart) {
	// Init DACE

	// Init
	double mu = 1;
	vectordb RTN_vector(SIZE_VECTOR/2, 0.2);
	vectordb cart_state_vector(SIZE_VECTOR + 2, 0.5);
	RTN_vector[1] = 5.0; cart_state_vector[4] = -1.0;

	// Evaluate
	vectordb cart_vector = RTN_2_cart(
		RTN_vector, cart_state_vector);

	// Tests
	EXPECT_EQ(RTN_vector.size(), cart_vector.size());
}
