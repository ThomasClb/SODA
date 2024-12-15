/**
	dynamics.cpp

	Purpose: Implementation of the Dynamics class
	methods.

	@author Thomas Caleb

	@version 2.0 13/12/2023
*/

#ifndef DEF_DYNAMIC
#define DEF_DYNAMIC

#pragma once

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <functional>

#include <dace/dace_s.h>

#include "settings.h"
#include "constants.h"
#include "parameters.h"
#include "integration.h"


// Define the 5 types of functions

// vDA-vDA-SpacecraftParameters-SolverParameters->vDA

// Returns the next state given
// the current state, the control and the parameters.
using dynFunction = std::function<DACE::vectorDA(
	DACE::vectorDA const&, DACE::vectorDA const&,
	SpacecraftParameters const&, Constants const&,
	SolverParameters const&)>;

// Returns the path inequality contraints given
// the current state, the control and the parameters.
using ineqFunction = std::function<DACE::vectorDA(
	DACE::vectorDA const&, DACE::vectorDA const&,
	SpacecraftParameters const&, Constants const&,
	SolverParameters const&)>;

// vDA-vDA-SpacecraftParameters-SolverParameters->DA

// Returns the cost-to-go given
// the current state, the control and the parameters.
using scFunction = std::function<DACE::DA(
	DACE::vectorDA const&, DACE::vectorDA const&, 
	SpacecraftParameters const&,
	Constants const&,
	SolverParameters const&)>;

// vDA-vdb-SpacecraftParameters-SolverParameters->DA

// Returns the terminal cost given
// the current state, the target state and the parameters.
using tcFunction = std::function<DACE::DA(
	DACE::vectorDA const&, DACE::vectordb const&,
	SpacecraftParameters const&,
	Constants const&,
	SolverParameters const&)>;

// vDA-vdb-SpacecraftParameters-SolverParameters->vDA

// Returns the terminal inequality constraints given
// the current state, the target state and the parameters.
using tineqFunction = std::function<DACE::vectorDA(
	DACE::vectorDA const&, DACE::vectordb const&,
	SpacecraftParameters const&,
	Constants const&,
	SolverParameters const&)>;

// Double versions

// vdb-vdb-SpacecraftParameters-SolverParameters->vdb

// Returns the next state given
// the current state, the control and the parameters.
using dynFunction_db = std::function<DACE::vectordb(
	DACE::vectordb const&,
	DACE::vectordb const&, SpacecraftParameters const&,
	Constants const&,
	SolverParameters const&)>;

// Returns the path inequality contraints given
// the current state, the control and the parameters.
using ineqFunction_db = std::function<DACE::vectordb(
	DACE::vectordb const&,
	DACE::vectordb const&, SpacecraftParameters const&,
	Constants const&, SolverParameters const&)>;

// vdb-vdb-SpacecraftParameters-SolverParameters->db

// Returns the cost-to-go given
// the current state, the control and the parameters.
using scFunction_db = std::function<double(
	DACE::vectordb const&, DACE::vectordb const&,
	SpacecraftParameters const&,
	Constants const&,
	SolverParameters const&)>;

// vdb-vdb-SpacecraftParameters-SolverParameters->db

// Returns the terminal cost given
// the current state, the target state and the parameters.
using tcFunction_db = std::function<double(
	DACE::vectordb const&, DACE::vectordb const&,
	SpacecraftParameters const&,
	Constants const&,
	SolverParameters const&)>;

// vdb-vdb-SpacecraftParameters-SolverParameters->vdb

// Returns the terminal inequality constraints given
// the current state, the target state and the parameters.
using tineqFunction_db = std::function<DACE::vectordb(
	DACE::vectordb const&, DACE::vectordb const&,
	SpacecraftParameters const&,
	Constants const&,
	SolverParameters const&)>;


/*

	DYNAMICS

*/

class Dynamics {

	// Attributes
protected:
	Constants constants_;

	// DA functions
	dynFunction dynamic_;
	scFunction stage_cost_;
	ineqFunction inequality_constraints_;
	tcFunction terminal_cost_;
	tineqFunction terminal_inequality_constraints_;

	// Double versions
	dynFunction_db dynamic_db_;
	scFunction_db stage_cost_db_;
	ineqFunction_db inequality_constraints_db_;
	tcFunction_db terminal_cost_db_;
	tineqFunction_db terminal_inequality_constraints_db_;

// Methods
public:
	// Constructors

	// Default constructors
	// No unit test.
	Dynamics();

	// Constructor
	// No unit test.
	Dynamics(
		Constants const& constants,
		dynFunction const& dynamic,
		scFunction const& stage_cost,
		ineqFunction const& inequality_constraints,
		tcFunction const& terminal_cost,
		tineqFunction const& terminal_inequality_constraints,
		dynFunction_db const& dynamic_db,
		scFunction_db const& stage_cost_db,
		ineqFunction_db const& inequality_constraints_db,
		tcFunction_db const& terminal_cost_db,
		tineqFunction_db const& terminal_inequality_constraints_db);

	// Copy constructor
	// No unit test.
	Dynamics(Dynamics const& dynamics);

	// Destructors
	~Dynamics();

	// Getters
	// No unit test.
	const Constants constants() const;

	const dynFunction dynamic() const;
	const scFunction stage_cost() const;
	const ineqFunction inequality_constraints() const;
	const tcFunction terminal_cost() const;
	const tineqFunction terminal_inequality_constraints() const;

	const dynFunction_db dynamic_db() const;
	const scFunction_db stage_cost_db() const;
	const ineqFunction_db inequality_constraints_db() const;
	const tcFunction_db terminal_cost_db() const;
	const tineqFunction_db terminal_inequality_constraints_db() const;
};


/*

	FUNCTIONS

*/

// Acceleration functions

// Computes the derivaties in an Sun-centered 2-body problem.
// Cartesian coordinates.
// With 3D continuous thrust.
// It takes at input [3*LU, 3*VU, MASSU, TU], [3*N], TU, SpacecraftParameters.
// It returns [3*VU, 3*VU/TU, MASSU/TU, 1].
// https://en.wikipedia.org/wiki/Two-body_problem
template<typename T>
DACE::AlgebraicVector<T> acceleration_tbp_SUN_lt(
	DACE::AlgebraicVector<T> const& x,
	DACE::AlgebraicVector<T> const& u, double const& t,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters) {
	// Unpack
	double v_e = spacecraft_parameters.ejection_velocity(); // [VU]

	// Init output
	DACE::AlgebraicVector<T> output(SIZE_VECTOR + 2);
	output[0] = x[3];
	output[1] = x[4];
	output[2] = x[5];
	output[7] = 0.0; // ToF

	// Acceleration kepler
	double mu = MU_SUN/ constants.mu();
	DACE::AlgebraicVector<T> r = x.extract(0, 2);
	T r_2 = r.dot(r);
	DACE::AlgebraicVector<T> acc_kep = -mu * r * pow(r_2, -1.5);

	// Thrust
	T inv_mass = 1 / x[6];
	DACE::AlgebraicVector<T> acc_thrust = inv_mass * u;

	// Assign
	DACE::AlgebraicVector<T> acc = acc_kep + acc_thrust;

	// dm [MASSU/LU]
	T thrust_norm = u.vnorm();
	T m_p = thrust_norm * pow(-1.0*v_e, -1);

	// Assign to output
	output[3] = acc[0];
	output[4] = acc[1];
	output[5] = acc[2];
	output[6] = m_p;
	
	return x[7] * output;
}

// Computes the derivaties in an Earth-centered 2-body problem.
// Equinoctial coordinates.
// With 3D continuous thrust using Gauss planetary equations.
// It takes at input [3*LU, 3*VU, MASSU, TU], [3*N], TU, SpacecraftParameters.
// It returns [3*VU, 3*VU/TU, MASSU/TU, 1].
// https://en.wikipedia.org/wiki/Two-body_problem
// Same formulation as in [Di Carlo et al. 2021]
// DOI: https://doi.org/10.1007/s10569-021-10007-x
template<typename T>
DACE::AlgebraicVector<T> acceleration_tbp_EARTH_lt(
	DACE::AlgebraicVector<T> const& x,
	DACE::AlgebraicVector<T> const& u, double const& t,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters) {
	// Unpack
	double v_e = spacecraft_parameters.ejection_velocity(); // [VU]
	double lu = constants.lu();
	double mu = MU_EARTH / constants.mu();
	bool with_J2 = solver_parameters.with_J2();
	T sma = x[0]; // Equinoctial elements
	T P_1 = x[1];
	T P_2 = x[2];
	T Q_1 = x[3];
	T Q_2 = x[4];
	T L = x[5];
	T u_R = u[0]; 
	T u_T = u[1];
	T u_N = u[2];

	// Init output
	DACE::AlgebraicVector<T> output(SIZE_VECTOR + 2);
	output[7] = 0.0; // ToF

	// Compute useful data
	T cos_L = cos(L);
	T sin_L = sin(L);
	T B = sqrt(1.0 - P_1 * P_1 - P_2 * P_2);
	T G = 1.0 + Q_1 * Q_1 + Q_2 * Q_2;
	T Phi_L = 1 + P_1 * cos_L + P_2 * sin_L;
	T inv_Phi_L = 1.0 / Phi_L;
	T sqrt_a_mu = sqrt(sma / mu);
	T Q_1_cos_L = Q_1 * cos_L;
	T Q_1_sin_L = Q_1 * sin_L;
	T Q_2_cos_L = Q_2 * cos_L;
	T Q_2_sin_L = Q_2 * sin_L;
	T Q_1_cos_L_m_Q_2_sin_L = Q_1_cos_L - Q_2_sin_L;
	T pert_R = u_R;
	T pert_T = u_T;
	T pert_N = u_N;

	// J2
	if (with_J2) {
		T G_2 = pow(G, 2.0);
		T J2_mag = (1.5 * J_2 * MU_EARTH * pow(R_EARTH, 2)) * pow(sma * lu * pow(B, 2) * inv_Phi_L, -4);
		J2_mag /= (constants.mu() / lu / lu) * G_2;
		T J2_R = (12 * pow(Q_1_cos_L_m_Q_2_sin_L, 2.0) - G_2) * J2_mag;
		T J2_buff = 4 * Q_1_cos_L_m_Q_2_sin_L * J2_mag;
		T J2_T = 2 * (Q_2_cos_L + Q_1_sin_L) * J2_buff;
		T J2_N = J2_buff * (2 - G);

		// Pertubating forces
		pert_R += J2_R;
		pert_T += J2_T;
		pert_N += J2_N;
	} 

	// Acceleration Gauss
	T d_sma = 2 * sma * sqrt_a_mu / B * (
		(P_2 * sin_L - P_1 * cos_L) * pert_R
		+ Phi_L * pert_T);
	T d_P_1 = B * sqrt_a_mu * (
		- cos_L * pert_R
		+ ((P_1 + sin_L) * inv_Phi_L + sin_L) * pert_T
		- P_2 * Q_1_cos_L_m_Q_2_sin_L * inv_Phi_L * pert_N);
	T d_P_2 = B * sqrt_a_mu * (
		- sin_L * pert_R
		+ ((P_2 + cos_L) * inv_Phi_L + cos_L) * pert_T
		- P_1 * Q_1_cos_L_m_Q_2_sin_L * inv_Phi_L * pert_N);
	T d_Q = 0.5 * B * sqrt_a_mu * G * inv_Phi_L * pert_N;
	T d_Q_1 = d_Q * sin_L;
	T d_Q_2 = d_Q * cos_L;
	T d_L = Phi_L * Phi_L * pow(B, -3.0) / sqrt_a_mu  // Mean motion
		- sma * sqrt_a_mu * B * inv_Phi_L * Q_1_cos_L_m_Q_2_sin_L * pert_N;

	// Thrust
	T inv_mass = 1 / x[6];
	DACE::AlgebraicVector<T> acc_thrust = inv_mass * u;

	// dm [MASSU/LU]
	T thrust_norm = u.vnorm();
	T m_p = thrust_norm * pow(-1.0 * v_e, -1);

	// Assign to output
	output[0] = d_sma;
	output[1] = d_P_1;
	output[2] = d_P_2;
	output[3] = d_Q_1;
	output[4] = d_Q_2;
	output[5] = d_L;
	output[6] = m_p; // Mass

	return x[7] * output;
}

// Computes the derivaties in an Earth-Moon circular restricted 3-body problem.
// With 3D continuous thrust.
// It takes at input [3*LU, 3*VU, MASSU, TU], [3*N], TU, SpacecraftParameters.
// It returns [3*VU, 3*VU/TU, MASSU/TU, 1].
// See [Poincaré 1892]:
// Les Méthodes Nouvelles de la Mécanique Céleste.
template<typename T>
DACE::AlgebraicVector<T>  acceleration_cr3bp_lt(
	DACE::AlgebraicVector<T>  const& state_vector,
	DACE::AlgebraicVector<T> const& u, double const& t,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters) {
	// Unpack
	T x = state_vector[0];
	T y = state_vector[1];
	T z = state_vector[2];
	T x_p = state_vector[3];
	T y_p = state_vector[4];
	T z_p = state_vector[5];
	T mass = state_vector[6];
	T s = state_vector[7]; // Period
	double v_e = spacecraft_parameters.ejection_velocity(); // [VU]

	// velocity
	DACE::AlgebraicVector<T>  output(SIZE_VECTOR + 2);
	output[0] = x_p;
	output[1] = y_p;
	output[2] = z_p;
	output[7] = 0.0; // Period is constant

	// Inversion
	double mu = constants.mu();
	T d_y_z_2 = y* y + z* z;
	T inv_r_1_3 = (1 - mu) * pow((x + mu)* (x + mu) + d_y_z_2, -1.5);
	T inv_r_2_3 = mu * pow(pow(x - (1 - mu), 2) + d_y_z_2, -1.5);

	// Thrust
	DACE::AlgebraicVector<T> acc_thrust = u / mass;

	// Acceleration
	output[3] = 2 * y_p + x - (x + mu) * inv_r_1_3 - (x - (1 - mu)) * inv_r_2_3 + acc_thrust[0];
	output[4] = -2 * x_p + y * (1 - inv_r_1_3 - inv_r_2_3) + acc_thrust[1];
	output[5] = -z * (inv_r_1_3 + inv_r_2_3) + acc_thrust[2];

	// dm [MASSU/TU]
	T thrust_norm = u.vnorm(); // [-]
	output[6] = (-1.0/ v_e) * thrust_norm;
	return s * output;
}

// Returns the next state given
// the current state, the control and the parameters.
// For a double integrator system.
// From [Lantoine Russell 2012]
// DOI: https://doi.org/10.1007/s10957-012-0038-1
// No unit test.
template<typename T>
DACE::AlgebraicVector<T> dynamic_double_integrator(
	DACE::AlgebraicVector<T> const& x, DACE::AlgebraicVector<T>const& u,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters) {
	// Init output
	DACE::AlgebraicVector<T> output(x.size(), 0.0);
	output[0] = x[0] + x[3];
	output[1] = x[1] + x[4];
	output[2] = x[2] + x[5];
	output[3] = x[3] + u[0];
	output[4] = x[4] + u[1];
	output[5] = x[5] + u[2];

	return output;
}

// Returns the next state given
// the current state, the control and the parameters.
// With acceleration acceleration_tbp_SUN_lt.
// No unit test.
template<typename T>
DACE::AlgebraicVector<T> dynamic_tbp_SUN_lt(
	DACE::AlgebraicVector<T> const& x, DACE::AlgebraicVector<T>const& u,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters) {
	return RK78(
		acceleration_tbp_SUN_lt, x, u, 0, 1.0,
		spacecraft_parameters, constants, solver_parameters);
}

// Returns the next state given
// the current state, the control and the parameters.
// With acceleration acceleration_tbp_EARTH_lt.
// No unit test.
template<typename T>
DACE::AlgebraicVector<T> dynamic_tbp_EARTH_lt(
	DACE::AlgebraicVector<T> const& x, DACE::AlgebraicVector<T>const& u,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters) {
	return RK78(
		acceleration_tbp_EARTH_lt, x, u, 0, 1.0,
		spacecraft_parameters, constants, solver_parameters);
}

// Returns the next state given
// the current state, the control and the parameters.
// With acceleration acceleration_cr3bp_lt.
// No unit test.
template<typename T>
DACE::AlgebraicVector<T> dynamic_cr3bp_lt(
	DACE::AlgebraicVector<T> const& x, DACE::AlgebraicVector<T> const& u,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters) {
	return RK78(
		acceleration_cr3bp_lt, x, u, 0, 1.0,
		spacecraft_parameters, constants, solver_parameters);
}

// Returns the cost-to-go given
// the current state, the control and the parameters.
// Double integrator system.
// No unit test.
template<typename T>
T stage_cost_double_integrator(
	DACE::AlgebraicVector<T> const& x, DACE::AlgebraicVector<T>const& u,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters) {
	return u.dot(u);
}

// Returns the cost-to-go given
// the current state, the control and the parameters.
// Homotopy between energy-optimal and fuel-optimal low-thrust.
// No unit test.
template<typename T>
T stage_cost(
	DACE::AlgebraicVector<T> const& x, DACE::AlgebraicVector<T>const& u,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters) {
	// Unpack parameters
	double tol = solver_parameters.DDP_tol();
	unsigned int N = solver_parameters.N();
	double gain = solver_parameters.stage_cost_gain(); // [-]
	double mass_leak = solver_parameters.mass_leak(); // [-]
	double homotopy_coefficient = solver_parameters.homotopy_coefficient(); // [-]
	double delta = solver_parameters.huber_loss_coefficient(); // [-]
	double ToF = solver_parameters.ToF(); // [-]
	double T_max = spacecraft_parameters.thrust(); // [THRUSTU]

	// Get time ratio
	double dt_avg = ToF / N; // [TU]
	T dt = x[SIZE_VECTOR + 1];

	// NRJ cost to go
	DACE::AlgebraicVector<T>  u_norm = u* (dt * (1 / (dt_avg * T_max))); // [-]
	T NRJ_stage_cost = u_norm.dot(u_norm); // [-]

	// Return cost
	if (homotopy_coefficient == 0.0)
		return (0.5 * gain) * NRJ_stage_cost;

	// Pseudo-Huber loss
	NRJ_stage_cost = NRJ_stage_cost / pow(delta, 2.0); // Normalise NRJ optimal term
	mass_leak = mass_leak / pow(delta, 2.0);
	T fuel_stage_cost = sqrt(1.0 + NRJ_stage_cost + mass_leak) - 1.0;

	// Homotopy
	double c_1 = 0.5 * gain * delta * (1.0 - homotopy_coefficient);
	double c_2 = gain * homotopy_coefficient;
	return c_1 * NRJ_stage_cost + c_2 * fuel_stage_cost;
}

// Returns the path inequality contraints given
// the current state, the control and the parameters.
// Null.
// No unit test.
template<typename T>
DACE::AlgebraicVector<T> inequality_constraints_null(
	DACE::AlgebraicVector<T> const& x, DACE::AlgebraicVector<T>const& u,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters) {	
	return DACE::AlgebraicVector<T>(0);
}

// Returns the path inequality contraints given
// the current state, the control and the parameters.
// Mass constraints, and thrust constraints.
// For low-thrust tbp
// Nineq = 2.
// No unit test.
template<typename T>
DACE::AlgebraicVector<T> inequality_constraints_tbp_lt(
	DACE::AlgebraicVector<T> const& x, DACE::AlgebraicVector<T>const& u,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters) {
	// Unpack parameters
	double T_max = spacecraft_parameters.thrust(); // [THRUSTU]
	double dry_mass = spacecraft_parameters.dry_mass(); // [MASSU]
	double initial_mass = spacecraft_parameters.initial_mass(); // [MASSU]

	// Init
	DACE::AlgebraicVector<T> output; output.reserve(1 + 1);

	// Thrust (1)
	T T_const = u.dot(u) - T_max * T_max; // [THRUSTU**2]
	output.push_back(T_const);

	// Mass (1)
	output.push_back(dry_mass - x[SIZE_VECTOR]); // [MASSU]
	
	return output;
}

// Returns the path inequality contraints given
// the current state, the control and the parameters.
// Mass constraints, and thrust constraints.
// For low-thrust cr3bp
// Nineq = 4.
// No unit test.
template<typename T>
DACE::AlgebraicVector<T> inequality_constraints_cr3bp_EARTH_MOON_lt(
	DACE::AlgebraicVector<T> const& x, DACE::AlgebraicVector<T>const& u,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters) {
	// Unpack parameters
	double mu = constants.mu(); // [MU]
	double lu = constants.lu(); // [lu]
	double T_max = spacecraft_parameters.thrust(); // [THRUSTU]
	double dry_mass = spacecraft_parameters.dry_mass(); // [MASSU]
	double initial_mass = spacecraft_parameters.initial_mass(); // [MASSU]
	T x_0(x[0]), x_1(x[1]), x_2(x[2]);

	// Init
	DACE::AlgebraicVector<T> output; output.reserve(1 + 1 + 2);

	// Thrust (1)
	T T_const = u.dot(u) - T_max*T_max; // [THRUSTU²]
	output.push_back(T_const);

	// Mass (1)
	output.push_back(dry_mass - x[SIZE_VECTOR]); // [MASSU]

	// Primaries (2)
	T x_1p2_sqr(DACE::sqr(x_1) + DACE::sqr(x_2));
	T r_1_sqr = DACE::sqr(x_0 - mu) + x_1p2_sqr;
	T r_2_sqr = DACE::sqr(x_0  + 1 - mu) + x_1p2_sqr;
	output.push_back(R_EARTH/lu - r_1_sqr); // Earth
	output.push_back(R_MOON/lu - r_2_sqr); // Moon
	
	return output;
}

// Returns the terminal cost given
// the current state, the target state and the parameters.
// Final difference.
// No unit test.
template<typename T>
T terminal_cost_double_integrator(
	DACE::AlgebraicVector<T> const& x, DACE::vectordb const& x_goal,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters) {
	// Get vectors
	DACE::AlgebraicVector<T> x_(x.extract(0, 2));
	DACE::vectordb x_goal_(x_goal.extract(0, 2));

	// Compute error
	DACE::AlgebraicVector<T> loss = x_ - x_goal_;
	T output = loss.dot(loss);
	return output*1e6;
}


// Returns the terminal cost given
// the current state, the target state and the parameters.
// Final difference.
// No unit test.
template<typename T>
T terminal_cost(
	DACE::AlgebraicVector<T> const& x, DACE::vectordb const& x_goal,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters) {
	// Unpack
	double gain = solver_parameters.terminal_cost_gain(); // [-]
	double homotopy_coefficient = solver_parameters.homotopy_coefficient(); // [-]

	// Get vectors
	DACE::AlgebraicVector<T> x_(x.extract(0, SIZE_VECTOR - 1));
	DACE::vectordb x_goal_(x_goal.extract(0, SIZE_VECTOR - 1));

	// Compute error
	DACE::AlgebraicVector<T> loss = x_ - x_goal_;
	T output = (0.5 * gain) * loss.dot(loss);

	return output;
}

// Returns the terminal cost given
// the current state, the target state and the parameters.
// Final difference, without true longitude.
// No unit test.
template<typename T>
T terminal_cost_equinoctial(
	DACE::AlgebraicVector<T> const& x, DACE::vectordb const& x_goal,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters) {
	// Unpack
	double gain = solver_parameters.terminal_cost_gain(); // [-]

	// Get vectors
	DACE::AlgebraicVector<T> x_(x.extract(0, SIZE_VECTOR - 2));
	DACE::vectordb x_goal_(x_goal.extract(0, SIZE_VECTOR - 2));

	// Compute error
	DACE::AlgebraicVector<T> loss = x_ - x_goal_;
	T output = (0.5 * gain) * loss.dot(loss);

	return output;
}

// Returns the terminal inequality constraints given
// the current state, the target state and the parameters.
// Null.
// No unit test.
template<typename T>
DACE::AlgebraicVector<T> terminal_inequality_constraints_null(
	DACE::AlgebraicVector<T> const& x, DACE::vectordb  const& x_goal,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters) {

	return DACE::AlgebraicVector<T>(0);
}

// Returns the terminal inequality constraints given
// the current state, the target state and the parameters.
// No unit test.
template<typename T>
DACE::AlgebraicVector<T> terminal_inequality_constraints(
	DACE::AlgebraicVector<T> const& x, DACE::vectordb  const& x_goal,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters) {
	// Unpack
	double gain = solver_parameters.terminal_cost_gain(); // [-]
	DACE::matrixdb terminal_cost_inv_covariance = solver_parameters.terminal_cost_inv_covariance().submat(
		SIZE_VECTOR - 1, SIZE_VECTOR - 1); // [-]

	// Get vectors
	DACE::AlgebraicVector<T> x_(x.extract(0, SIZE_VECTOR - 1));
	DACE::vectordb x_goal_(x_goal.extract(0, SIZE_VECTOR - 1));
	DACE::AlgebraicMatrix<T> difference(x_goal_.size(),1);
	difference.setcol(0, x_ - x_goal_);

	// Get Mahalanobis distance
	T mahalanobis_distance = ((difference.transpose()*terminal_cost_inv_covariance)*difference).at(0,0)
		-inv_chi_2_cdf(difference.size(), 1 - 0.05); // 95% sphere 

	// Return 
	return DACE::AlgebraicVector<T>{mahalanobis_distance/gain};
}


// Transformations

// Transforms a keplerian state vector into an equinoctial one.
DACE::vectordb kep_2_equi(DACE::vectordb const& kep_state_vector);

// Transforms a keplerian state vector into an equinoctial one.
DACE::vectordb equi_2_kep(DACE::vectordb const& equi_state_vector);

// Transforms a keplerian state vector into a cartesian one.
DACE::vectordb kep_2_cart(
	DACE::vectordb const& kep_state_vector,
	double const& mu);

// Transforms coordinates in the RTN reference frame into cartesian coordinates.
DACE::vectordb RTN_2_cart(
	DACE::vectordb const& RTN_vector,
	DACE::vectordb const& cart_state_vector);

// Returns dynamics with accelerations.
// Terminal constraints and thrust constraints.
// No unit test.
Dynamics get_double_integrator_dynamics();
Dynamics get_tbp_SUN_lt_dynamics();
Dynamics get_tbp_EARTH_lt_dynamics();
Dynamics get_cr3bp_EARTH_MOON_lt_dynamics();

#endif
