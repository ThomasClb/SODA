/**
	pn_solver.h

	Purpose: Implementation of the PNSolver class.

	@author Thomas Caleb

	@version 2.0 02/09/2023
*/

#ifndef DEF_PN
#define DEF_PN

#pragma once

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <chrono>

#include <dace/dace_s.h>

#include "aul_solver.h"
#include "stats.h"

// Stores the list of active index
// N concatenated [active_index_ineq, active_index_cont]
// concatenated with [active_index_tineq]
using bi_vector_size_t = std::vector<std::vector<std::size_t>>;

// Stores the linearised constraint
// [constraints, Jacobian, active constraints index]
using linearised_constraints = std::tuple<
	DACE::vectordb, std::vector<DACE::matrixdb>, 
	std::vector<std::vector<std::size_t>>>;

class PNSolver {

	// Attributes
protected:
	AULSolver AULsolver_;
	SolverParameters solver_parameters_;
	SpacecraftParameters spacecraft_parameters_;
	Dynamics dynamics_;
	std::vector<statedb> list_x_; // List of states
	std::vector<controldb> list_u_; // List of controls
	double cost_; // Output cost [-]
	double violation_; // Constraints violation [-]
	double d_th_order_failure_risk_; // Estimation of the failure risk [-]
	std::size_t n_iter_; // Number of iterations

	// Looping attributes
	DACE::vectordb X_U_; // Concatenated states and controls
	DACE::vectordb INEQ_; // Concatenated constraints
	DACE::vectordb correction_; // Vector of all corrections
	std::vector<DACE::matrixdb> der_INEQ_; // Concatenated constraints derivatives
	std::vector<std::vector<DACE::matrixdb>> list_der_cost_; // Output cost derivatives [-]
	std::vector<DACE::vectorDA> list_dynamics_eval_;
	std::vector<DACE::matrixdb> list_feedback_gain_;
	std::vector<DACE::matrixdb> list_Sigma_;
	std::vector<DACE::vectorDA> list_constraints_eval_;

// Methods
public:
	// Empty constructors
	// No unit test.
	PNSolver();

	// Constructors
	// No unit test.
	PNSolver(AULSolver const& AULSolver);

	// Copy constructor
	// No unit test.
	PNSolver(PNSolver const& solver); 

	// Destructors
	// No unit test.
	~PNSolver();

	// Getters
	// No unit test.
	const AULSolver AULsolver() const;
	const DDPSolver DDPsolver() const;
	const std::vector<statedb> list_x() const;
	const std::vector<controldb> list_u() const;
	const double cost() const;
	const double violation() const;
	const double d_th_order_failure_risk() const;
	const std::vector<DACE::vectorDA> list_dynamics_eval() const;
	const std::size_t n_iter() const;

	// Setters
	// No unit test.
	void set_list_x_u();

	// Solves the optimisation problem with a projected Newton method
	// Inspired by ALTRO (Julia).
	// See: https://github.com/RoboticExplorationLab/Altro.jl
	// No unit test.
	void solve(statedb const& x_goal);
	
	// Iterative line search for PN.
	// Inspired by ALTRO (Julia).
	// See: https://github.com/RoboticExplorationLab/Altro.jl
	// No unit test.
	double line_search_(
		statedb const& x_goal,
		sym_tridiag_matrixdb const& tridiag_L,
		std::vector<DACE::matrixdb> const& block_D,
		DACE::vectordb const& d_0,
		double const& violation_0);

	// Evaluates the dth order risk.
	// DOI: WIP.
	// No unit test.
	double evaluate_risk();

	// Computes the maximum continuity constraints.
	// No unit test.
	double get_max_continuity_constraint_(
		DACE::vectordb const& INEQ);

	// Computes the maximum constraints given eq and ineq constraints.
	// No unit test.
	double get_max_constraint_(
		DACE::vectordb const& INEQ);

	// Computes the new constraints given states and controls.
	// Can recompute all DA maps and the dynamics.
	// Updates the derivatives.
	// No unit test.
	void update_constraints_(
		statedb const& x_goal, bool const& force_DA);

	// Computes the new constraints given states and controls.
	// Uses the DA mapping of the dynamics.
	// Faster than update_constraints_.
	// No unit test.
	DACE::vectordb update_constraints_double_(
			statedb const& x_goal,
			DACE::vectordb const& X_U,
			DACE::vectordb const& correction);

	// Returns the vector of active constraints and their gradients
	// - first it the active constraints vector.
	// - second is a pair with the list of gradients of constraints first.
	// - third is the list of active constraints.
	// No unit test.
	linearised_constraints get_linearised_constraints_();

	// Return the matrix S = D_a * D_a^t 
	// Where D_a is the gradient of the active linearized constraints.
	// Using tridiagonal symetric block computation.
	// No unit test.
	sym_tridiag_matrixdb get_block_S_sq_(
		std::vector<DACE::matrixdb> const& block_Delta,
		std::vector<std::vector<std::size_t>> const& list_active_index);
};

#endif
