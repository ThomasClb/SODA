/**
	ddp_solver.h

	Purpose: Implementation of the DDPSolver class.

	@author Thomas Caleb

	@version 2.0 02/09/2024
*/

#ifndef DEF_DDP
#define DEF_DDP

#pragma once

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <chrono>

#include <dace/dace_s.h>

#include "state_t.h"
#include "control_t.h"
#include "parameters.h"
#include "derivation.h"
#include "dynamics.h"
#include "loads.h"
#include "linalg.h"
#include "settings.h"
#include "transcription.h"

class DDPSolver {

// Attributes
protected:
	SolverParameters solver_parameters_; // Solver parameters
	SpacecraftParameters spacecraft_parameters_; // Spacecraft parameters
	Dynamics dynamics_; // Dynamics
	std::vector<statedb> list_x_; // List of states
	std::vector<controldb> list_u_; // List of controls
	std::vector<DACE::vectordb> list_ineq_; // List of inequality constraints
	DACE::vectordb tineq_; // List of terminal inequality constraints
	double cost_; // Output cost [-]

	// Internal variables
	unsigned int n_iter_;
	bool recompute_dynamics_;
	double rho_;
	double alpha_;
	double d_rho_;
	std::vector<DACE::vectorDA> list_dynamic_eval_;
	std::vector<DACE::vectorDA> list_constraints_eval_;
	std::vector<DACE::vectorDA> list_deterministic_constraints_eval_;
	std::vector<DACE::matrixdb> list_Qu_;
	std::vector<DACE::matrixdb> list_k_;
	std::vector<DACE::matrixdb> list_K_;

// Methods
public:
	// Empty constructor
	// No unit test.
	DDPSolver();

	// Constructor
	// No unit test.
	DDPSolver(
		SolverParameters const& solver_parameters,
		SpacecraftParameters const& spacecraft_parameters,
		Dynamics const& dynamics);
	
	// Copy constructor
	// No unit test.
	DDPSolver(DDPSolver const& solver);

	// Destructors
	// No unit test.
	~DDPSolver();

	// Getters
	// No unit test.
	const SolverParameters solver_parameters() const;
	const SpacecraftParameters spacecraft_parameters() const;
	const Dynamics dynamics() const;
	const std::vector<statedb> list_x() const;
	const std::vector<controldb> list_u() const;
	const std::vector<DACE::vectordb> list_ineq() const;
	const std::vector<DACE::vectorDA> list_dynamic_eval() const;
	const std::vector<DACE::vectorDA> list_deterministic_constraints_eval() const;
	const DACE::vectordb tineq() const;
	const double cost() const;
	const unsigned int n_iter() const;

	// Setters
	void set_list_lambda(std::vector<DACE::vectordb> const& list_lambda);
	void set_list_mu(std::vector<DACE::vectordb> const& list_mu);
	void set_ToF(double const& ToF);
	void set_homotopy_coefficient(double const& homotopy_coefficient);
	void set_huber_loss_coefficient(double const& huber_loss_coefficient);
	void set_recompute_dynamics(bool const& recompute_dynamics);
	void set_path_quantile(double const& path_quantile);
	void set_terminal_quantile(double const& terminal_quantile);
	void set_navigation_error_covariance(DACE::matrixdb const& navigation_error_covariance);

	// Returns the Augmented lagrangian cost-to-go:
	// AUL_sc = sc + Lambda^T * c(x) + 0.5 * c(x)^T * I_mu * c(x)
	// No unit test.
	DACE::vectorDA get_AUL_stage_cost(
		stateDA const& x_star_DA, controlDA const& u_star_DA, std::size_t const& index);
		
	// Returns the Augmented lagrangian cost-to-go:
	// AUL_sc = sc + Lambda^T * c(x) + 0.5 * c(x)^T * I_mu * c(x)
	// No unit test.
	DACE::vectorDA get_AUL_terminal_cost(
		stateDA const& x_star, statedb const& x_goal);

	// Increases the regulation to ensure Quu + rho*I
	// is symeteric positive definite.
	// No unit test.
	void increase_regularisation_();

	// Increases the regulation to ensure Quu + rho*I
	// is symeteric positive definite.
	// With a safe guard.
	// No unit test.
	void increase_regularisation_(DACE::matrixdb const& Quu);

	// Decreases the regulation.
	// No unit test.
	void decrease_regularisation_();

	// Computes the expected cost after backward sweep for linesearch.
	// No unit test.
	double expected_cost_(
		double const& alpha);

	// Evaluates the convergence of DDP optimisation.
	// No unit test.
	bool evaluate_convergence_(
		double const& d_cost);

	// Compute the value of the max constraint.
	// No unit test.
	double get_max_constraint_();

	// Returns a completed state with STM, derivatives, and covariance.
	// No unit test.
	statedb make_state(
		unsigned int const& Nx, unsigned int const& Nu,
		DACE::vectorDA const& x_k_DA, statedb const& x_km1, controldb const& u_km1);

	// Performs the DDP backward sweep, that consists in the computation
	// of the gains corrections.
	// Inspired from ALTRO (Julia).
	// See: https://github.com/RoboticExplorationLab/Altro.jl
	// No unit test.
	void backward_sweep_();

	// Performs the DDP forward pass, that consists in the computation
	// of the new states and control after correction.
	// Inspired from ALTRO (Julia).
	// DA only for automatic differentiation.
	// See: https://github.com/RoboticExplorationLab/Altro.jl
	// No unit test.
	void forward_pass_(
		std::vector<statedb> const& list_x,
		std::vector<controldb> const& list_u,
		statedb const& x_goal);

	// Performs the DDP forward pass, that consists in the computation
	// of the new states and control after correction using the DA mapping
	// The linesearch is tweaked to implement a memory from one iteration to the other.
	// No unit test.
	void forward_pass_DA_(
		std::vector<statedb> const& list_x,
		std::vector<controldb> const& list_u,
		statedb const& x_goal);

	// Performs DDP solving given a starting point,
	// initial controls and a final state.
	// Inspired from ALTRO (Julia).
	// See: https://github.com/RoboticExplorationLab/Altro.jl
	// No unit test.
	void solve(
		statedb const& x0,
		std::vector<controldb> const& list_u_init,
		statedb const& x_goal);
};

// Evaluates the convergence radius of a DA vector.
// No unit test.
double convRadius(DACE::vectorDA const& x, double const& tol);

#endif
