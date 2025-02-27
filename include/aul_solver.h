/**
	aul_solver.h

	Purpose: Implementation of the AULSolver class.

	@author Thomas Caleb

	@version 2.0 02/09/2024
*/

#ifndef DEF_AUL
#define DEF_AUL

#pragma once

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <chrono>

#include <dace/dace_s.h>

#include "ddp_solver.h"
#include "trajectory_split.h"

class AULSolver {

// Attributes
protected:
	DDPSolver DDPsolver_; // Solver parameters
	TrajectorySplit trajectory_split_; // Compiles a list of states, controls and a splitting history
	double cost_; // Output cost [-]
	double violation_; // Constraints violation [-]
	double d_th_order_failure_risk_; // Estimation of the failure risk [-]
	std::vector<DACE::vectordb> list_ineq_; // List of inequality constraints
	DACE::vectordb tineq_; // List of terminal inequality constraints
	std::vector<DACE::vectordb> list_lambda_; // List of Lagrange multiplicator
	std::vector<DACE::vectordb> list_mu_; // List of penalty factors
	std::deque<std::pair<
		std::vector<DACE::vectordb>,std::vector<DACE::vectordb>>> list_lambda_mu_; // List of penalty factors and Lagrange multiplicator


	// Iterations
	unsigned int AUL_n_iter_; // Number of AUL iterations
	unsigned int DDP_n_iter_; // Total number of DDP iterations

// Methods
public:
	// Empty constructor
	// No unit test.
	AULSolver();

	// Constructor
	// No unit test.
	AULSolver(
		SolverParameters const& solver_parameters,
		SpacecraftParameters const& spacecraft_parameters,
		Dynamics const& dynamics);

	// Copy constructor
	// No unit test.
	AULSolver(AULSolver const& solver);

	// Destructors
	// No unit test.
	~AULSolver();

	// Getters
	// No unit test.
	const DDPSolver DDPsolver() const;
	const TrajectorySplit trajectory_split() const;
	const double cost() const;
	const double violation() const;
	const double d_th_order_failure_risk() const;
	const std::vector<DACE::vectordb> list_ineq() const;
	const DACE::vectordb tineq() const;
	const std::vector<DACE::vectordb> list_lambda() const;
	const std::vector<DACE::vectordb> list_mu() const;
	const unsigned int AUL_n_iter() const;
	const unsigned int DDP_n_iter() const;

	// Setters
	void set_ToF(double const& ToF);
	void set_homotopy_coefficient(double const& homotopy_coefficient);
	void set_huber_loss_coefficient(double const& huber_loss_coefficient);
	void set_path_quantile(double const& path_quantile);
	void set_terminal_quantile(double const& terminal_quantile);
	void set_navigation_error_covariance(DACE::matrixdb const& navigation_error_covariance);

	// Update dual state in Augmented Lagrangian formulation.
	// No unit test.
	void update_lambda_();

	// Update penalty in Augmented Lagrangian formulation.
	// No unit test.
	void update_mu_();

	// Computes the d-th order risk.
	// DOI: WIP
	// No unit test.
	double evaluate_risk();

	// Performs AUL solving given a starting point,
	// initial controls and a final state.
	// Inspired from ALTRO (Julia).
	// See: https://github.com/RoboticExplorationLab/Altro.jl
	// No unit test.
	void solve(
		std::deque<TrajectorySplit>* const& p_list_trajectory_split,
		statedb const& x_goal);
};

#endif
