/**
	SODA.h

	Purpose: Implementation of the SODA solver class.

	@author Thomas Caleb

	@version 1.0 02/09/2024
*/

#ifndef DEF_SODA
#define DEF_SODA

#pragma once

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <chrono>

#include <dace/dace_s.h>

#include "settings.h"
#include "state.h"
#include "control.h"
#include "state_t.h"
#include "control_t.h"
#include "constants.h"
#include "dynamics.h"
#include "IO.h"
#include "aul_solver.h"
#include "pn_solver.h"
#include "loads.h"


class SODA {

// Attributes
protected:
	AULSolver AULsolver_; // Augmented Lagrangian solver
	SolverParameters solver_parameters_;
	SpacecraftParameters spacecraft_parameters_;
	Dynamics dynamics_;
	PNSolver PNsolver_; // Projected Newton solver
	std::deque<TrajectorySplit> list_trajectory_split_;
	std::vector<statedb> list_x_;
	std::vector<controldb> list_u_;

	// Iterations
	double PN_runtime_; // Runtime of PN
	double AUL_runtime_; // Runtime of AUL
	double runtime_; // Runtime
	std::size_t PN_n_iter_; // Number of PN iterations
	std::size_t AUL_n_iter_; // Number of AUL iterations
	std::size_t DDP_n_iter_; // Total number of DDP iterations

// Methods
public:
	// Empty constructor
	// No unit test.
	SODA();

	// Constructor
	// No unit test.
	SODA(
		SolverParameters const& solver_parameters,
		SpacecraftParameters const& spacecraft_parameters,
		Dynamics const& dynamics);

	// Copy constructor
	// No unit test.
	SODA(SODA const& solver);

	// Destructors
	// No unit test.
	~SODA();

	// Getters
	// No unit test.
	const AULSolver AULsolver() const;
	const PNSolver PNsolver() const;
	const std::vector<statedb> list_x() const;
	const std::vector<controldb> list_u() const;
	const std::deque<TrajectorySplit> list_trajectory_split() const;
	const double cost() const;
	const double violation() const;
	const double d_th_order_failure_risk() const;
	const double PN_runtime() const;
	const double AUL_runtime() const;
	const double runtime() const;
	const std::size_t PN_n_iter() const;
	const std::size_t DDP_n_iter() const;
	const std::size_t AUL_n_iter() const;

	// Performs solving given a starting point,
	// initial controls and a final state.
	// Inspired from ALTRO (Julia).
	// See: https://github.com/RoboticExplorationLab/Altro.jl
	// No unit test.
	void solve(
		statedb const& x0,
		std::vector<controldb> const& list_u_init,
		statedb const& x_goal,
		bool const& robust_solving,
		bool const& fuel_optimal,
		bool const& pn_solving);
};

#endif
