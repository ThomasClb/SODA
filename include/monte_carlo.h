/**
	monte_carlo.h

	Purpose: Implementation of the Monte-Carlo validation methods.

	@author Thomas Caleb

	@version 1.0 26/09/2024
*/

#ifndef DEF_MONTE_CARLO
#define DEF_MONTE_CARLO

#pragma once

#include <vector>
#include <random>

#include <dace/dace_s.h>

#include "dynamics.h"
#include "linalg.h"
#include "soda.h"
#include "robust_trajectory.h"

// Returns a sample of a multivariate normal distribution.
DACE::matrixdb generate_normal_sample(
	std::size_t const& size_sample,
	DACE::vectordb const& mean, DACE::matrixdb const& covariance);

// Returns the empirical quantile value for a given probability.
// No unit test.
double get_quantile(DACE::vectordb const& sample, double const& probability);

// Propagate a trajectory and returns
// the remaining fuel and the max constraint.
// No unit test.
DACE::vectordb propagate_trajectory(
	DACE::vectordb const& x_0_sample,
	statedb const& x_goal,
	RobustTrajectory const& robust_trajectory,
	DACE::matrixdb const& navigation_error_sample,
	SolverParameters const& solver_parameters,
	SpacecraftParameters const& spacecraft_parameters,
	Dynamics const& dynamics, int const& test_case_id,
	DACE::matrixdb* const& p_mat_state,
	DACE::matrixdb* const& p_mat_control,
	DACE::matrixdb* const& p_mat_path_constraints,
	DACE::matrixdb* const& p_mat_terminal_constraints,
	bool const& save_trajectory);

// Tests a robust command law with MC.
// No unit test.
std::vector<std::vector<DACE::matrixdb>> test_trajectory(
	RobustTrajectory const& robust_trajectory,
	statedb const& x_0,
	statedb const& x_goal, std::size_t const& size_sample,
	SODA const& solver, int const& test_case_id,
	bool const& robust_optimisation,
	SolverParameters const& solver_parameters,
	SpacecraftParameters const& spacecraft_parameters,
	Dynamics const& dynamics, size_t const& size_sample_saved_);


#endif
