/**
	IO.h

	Purpose: Implementation of the inputs and outputs of data.

	@author Thomas Caleb

	@version 2.0 13/12/2024
*/

#ifndef DEF_IO
#define DEF_IO

#pragma once

#include <vector>
#include <string>

#include <dace/dace_s.h>
#include "parameters.h"
#include "dynamics.h"
#include "linalg.h"
#include "state_t.h"
#include "control_t.h"
#include "trajectory_split.h"

// Splits a string into substring.
std::vector<std::string> split(std::string s, std::string delimiter);

// Function to print a dataset at a given name in order to
// produce python visuals.
// No unit test.
void print_dataset(
	std::string const& file_name,
	std::string const& system_name,
	SpacecraftParameters const& spacecraft_parameters,
	std::vector<std::vector<std::string>> const& list_title,
	std::vector<std::vector<DACE::vectordb>> const& list_data);

// Function to propagate a vector without control.
// No unit test.
DACE::matrixdb get_mat_reference_trajectory(
	DACE::vectordb const& x_0,
	Dynamics const& dynamics,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters,
	int const& nb_point);

// Prints a transfer dataset on a standardised format.
// No unit test.
void print_sample_trajectory_dataset(
	std::string const& file_name, std::string const& system_name,
	std::vector<DACE::matrixdb> const& list_mat_state,
	std::vector<DACE::matrixdb> const& list_mat_control,
	std::vector<DACE::matrixdb> const& list_mat_path_constraints,
	std::vector<DACE::matrixdb> const& list_mat_terminal_constraints,
	double const& ToF, bool const& robust_solving,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters);

// Prints a transfer dataset on a standardised format
// with reference orbits.
// No unit test.
void print_robust_trajectory_dataset(
	std::string const& file_name,
	std::string const& system_name,
	std::deque<TrajectorySplit> const& robust_trajectory,
	DACE::vectordb const& x_0, DACE::vectordb const& x_f,
	double const& ToF, bool const& robust_solving,
	Dynamics const& dynamics,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters);

// Loads a printed robust trajectory from a printed file.
// No unit test.
std::pair<
	std::vector<statedb>,
	std::vector<controldb>> load_robust_trajectory(
	std::string const& file_name,
	double const& ToF, bool const& robust_solving,
	Dynamics const& dynamics,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters);

#endif
