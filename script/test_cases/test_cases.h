/**
	test_cases.h

	Purpose: Implementation of the test cases for SODA.

	@author Thomas Caleb

	@version 1.0 07/12/2023
*/

#ifndef DEF_TEST_CASES
#define DEF_TEST_CASES

#pragma once

#include <dace/dace_s.h>
#include <chrono>

#include "soda.h"
#include "monte_carlo.h"

// Test cases
void double_integrator(int argc, char** argv);
void tbp_SUN_lt_earth_to_mars(int argc, char** argv);
void cr3bp_EARTH_MOON_lt_haloL2_to_haloL1(int argc, char** argv);
void cr3bp_EARTH_MOON_lt_nrho_to_dro(int argc, char** argv);
void cr3bp_EARTH_MOON_lt_dro_to_dro(int argc, char** argv);
void cr3bp_EARTH_MOON_lt_lyapunovL1_to_lyapunovL2(int argc, char** argv);

// Misc
void save_control(
	std::string const& file_name,
	std::vector<DACE::vectordb> const& list_u);
std::vector<DACE::vectordb> load_control(std::string const& file_name);

// Runs the selected test case.
void run_test_cases(int argc, char** argv);

#endif
