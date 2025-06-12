/**
	transcription.h

	Purpose: Implementation of the transcription methods.

	@author Thomas Caleb

	@version 1.0 09/09/2024
*/

#ifndef DEF_TRANSCRIPTION
#define DEF_TRANSCRIPTION

#pragma once

#include <vector>
#include <cmath>

#include <dace/dace_s.h>

#include "dynamics.h"
#include "state_t.h"
#include "derivation.h"
#include "linalg.h"
#include "stats.h"


// Adaptive first order transcription for Newton method.
// DOI: 10.48550/arXiv.2502.15949
DACE::vectorDA first_order_transcription(
	DACE::vectorDA const& constraints_eval,
	std::vector<DACE::matrixdb> const& list_Sigma, std::vector<DACE::matrixdb> const& list_feedback_gain,
	double const& fact_conservatism,
	SpacecraftParameters const& spacecraft_parameters, Constants const& constants,
	SolverParameters const& solver_parameters);

// First order method for path constraints.
// DOI: 10.48550/arXiv.2502.15949
DACE::vectorDA first_order_path_transcription(
	DACE::vectorDA const& constraints_eval, stateDA const& x_DA, controlDA const& u_DA,
	SpacecraftParameters const& spacecraft_parameters, Constants const& constants,
	SolverParameters const& solver_parameters);

// First order method for terminal constraints.
// DOI: 10.48550/arXiv.2502.15949
DACE::vectorDA first_order_terminal_transcription(
	DACE::vectorDA const& constraints_eval, stateDA const& x_DA, statedb const& x_goal,
	SpacecraftParameters const& spacecraft_parameters, Constants const& constants,
	SolverParameters const& solver_parameters);


// Spectral radius method for path constraints.
// Generalisation of [Ridderhof et al. 2021]
// DOI: 10.48550/arXiv.2502.15949
DACE::vectorDA spectral_radius_path_transcription( 
	DACE::vectorDA const& constraints_eval, stateDA const& x_DA, controlDA const& u_DA,
	SpacecraftParameters const& spacecraft_parameters, Constants const& constants,
	SolverParameters const& solver_parameters);

// Spectral radius method for terminal constraints.
// Generalisation of [Ridderhof et al. 2021]
// DOI: 10.48550/arXiv.2502.15949
DACE::vectorDA spectral_radius_terminal_transcription( 
	DACE::vectorDA const& constraints_eval, stateDA const& x_DA, statedb const& x_goal,
	SpacecraftParameters const& spacecraft_parameters, Constants const& constants,
	SolverParameters const& solver_parameters);


#endif
