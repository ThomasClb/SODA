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




DACE::vectorDA dth_order_path_inequality_transcription( // PN version
	DACE::vectorDA const& constraints_eval,
	std::vector<DACE::matrixdb> const& list_Sigma, std::vector<DACE::matrixdb> const& list_feedback_gain,
	SpacecraftParameters const& spacecraft_parameters, Constants const& constants,
	SolverParameters const& solver_parameters);

// First order method.
// DOI: WIP

// Path constraints.
DACE::vectorDA first_order_path_inequality_transcription( // DDP version
	DACE::vectorDA const& constraints_eval, stateDA const& x_DA, controlDA const& u_DA,
	SpacecraftParameters const& spacecraft_parameters, Constants const& constants,
	SolverParameters const& solver_parameters);
DACE::vectorDA first_order_path_inequality_transcription( // PN version
	DACE::vectorDA const& constraints_eval,
	DACE::matrixdb const& Sigma_k, DACE::matrixdb const& feedback_gain,
	SpacecraftParameters const& spacecraft_parameters, Constants const& constants,
	SolverParameters const& solver_parameters);
DACE::vectorDA first_order_path_inequality_transcription( // PN version
	DACE::vectorDA const& constraints_eval,
	std::vector<DACE::matrixdb> const& list_Sigma, std::vector<DACE::matrixdb> const& list_feedback_gain,
	SpacecraftParameters const& spacecraft_parameters, Constants const& constants,
	SolverParameters const& solver_parameters);

// Terminal constraints
DACE::vectorDA first_order_terminal_inequality_transcription( // DDP version
	DACE::vectorDA const& constraints_eval, stateDA const& x_DA, statedb const& x_goal,
	SpacecraftParameters const& spacecraft_parameters, Constants const& constants,
	SolverParameters const& solver_parameters);
DACE::vectorDA first_order_terminal_inequality_transcription( // PN version
	DACE::vectorDA const& constraints_eval, DACE::matrixdb const& Sigma_k,
	SpacecraftParameters const& spacecraft_parameters, Constants const& constants,
	SolverParameters const& solver_parameters);


// Spectral radius method.
// Generalisation of [Ozaki et al. 2020]
// DOI: https://doi.org/10.2514/1.G004363

// Path constraints.
DACE::vectorDA spectral_radius_path_inequality_transcription( // DDP version
	DACE::vectorDA const& constraints_eval, stateDA const& x_DA, controlDA const& u_DA,
	SpacecraftParameters const& spacecraft_parameters, Constants const& constants,
	SolverParameters const& solver_parameters);
DACE::vectorDA spectral_radius_path_inequality_transcription( // PN version
	DACE::vectorDA const& constraints_eval,
	DACE::matrixdb const& Sigma_k, DACE::matrixdb const& feedback_gain,
	SpacecraftParameters const& spacecraft_parameters, Constants const& constants,
	SolverParameters const& solver_parameters);

// Terminal constraints
DACE::vectorDA spectral_radius_terminal_inequality_transcription( // DDP version
	DACE::vectorDA const& constraints_eval, stateDA const& x_DA, statedb const& x_goal,
	SpacecraftParameters const& spacecraft_parameters, Constants const& constants,
	SolverParameters const& solver_parameters);
DACE::vectorDA spectral_radius_terminal_inequality_transcription( // PN version
	DACE::vectorDA const& constraints_eval, DACE::matrixdb const& Sigma_k,
	SpacecraftParameters const& spacecraft_parameters, Constants const& constants,
	SolverParameters const& solver_parameters);



#endif
