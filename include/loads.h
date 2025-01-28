/**
	loads.h

	Purpose: Implementation of the Low-Order Automatic Domain Splitting methods.

	@author Thomas Caleb

	@version 1.0 23/10/2024
*/

#ifndef DEF_LOADS
#define DEF_LOADS

#pragma once

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>

#include <dace/dace_s.h>

#include "settings.h"
#include "linalg.h"
#include "state_t.h"
#include "control_t.h"
#include "stats.h"
#include "splitting_history.h"
#include "trajectory_split.h"

// Shortcut function to retrieve eigenvalues knowing robust trajectory.
std::pair<DACE::vectordb, DACE::matrixdb> get_eigenvalues(
	DACE::matrixdb const& Sigma_x, DACE::matrixdb const& feedback_gain);

// Scales a vector based on a robust trajectory.
std::pair<DACE::vectorDA, DACE::vectordb> scale(
	DACE::vectorDA const& y, DACE::matrixdb const& Sigma_x,
	DACE::matrixdb const& feedback_gain, double const& transcription_beta);

// Return the GMM decomposition of a given multivariate Gaussian
// Distribution along a specified direction.
// From [DeMars et al. 2013]
// DOI: https://doi.org/10.2514/1.58987
std::vector<std::tuple<double, DACE::vectordb, DACE::matrixdb>> split_gmm(
    DACE::vectordb const& y_mean,
    DACE::matrixdb const& Sigma,
    std::size_t const& direction);

// Return the merged GMM decomposition of a given central part of a GMM.
// Distribution along a specified direction.
// From [Losacco et al. 2024]
// DOI: https://doi.org/10.2514/1.G007271
std::pair<DACE::vectordb, DACE::matrixdb> merge_gmm(
    DACE::vectordb const& y_mean,
    DACE::matrixdb const& Sigma,
    std::size_t const& direction);

// Computes the NonLinearity Index (NLI) of a vector given a scaling.
// From [Losacco et al. 2024]
// DOI: https://doi.org/10.2514/1.G007271
double nl_index(DACE::vectorDA const& y, DACE::vectordb const& lambda);

// Computes the NonLinearity Index (NLI) of a vector given a robust trajectory step.
// From [Losacco et al. 2024]
// DOI: https://doi.org/10.2514/1.G007271
double nl_index(DACE::vectorDA const& y, DACE::matrixdb const& Sigma,
	DACE::matrixdb const& feedback_gain, double const& transcription_beta);

// Computes the NonLinearity Index (NLI) of a vector given a robust trajectory split.
// From [Losacco et al. 2024]
// DOI: https://doi.org/10.2514/1.G007271
DACE::vectordb nl_index(
    std::vector<DACE::vectorDA> const& list_dynamics_eval, std::vector<statedb> const& list_x,
    std::vector<controldb> const& list_u, double const& transcription_beta);

// Computes the NonLinearity Index (NLI) of a vector given a scaling along each direction.
// From [Losacco et al. 2024]
// DOI: https://doi.org/10.2514/1.G007271
DACE::vectordb nl_index_dir(DACE::vectorDA const& y, DACE::vectordb const& lambda);

// Computes the NonLinearity Index (NLI) along each direction
// From t=0 to a given t_index given a robust trajectory split.
// From [Losacco et al. 2024]
// DOI: https://doi.org/10.2514/1.G007271 
DACE::vectordb nl_index_dir(
	std::size_t const& index,
    std::vector<DACE::vectorDA> const& list_dynamics_eval, std::vector<statedb> const& list_x,
    std::vector<controldb> const& list_u, double const& transcription_beta);


#endif
