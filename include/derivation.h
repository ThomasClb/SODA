/**
	derivation.h

	Purpose: Implementation of the DA-based
	automatic differentiation methods.

	@author Thomas Caleb

	@version 2.0 02/09/2024
*/

#ifndef DEF_DERIVATION
#define DEF_DERIVATION

#pragma once

#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include <dace/dace_s.h>

#include "state_t.h"
#include "control_t.h"

double cons(double const& value);
double abs_cons(double const& value);
double abs_cons(DACE::DA const& value);

// DACE vector conversion

// Turns a vectordb into a vectorDA without perturbation.
DACE::vectorDA db_2_DA(
	DACE::vectordb const& vector_db);

// Turns a vectordb into a vectorDA with identity perturbation.
DACE::vectorDA id_vector(
	DACE::vectordb const& vector_db,
	unsigned int const& begin_index,
	unsigned int const& begin_da_index,
	unsigned int const& end_da_index);
stateDA id_vector(
	statedb const& vector_db,
	unsigned int const& begin_index,
	unsigned int const& begin_da_index,
	unsigned int const& end_da_index);
controlDA id_vector(
	controldb const& vector_db,
	unsigned int const& begin_index,
	unsigned int const& begin_da_index,
	unsigned int const& end_da_index);

// Differentiation functions

// Differentiates with respect to x, and u.
std::vector<DACE::matrixdb> deriv_xu(
	DACE::vectorDA const& f_eval,
	unsigned int const& Nx, unsigned int const& Nu,
	bool const& hessian_computation);

// Differentiates with respect to x, and u. Returns DA matrices.
std::vector<DACE::matrixDA> deriv_xu_DA(
	DACE::vectorDA const& f_eval,
	unsigned int const& Nx, unsigned int const& Nu,
	bool const& hessian_computation);

// Differentiates with respect to x.
std::vector<DACE::matrixdb> deriv_x(
	DACE::vectorDA const& f_eval,
	unsigned int const& Nx,
	bool const& hessian_computation);

// Differentiates with respect to x. Returns DA matrices.
std::vector<DACE::matrixDA> deriv_x_DA(
	DACE::vectorDA const& f_eval,
	unsigned int const& Nx,
	bool const& hessian_computation);

#endif
