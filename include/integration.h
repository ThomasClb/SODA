/**
	integration.h

	Purpose: Implementation of the integration methods.

	@author Thomas Caleb

	@version 2.0 02/09/2024
*/

#ifndef DEF_INTEGRATION
#define DEF_INTEGRATION

#pragma once

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>

#include <dace/dace_s.h>

#include "settings.h"
#include "state.h"
#include "control.h"
#include "constants.h"
#include "parameters.h"

// Declare integration constants
const std::vector<double> get_A();
const std::vector<double> get_C();
const std::vector<double> get_D();
const DACE::matrixdb get_B();

// Instantiate integration constants
const std::vector<double> A = get_A();
const DACE::matrixdb B = get_B();
const std::vector<double> C = get_C();
const std::vector<double> D = get_D();

// Utilities

// Computes the minimum
double min(double const& a, double const& b);

// Computes the maximum
double max(double const& a, double const& b);

// Computes the infinite norm
double normtmp(int const& N, DACE::vectordb const& X);

// Integrators

// Propagates a state using a order 7(8) Runge - Kutta - Fehlberg embedded method(RK78)
// See [Fehlberg 1968] Classical fifth-, sixth-, seventh-, and eighthorder
// Runge - Kutta formulas with stepsize control.
// The integrated function takes a DA state, a DA control, a time and a SpacecraftParameters
// as inputs and provides a DA derivative of the same size.
// Adapted from: https://github.com/dacelib/dace/blob/master/Tutorials/Tutorial1/Example11DARK78.cpp
DACE::vectorDA RK78(
	DACE::vectorDA(*f)(
		DACE::vectorDA const&, DACE::vectorDA const&,
		double const&,
		SpacecraftParameters const&, Constants const&,
		SolverParameters const&),
	DACE::vectorDA Y0, DACE::vectorDA const& U,
	double const& X0, double const& DX,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters);

// Propagates a state using a order 7(8) Runge - Kutta - Fehlberg embedded method(RK78)
// See [Fehlberg 1968] Classical fifth-, sixth-, seventh-, and eighthorder
// Runge - Kutta formulas with stepsize control.
// The integrated function takes a double state, a double control, a time and a SpacecraftParameters
// as inputs and provides a double derivative of the same size.
// Adapted from: https://github.com/dacelib/dace/blob/master/Tutorials/Tutorial1/Example11DARK78.cpp
DACE::vectordb RK78(
	DACE::vectordb(*f)(
		DACE::vectordb  const&, DACE::vectordb  const&,
		double const&,
		SpacecraftParameters const&, Constants const&,
		SolverParameters const&),
	DACE::vectordb Y0, DACE::vectordb const& U,
	double const& X0, double const& DX,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters);


// Propagates a state using Runge-Kutta 4th order method
// See: https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods
// The integrated function takes a state, a control, a time and a SpacecraftParameters
// as inputs and provides a derivative of the same size.
template<typename T>
DACE::AlgebraicVector<T>  RK4(
	DACE::AlgebraicVector<T>(*acc)(
		DACE::AlgebraicVector<T> const&,
		DACE::AlgebraicVector<T>const&,
		double const&,
		SpacecraftParameters const&, Constants const&,
		SolverParameters const&),
	DACE::AlgebraicVector<T>  const& x,
	DACE::AlgebraicVector<T> const& u,
	double const& t, double const& dt,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters) {

	DACE::AlgebraicVector<T>  k1 = acc(
		x, u, t,
		spacecraft_parameters, constants, solver_parameters);
	DACE::AlgebraicVector<T>  k2 = acc(
		x + (dt / 2) * k1, u, t + dt / 2,
		spacecraft_parameters, constants, solver_parameters);
	DACE::AlgebraicVector<T>  k3 = acc(
		x + (dt / 2) * k2, u, t + dt / 2,
		spacecraft_parameters, constants, solver_parameters);
	DACE::AlgebraicVector<T>  k4 = acc(
		x + dt * k3, u, t + dt,
		spacecraft_parameters, constants, solver_parameters);

	return x + (dt / 6) * (k1 + 2 * (k2 + k3) + k4);
}

#endif
