/**
	state.h

	Purpose: Implementation of the State class.

	@author Thomas Caleb

	@version 1.0 03/05/2024
*/

#ifndef DEF_STATE
#define DEF_STATE

#pragma once

#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include <dace/dace_s.h>

#include "settings.h"
#include "constants.h"


/*

	STATE

*/
template<typename T> class State {

// Attributes
protected:
	DACE::AlgebraicVector<T> nominal_state_;
	DACE::AlgebraicMatrix<double> Sigma_;
	DACE::AlgebraicMatrix<double> der_dynamics_;

// Methods
public:

	// Constructors
	State();																	// Default constructor
    State(size_t const& size);                                              	// Constructor with size
    State(size_t const& size, T const& d);                                      // Constructor with size and elements value
    State(std::vector<T> const& v);                                             // From list    
    State(State<T> const& v);                                                   // Copy constructor

	// Destructors
	~State();

	// Getters
	const DACE::AlgebraicVector<T> nominal_state() const;
	const DACE::AlgebraicMatrix<double> Sigma() const;
	const DACE::AlgebraicMatrix<double> der_dynamics() const;

	// Setters
	void set_nominal_state(DACE::AlgebraicVector<T> const& v);
	void set_Sigma(DACE::AlgebraicMatrix<double> const& M);
	void set_der_dynamics(DACE::AlgebraicMatrix<double> const& M);

	// Operators
	State& operator=(State const& obj);
	T& operator[](unsigned int const& index);
	State<double> const cons() const;

	// IO operator
	template<typename U> friend std::ostream& operator<<(std::ostream& os, const State<U>& state);
	template<typename U> friend std::istream& operator>>(std::istream& is, State<U>& state);

	// Building shortcut
	State<T> make_initial_state(
		DACE::AlgebraicVector<T> const& nominal_state,
		DACE::vectordb const& covariance_diag,
		std::size_t const& Nu);
};

// shortcuts for common state types
typedef State<DACE::DA> stateDA;       //!< Shorthand notation for State<DA>.
typedef State<double> statedb;   //!< Shorthand notation for State<double>.


#endif
