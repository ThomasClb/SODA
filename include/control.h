/**
	control.h

	Purpose: Implementation of the Control class.

	@author Thomas Caleb

	@version 1.0 03/05/2024
*/

#ifndef DEF_CONTROL
#define DEF_CONTROL

#pragma once

#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include <dace/dace_s.h>

#include "settings.h"
#include "constants.h"


/*

	CONTROL

*/
template<typename T> class Control {

// Attributes
protected:
	DACE::AlgebraicVector<T> nominal_control_;
	DACE::AlgebraicMatrix<double> feedback_gain_;	

// Methods
public:

	// Constructors
	Control();																	  // Default constructor
    Control(size_t const& size);                                              	  // Constructor with size
    Control(size_t const& size, T const& d);                                      // Constructor with size and elements value
    Control(std::vector<T> const& v);                                             // From list 
    Control(Control<T> const& v);                                                 // Copy constructor

	// Destructors
	~Control();

	// Getters
	const DACE::AlgebraicVector<T> nominal_control() const;
	const DACE::AlgebraicMatrix<double> feedback_gain() const;

	// Setters
	void set_nominal_control(DACE::AlgebraicVector<T> const& v);
	void set_feedback_gain(DACE::AlgebraicMatrix<double> const& M);

	// Operators
	Control& operator=(Control const& obj);
	T& operator[](unsigned int const& index);
	Control<double> const cons() const;

	// IO operator
	template<typename U> friend std::ostream& operator<<(std::ostream& os, const Control<U>& control);
	template<typename U> friend std::istream& operator>>(std::istream& is, Control<U>& control);

	// Building shortcut
	Control<T> make_first_guess(
		DACE::AlgebraicVector<T> const& nominal_control, double const& initial_gain_value, std::size_t const& Nx);
};

// shortcuts for common control types
typedef Control<DACE::DA> controlDA;       //!< Shorthand notation for Control<DA>.
typedef Control<double> controldb;   //!< Shorthand notation for Control<double>.


#endif
