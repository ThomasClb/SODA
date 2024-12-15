/**
	control_t.h

	Purpose: Implementation of the Control class.

	@author Thomas Caleb

	@version 1.0 03/05/2024
*/


#pragma once

#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include <dace/dace_s.h>

#include "settings.h"
#include "constants.h"
#include "control.h"

// Constructors

// Default constructors.
template<typename T> Control<T>::Control() :
	nominal_control_(), feedback_gain_() {}

// Constructor with size.
template<typename T> Control<T>::Control(size_t const& size) :
	nominal_control_(size), feedback_gain_() {}

// Constructor with size and elements value.
template<typename T> Control<T>::Control(size_t const& size, T const& d) :
	nominal_control_(size, d), feedback_gain_() {}

// From list.
template<typename T> Control<T>::Control(std::vector<T> const& v) :
	nominal_control_(v), feedback_gain_(){}

// Copy constructor.
template<typename T> Control<T>::Control(Control<T> const& v) :
	nominal_control_(v.nominal_control_), feedback_gain_(v.feedback_gain_) {}

// Destructors.
template<typename T> Control<T>::~Control() {}

// Getters.
template<typename T> const DACE::AlgebraicVector<T> Control<T>::nominal_control() const {return nominal_control_;}
template<typename T> const DACE::AlgebraicMatrix<double> Control<T>::feedback_gain() const {return feedback_gain_;}

// Setters.
template<typename T> void Control<T>::set_nominal_control(DACE::AlgebraicVector<T> const& v) { nominal_control_=v; };
template<typename T> void Control<T>::set_feedback_gain(DACE::AlgebraicMatrix<double> const& M) { feedback_gain_=M; };

// Operators
template<typename T> Control<T>& Control<T>::operator=(Control<T> const& obj) {
	this->nominal_control_ = obj.nominal_control_;
	this->feedback_gain_ = obj.feedback_gain_;
    return *this;
}
template<class T> T& Control<T>::operator[](unsigned int const& index) {
    return this->nominal_control_[index];
}
template<class T> Control<double> const Control<T>::cons() const {
	Control<double> output(this->nominal_control_.cons());
	output.set_feedback_gain(this->feedback_gain_);
	return output;
}

// IO operator
template<typename U> std::ostream& operator<<(std::ostream& os, const Control<U>& control) {
	os << control.nominal_control_;
	os << control.feedback_gain_;
	return os;
}
template<typename U> std::istream& operator>>(std::istream& is, Control<U>& control) {
	is >> control.nominal_control_;
	is >> control.feedback_gain_;
	return is;
}

// Building shortcut
template<typename T> Control<T> make_first_guess(
	DACE::AlgebraicVector<T> const& nominal_control, double const& initial_gain_value, std::size_t const& Nx) {
	std::size_t Nu = nominal_control.size();

	Control<T> u_init(nominal_control);
	DACE::matrixdb feedback_gain(Nu, Nx, initial_gain_value);
	u_init.set_feedback_gain(feedback_gain);
	return u_init;
}
