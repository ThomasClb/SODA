/**
	sstate_t.h

	Purpose: Implementation of the State class.

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
#include "state.h"

// Constructors

// Default constructors
template<typename T> State<T>::State() :
	nominal_state_(), Sigma_(), der_dynamics_() {}

// Constructor with size
template<typename T> State<T>::State(size_t const& size) :
	nominal_state_(size), Sigma_(), der_dynamics_()  {}

// Constructor with size and elements value
template<typename T> State<T>::State(size_t const& size, T const& d) : 
	nominal_state_(size, d), Sigma_(), der_dynamics_()  {}

// From list
template<typename T> State<T>::State(std::vector<T> const& v) :
	nominal_state_(v), Sigma_(), der_dynamics_()  {}

// Copy constructor
template<typename T> State<T>::State(State<T> const& v) :
	nominal_state_(v.nominal_state_),
	Sigma_(v.Sigma_), der_dynamics_(v.der_dynamics_)  {}

// Destructors
template<typename T> State<T>::~State() {}


// Getters
template<typename T> const DACE::AlgebraicVector<T> State<T>::nominal_state() const {return nominal_state_;}
template<typename T> const DACE::AlgebraicMatrix<double> State<T>::Sigma() const {return Sigma_;}
template<typename T> const DACE::AlgebraicMatrix<double> State<T>::der_dynamics() const {return der_dynamics_;}


// Setters
template<typename T> void State<T>::set_nominal_state(DACE::AlgebraicVector<T> const& v) {nominal_state_=v;};
template<typename T> void State<T>::set_Sigma(DACE::AlgebraicMatrix<double> const& M) {Sigma_=M;};
template<typename T> void State<T>::set_der_dynamics(DACE::AlgebraicMatrix<double> const& M) {der_dynamics_=M;};

// Operators
template<typename T> State<T>& State<T>::operator=(State<T> const& obj) {
	this->nominal_state_ = obj.nominal_state_;
	this->Sigma_ = obj.Sigma_;
	this->der_dynamics_ = obj.der_dynamics_;
    return *this;
}
template<typename T> T& State<T>::operator[](unsigned int const& index) {
    return this->nominal_state_[index];
}
template<typename T> State<double> const State<T>::cons() const {
	State<double> output(this->nominal_state_.cons());
	output.set_Sigma(this->Sigma_);
	output.set_der_dynamics(this->der_dynamics_);
	return output;
}

// IO operator
template<typename U> std::ostream& operator<<(std::ostream& os, const State<U>& state) {
	os << state.nominal_state_;
	os << state.Sigma_;
	os << state.der_dynamics_;
	return os;
}
template<typename U> std::istream& operator>>(std::istream& is, State<U>& state) {
	is >> state.nominal_state_;
	is >> state.Sigma_;
	is >> state.der_dynamics_;
	return is;
}

// Building shortcut
template<typename T> State<T> make_initial_state(
	DACE::AlgebraicVector<T> const& nominal_state, DACE::vectordb const& covariance_diag, std::size_t const& Nu) {
	std::size_t Nx = nominal_state.size();

	State<T> x0(nominal_state);
	DACE::matrixdb L_sigma(Nx, Nx, 0.0);
	for (std::size_t i=0; i<Nx; i++) {
		L_sigma.at(i,i) = covariance_diag[i];
	}
	x0.set_Sigma(L_sigma*L_sigma.transpose());
	x0.set_der_dynamics(DACE::matrixdb(Nx, Nx+Nu, 0.0));
	return x0;
}