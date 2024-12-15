/**
	linalg.cpp

	Purpose: Implementation of the linear algebra methods.

	@author Thomas Caleb

	@version 2.0 12/12/2024
*/

#ifndef DEF_LINALG
#define DEF_LINALG

#pragma once

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>

#include <dace/dace_s.h>

#include "settings.h"
#include "derivation.h"

// Define a symmetric tridagonal matxix type
// The first vector contains the diagonal
// The second vector contains the sub-diagonal
using sym_tridiag_matrixdb = std::pair<
	std::vector<DACE::matrixdb>,
	std::vector<DACE::matrixdb>>;

// Wraps a value between 0 and mod > 0.
double wrap_mod(double const& value, double const& mod);
DACE::DA wrap_mod(DACE::DA const& value, double const& mod); // DA version

// Turns a sym_tridiag_matrixdb to a matrixdb
// upperdiag = true (default) means the sub-diagonal is copied on the upper diagonal
// upperdiag = false means the upper diagonal is 0, therefore, the output is not symmetric
DACE::matrixdb sym_tridiag_matrixdb_2_matrixdb_(
	sym_tridiag_matrixdb const& tridiag,
	bool const& upperdiag = true);

// Gets an identity matrix of size n.
DACE::matrixdb identity_(std::size_t const& n);

// Determines the norm of A - diag(A).
double diagonal_error_(
	DACE::matrixdb const& A);

// Returns a vector equal to the diagonal of A.
template<typename T> DACE::AlgebraicVector<T> get_diag_vector_(
	DACE::AlgebraicMatrix<T> const& A) {
	// Unpack and init
    std::size_t n = A.ncols();
    DACE::AlgebraicVector<T>  v(n);

    // Assign and return
    for (std::size_t i = 0; i < n; i++) { v[i] = A.at(i, i); }
    return v;
}

// Returns a diagonal matrix equal to a given vector.
template<typename T> DACE::AlgebraicMatrix<T> make_diag_matrix_(
	DACE::AlgebraicVector<T> const& diag) {
	// Unpack and init
    std::size_t n = diag.size();
    DACE::AlgebraicMatrix<T>  A(n, n, 0.0);

    // Assign and return
    for (std::size_t i = 0; i < n; i++) { A.at(i,i) = diag[i]; }
    return A;
}

// Turns a matrix to a line.
DACE::vectordb matrix_to_vector(DACE::matrixdb const& M);

// Computes the trace of a matrix
double trace_(DACE::matrixdb const& A);

// Checks if a given symmetric matrix is positive-definite
// Using the eigenvalue criterion: Sp(S) > 0 and S is sym <=> S is sym def pos.
bool is_def_pos_(DACE::matrixdb const& S);

// Computes the Frobenius norm of a matrix
// That is sqrt(trace(A * A^t))
double frobenius_norm_(DACE::matrixdb const& A);

// Inverts a lower triangular matrix.
DACE::matrixdb inv_traingluar_matrix_(
	DACE::matrixdb const& L);

// Computes the Cholesky factorization of
// a symetric positive-definite matrix.
// That is the unique Lower triangular matrix such that:
// S = L * L^t
// See https://en.wikipedia.org/wiki/Cholesky_decomposition
DACE::matrixdb cholesky_(DACE::matrixdb const& S);

// Computes the Cholesky factorization of
// a symetric positive-definite matrix.
// That is the unique Lower triangular matrix such that:
// S = L * L^t
// See https://en.wikipedia.org/wiki/Cholesky_decomposition
// DA version.
DACE::matrixDA cholesky_(DACE::matrixDA const& S);

// Computes the Cholesky factorization of
// a symetric positive-definite tridiagonal matrix.
// That is the unique Lower triangular matrix such that:
// S = L * L^t
// See https://en.wikipedia.org/wiki/Cholesky_decompositionsym_tridiag_matrix 
// This algorithm adapted to triadiagonal matrices was copied from [Cao et al. 2002]
// DOI: https://doi.org/10.1109/ICPPW.2002.1039748
sym_tridiag_matrixdb cholesky_(sym_tridiag_matrixdb const& tridiag_S);

// Perfoms lower triangular system solving.
DACE::vectordb forward_substitution_(
	DACE::matrixdb const& L, DACE::vectordb const& b,
	std::size_t const& starting_index=0);

// Perfoms upper triangular system solving.
// L is a lower triangular matrix, to avoid costly transposition.
DACE::vectordb backward_substitution_(
	DACE::matrixdb const& L, DACE::vectordb const& b);

// Solves the system: L * L^t * x = b
// Where L is lower triangular, and b is a vector
DACE::vectordb solve_cholesky_(
	DACE::matrixdb const& L,
	DACE::vectordb const& b);

// Solves the system: L * L^t * X = B
// Where L is lower triangular, and B is a matrix
DACE::matrixdb solve_cholesky_(
	DACE::matrixdb const& L,
	DACE::matrixdb const& B);

// Solves the system: L * L^t * x = b
// Where L is lower triangular and tridiagonal, and b is a vector
DACE::vectordb solve_cholesky_(
	sym_tridiag_matrixdb const& tridiag_L,
	DACE::vectordb const& b);

// Jacobi eigenvalue algorithm
// Only for symetric matrices
// See https://en.wikipedia.org/wiki/Jacobi_eigenvalue_algorithm
std::pair<DACE::vectordb, DACE::matrixdb> jacobi_eigenvalue_(
	DACE::matrixdb const& S_);

// Jacobi eigenvalue algorithm
// Only for symetric matrices
// See https://en.wikipedia.org/wiki/Jacobi_eigenvalue_algorithm
// DA version
std::pair<DACE::vectorDA, DACE::matrixDA> jacobi_eigenvalue_(
	DACE::matrixDA const& S_);

#endif
