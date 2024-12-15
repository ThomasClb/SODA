/**
	test_linalg.cpp

	Purpose: Test of the implementation of the linear algebra
    methods.
    
	@author Thomas Caleb

	@version 2.0 12/12/2024
*/

#include "gtest/gtest.h"
#include "../include/soda.h"

using namespace DACE;
using namespace std;

// Basic functions
TEST(TestLinalg, WrapMod) {

	// Init DACE
	size_t nb_variables = 3;
	DA::init(2, nb_variables);

	// Init
	double mod = 2*PI;
	double value = 1.0 - 5*mod;
	double value_wrapped = 1.0;
	DA value_DA = value + DA(1);
	DA value_wrapped_DA = value_wrapped + DA(1);

	// Test
	EXPECT_EQ(wrap_mod(value_DA, mod).cons(), value_wrapped_DA.cons());
	EXPECT_EQ(wrap_mod(value, mod), value_wrapped);
}
TEST(TestLinalg, SymTridiagMatrixdb2Matrixdb) {
	// Init
	size_t n = 3;
	matrixdb diag(n, n, 0.0);
	for (size_t i = 0; i < n; i++) {
		diag.at(i, i) = 1.0;
	}
	for (size_t i = 0; i < n - 1; i++) {
		diag.at(i, i + 1) = 0.1;
		diag.at(i + 1, i) = 0.1;
	}
	vector<matrixdb> list_D{ diag, diag };
	vector<matrixdb> list_E{ 0.001 * matrixdb(n, n, 1.0) };
	sym_tridiag_matrixdb tridiag_S(list_D, list_E);
	matrixdb S = sym_tridiag_matrixdb_2_matrixdb_(tridiag_S);

	// Test
	// D0
	for (size_t i = 0; i < n; i++) {
		EXPECT_EQ(S.at(i, i), 1.0);
	}
	for (size_t i = 0; i < n - 1; i++) {
		EXPECT_EQ(S.at(i + 1, i), 0.1);
		EXPECT_EQ(S.at(i, i + 1), 0.1);
	}
	// D1
	for (size_t i = 0; i < n; i++) {
		EXPECT_EQ(S.at(i + n, i + n), 1.0);
	}
	for (size_t i = 0; i < n - 1; i++) {
		EXPECT_EQ(S.at(i + 1 + n, i + n), 0.1);
		EXPECT_EQ(S.at(i + n, i + 1 + n), 0.1);
	}
	// E0
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			EXPECT_EQ(S.at(i + n, j), 0.001);
			EXPECT_EQ(S.at(i, j + n), 0.001);
		}
	}

	// Upperdiag = false
	S = sym_tridiag_matrixdb_2_matrixdb_(tridiag_S, false);

	// Test
	// D0
	for (size_t i = 0; i < n; i++) {
		EXPECT_EQ(S.at(i, i), 1.0);
	}
	for (size_t i = 0; i < n - 1; i++) {
		EXPECT_EQ(S.at(i + 1, i), 0.1);
	}
	// D1
	for (size_t i = 0; i < n; i++) {
		EXPECT_EQ(S.at(i + n, i + n), 1.0);
	}
	for (size_t i = 0; i < n - 1; i++) {
		EXPECT_EQ(S.at(i + 1 + n, i + n), 0.1);
	}
	// E0
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			EXPECT_EQ(S.at(i + n, j), 0.001);
		}
	}
}
TEST(TestLinalg, Indentity) {
	// Init
	size_t n = 10;
	matrixdb A = identity_(n);

	// Test
	for (size_t i = 0; i < n; i++) {
		EXPECT_EQ(A.at(i, i), 1);
	}
}
TEST(TestLinalg, DiagonalError) {
	// Init
	size_t n = 10;
	matrixdb A = identity_(n);
	matrixdb B(n, n, 1.0);

	// Test
	EXPECT_EQ(diagonal_error_(A), 0.0);
	EXPECT_EQ(diagonal_error_(B), sqrt(n*(n - 1.0)));
}
TEST(TestLinalg, GetDiagVector) {
	// Init
	size_t n = 10;
	matrixdb A = identity_(n);
	vectordb diag_A = get_diag_vector_(A);

	// Test
	for (size_t i = 0; i < n; i++) {
		EXPECT_EQ(A.at(i, i), diag_A[i]);
	}
}
TEST(TestLinalg, MakeDiagMatrix) {
	// Init
	size_t n = 10;
	vectordb diag(n, 1.0);
	matrixdb A = make_diag_matrix_(diag);

	// Test
	for (size_t i = 0; i < n; i++) {
		EXPECT_EQ(A.at(i, i), diag[i]);
	}
}
TEST(TestLinalg, MatrixtoVector) {
	// Init
	size_t n = 10;
	size_t m = 12;
	matrixdb A(n, m);
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < m; j++) {
			A.at(i, j) = (1.0*n)*(1.0*m);
		}
	}
	vectordb A_vector = matrix_to_vector(A);

	// Test
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < m; j++) {
			EXPECT_EQ(A.at(i, j), A_vector[i*m + j]);
		}
	}
}
TEST(TestLinalg, Trace) {
	// Init
	size_t n = 10;
	matrixdb A(n, n, 1.0);

	// Test
	EXPECT_EQ(trace_(A), 10);
}
TEST(TestLinalg, IsDefPos) {
	// Init
	size_t n = 10;
	matrixdb S_pos(n, n, 0.0);
	for (size_t i = 0; i < n; i++) {
		S_pos.at(i, i) = 1.0;
	}
	matrixdb S_neg = S_pos; S_neg.at(3, 3) = -100;
	 
	// Test
	EXPECT_TRUE(is_def_pos_(S_pos));
	EXPECT_TRUE(!is_def_pos_(S_neg));
}
TEST(TestLinalg, FrobeniusNorm) {
	// Init
	size_t n = 10;
	size_t m = 5;
	matrixdb A(n, m, 1.0);
	double expected_norm = sqrt(1.0 * m * n);

	// Test
	EXPECT_EQ(frobenius_norm_(A), expected_norm);
}
TEST(TestLinalg, InvTriangularMatrix) {
	// Init
	size_t n = 10;
	matrixdb L(n, n, 0.0);
	for (size_t i = 0; i < n; i++) {
		L.at(i, i) = 1.0;
	}
	for (size_t i = 0; i < n - 1; i++) {
		L.at(i + 1, i) = 0.1;
	}
	matrixdb inv_L = inv_traingluar_matrix_(L);

	// Test
	EXPECT_NEAR(frobenius_norm_(L * inv_L - identity_(n)), 0.0, EPS);
}
TEST(TestLinalg, Cholesky) {
	// Init DACE
	size_t nb_variables = 3;
	DA::init(2, nb_variables);

	// Init
	size_t n = 10;
	matrixdb S(n, n, 0.0);
	for (size_t i = 0; i < n; i++) {
		S.at(i, i) = 1.0;
	}
	for (size_t i = 0; i < n - 1; i++) {
		S.at(i, i + 1) = 0.1;
		S.at(i + 1, i) = 0.1;
	}
	matrixdb L = cholesky_(S);
	matrixDA S_DA(S + DA(1));
	matrixDA L_DA = cholesky_(S_DA);

	// Test
	EXPECT_NEAR(frobenius_norm_(L * L.transpose() - S), 0.0, EPS);
	EXPECT_NEAR(frobenius_norm_((L_DA * L_DA.transpose() - S_DA).cons()), 0.0, EPS);
}
TEST(TestLinalg, SymTridiagMatrixdbCholesky) {
	// Init
	size_t n = 10;
	matrixdb diag(n, n, 0.0);
	for (size_t i = 0; i < n; i++) {
		diag.at(i, i) = 1.0;
	}
	for (size_t i = 0; i < n - 1; i++) {
		diag.at(i, i + 1) = 0.1;
		diag.at(i + 1, i) = 0.1;
	}
	vector<matrixdb> list_D{ diag, diag };
	vector<matrixdb> list_E{ 0.001* matrixdb(n, n, 1.0) };
	sym_tridiag_matrixdb tridiag_S(list_D, list_E);
	matrixdb S = sym_tridiag_matrixdb_2_matrixdb_(tridiag_S);
	sym_tridiag_matrixdb tridiag_L = cholesky_(tridiag_S);
	matrixdb L = sym_tridiag_matrixdb_2_matrixdb_(tridiag_L, false);

	// Test
	EXPECT_NEAR(frobenius_norm_(L * L.transpose() - S), 0.0, EPS);
}
TEST(TestLinalg, ForwardSubstitution) {
	// Init
	size_t n = 10;
	matrixdb S(n, n, 0.0);
	vectordb b(n);
	for (size_t i = 0; i < n; i++) {
		S.at(i, i) = 1.0;
		b[i] = 10.0 * (i + 1);
	}
	for (size_t i = 0; i < n - 1; i++) {
		S.at(i, i + 1) = 0.1;
		S.at(i + 1, i) = 0.1;
	}
	matrixdb L = cholesky_(S);
	vectordb x = forward_substitution_(L, b);

	// Test
	EXPECT_NEAR((L * x - b).vnorm(), 0.0, EPS);

	// With strating index
	size_t strating_index = 4;
	vectordb b_red = b.extract(strating_index, n - 1);
	matrixdb L_red = L.submat(
		strating_index, strating_index, n - 1, n - 1);
	vectordb x_red_expected = forward_substitution_(L_red, b_red, 0);
	vectordb x_red = forward_substitution_(L, b, strating_index);

	// Test
	EXPECT_NEAR((x_red_expected - x_red.extract(strating_index, n - 1)).vnorm(), 0.0, EPS);
	EXPECT_EQ(x_red.extract(0, strating_index- 1).vnorm(), 0.0);
}
TEST(TestLinalg, BackwardSubstitution) {
	// Init
	size_t n = 10;
	matrixdb S(n, n, 0.0);
	vectordb b(n);
	for (size_t i = 0; i < n; i++) {
		S.at(i, i) = 1.0;
		b[i] = 10.0 * (i + 1);
	}
	for (size_t i = 0; i < n - 1; i++) {
		S.at(i, i + 1) = 0.1;
		S.at(i + 1, i) = 0.1;
	}
	matrixdb L = cholesky_(S);
	vectordb x = backward_substitution_(L, b);

	// Test
	EXPECT_NEAR((L.transpose() * x - b).vnorm(), 0.0, EPS);
}
TEST(TestLinalg, SolveCholeskyVector) {
	// Init
	size_t n = 10;
	matrixdb S(n, n, 0.0);
	vectordb b(n);
	for (size_t i = 0; i < n; i++) {
		S.at(i, i) = 1.0;
		b[i] = 10.0 * (i + 1);
	}
	for (size_t i = 0; i < n - 1; i++) {
		S.at(i, i + 1) = 0.1;
		S.at(i + 1, i) = 0.1;
	}
	matrixdb L = cholesky_(S);
	vectordb x = solve_cholesky_(L, b);

	// Test
	EXPECT_NEAR((S * x - b).vnorm(), 0.0, EPS);
}
TEST(TestLinalg, SolveCholeskyMatrix) {
	// Init
	size_t n = 10;
	size_t m = 6;
	matrixdb S(n, n, 0.0);
	matrixdb B(n, m);
	for (size_t i = 0; i < n; i++) {
		S.at(i, i) = 1.0;
		for (size_t j = 0; j < m; j++) {
			B.at(i, j) = 10.0 * (i + 1) * m + 3.0 * m;
		}
	}
	for (size_t i = 0; i < n - 1; i++) {
		S.at(i, i + 1) = 0.1;
		S.at(i + 1, i) = 0.1;
	}
	matrixdb L = cholesky_(S);
	matrixdb X = solve_cholesky_(L, B);

	// Test
	EXPECT_NEAR(frobenius_norm_(S * X - B) / (n * m), 0.0, EPS);
}
TEST(TestLinalg, SolveCholeskySymTridiagMatrix) {
	// Init
	size_t n = 3;
	matrixdb diag(n, n, 0.0);
	vectordb b(2 * n);
	for (size_t i = 0; i < n; i++) {
		diag.at(i, i) = 1.0;
		b[i] = 10.0 * (i + 1);
		b[i + n] = 10.0 * (i + 1 + n);
	}
	for (size_t i = 0; i < n - 1; i++) {
		diag.at(i, i + 1) = 0.1;
		diag.at(i + 1, i) = 0.1;
	}
	vector<matrixdb> list_D{ diag, diag };
	vector<matrixdb> list_E{ 0.001 * matrixdb(n, n, 1.0) };
	sym_tridiag_matrixdb tridiag_S(list_D, list_E);
	sym_tridiag_matrixdb tridiag_L = cholesky_(tridiag_S);
	vectordb x = solve_cholesky_(tridiag_L, b);
	matrixdb S = sym_tridiag_matrixdb_2_matrixdb_(tridiag_S);
	matrixdb L = cholesky_(S);
	vectordb x_expected = solve_cholesky_(cholesky_(S), b);

	// Test
	EXPECT_NEAR((x_expected - x).vnorm(), 0.0, EPS);
}
TEST(TestLinalg, JacobiEigenvalue) {
	// Init DACE
	size_t nb_variables = 3;
	DA::init(2, nb_variables);

	// Init
	size_t n = 10;
	matrixdb S(n, n, 0.0);
	for (size_t i = 0; i < n; i++) {
		S.at(i, i) = 1.0;
	}
	for (size_t i = 0; i < n - 1; i++) {
		S.at(i, i + 1) = 0.1;
		S.at(i + 1, i) = 0.1;
	}
	matrixDA S_DA = S + 0*DA(1);
	pair<vectordb, matrixdb> eig_pair = jacobi_eigenvalue_(S);
	pair<vectorDA, matrixDA> eig_pair_DA = jacobi_eigenvalue_(S_DA);
	matrixdb D = make_diag_matrix_(eig_pair.first);
	matrixDA D_DA = make_diag_matrix_(eig_pair_DA.first);
	matrixdb P = eig_pair.second;
	matrixDA P_DA = eig_pair_DA.second;

	// Test
	EXPECT_NEAR(frobenius_norm_(D * P - P * S), 0.0, EPS);
	EXPECT_NEAR(frobenius_norm_((D_DA * P_DA - P_DA * S_DA).cons()), 0.0, EPS);
}
