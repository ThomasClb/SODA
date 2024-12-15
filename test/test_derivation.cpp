/**
	test_derivation.cpp

	Purpose: Test of the implementation of the DA-based
	automatic differentiation methods.

	@author Thomas Caleb

	@version 1.0 05/09/2024
*/

#include "gtest/gtest.h"
#include "../include/soda.h"

using namespace DACE;
using namespace std;

// Test of db_2_DA
TEST(TestDerivation, db2DA) {
	// Init DACE
	size_t nb_variables = 3;
	DA::init(2, nb_variables);

	// Init
	size_t n = 10;
	vectordb vector_db(n, 0.0);
	for (size_t i = 0; i < n; i++) {
		vector_db[i] = 10.0 * i;
	}
	vectorDA vector_DA = db_2_DA(vector_db);
	
	// Tests
	for (size_t i = 0; i < n; i++) {
		double norm_i = (vector_db[i] - vector_DA[i]).norm();
		EXPECT_EQ(norm_i, 0.0);
	}
}

// Test of id_vector
TEST(TestDerivation, idvector) {
	// Init DACE
	size_t nb_variables = 6;
	DA::init(2, nb_variables);

	// Init
	size_t n = 10;
	vectordb vector_db(n, 0.0);
	size_t begin_index = 3;
	size_t begin_da_index = 2;
	size_t end_da_index = 5;
	for (size_t i = 0; i < n; i++) {
		vector_db[i] = 10.0 * i;
	}
	vectorDA vector_DA = id_vector(
		vector_db, begin_index,
		begin_da_index, end_da_index);

	// Tests
	vectorDA diff = vector_DA - vector_db;
	vectordb dx_eval(nb_variables, 1.0);
	for (size_t i = 0; i < nb_variables; i++) {
		dx_eval[i] = 1.0 * i;
	}
	vectordb diff_eval = diff.eval(dx_eval);
	for (size_t i = 0; i < n; i++) {
		if (i < begin_index || i >= begin_index + end_da_index - begin_da_index)
			EXPECT_EQ(diff_eval[i], 0.0);
		else 
			EXPECT_EQ(diff_eval[i], i - begin_index + begin_da_index);
	}
}
TEST(TestDerivation, idvectorState) {
	// TO DO
}
TEST(TestDerivation, idvectorControl) {
	// TO DO
}

// Differentiation functions

// Test of deriv_xu
TEST(TestDerivation, derivxu) {
	// Init DACE
	size_t nb_variables = 5;
	DA::init(2, nb_variables);

	// Init
	size_t Nx = 3;
	size_t Nu = 2;
	size_t n = nb_variables;
	vectordb f(n);
	for (size_t i = 0; i < n; i++) {
		f[i] = 1.0 * (i + 1.0);
	}
	vectorDA f_eval = id_vector(f, 0, 0, nb_variables);
	f_eval = sqr(f_eval);
	bool hessian_computation = false;
	vector<matrixdb> list_der = deriv_xu(f_eval, Nx, Nu, hessian_computation);	

	// Tests

	// Check dimensions
	matrixdb der_x = list_der[0];
	matrixdb der_u = list_der[1];
	EXPECT_EQ(der_x.nrows(), n);
	EXPECT_EQ(der_x.ncols(), Nx);
	EXPECT_EQ(der_u.nrows(), n);
	EXPECT_EQ(der_u.ncols(), Nu);

	// Check values
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < Nx; j++) {
			if (i == j) {
				EXPECT_EQ(der_x.at(i, j), 2*(i + 1.0));
			}
			else {
				EXPECT_EQ(der_x.at(i, j), 0.0);
			}
		}
	}
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < Nu; j++) {
			if (i - Nx == j) {
				EXPECT_EQ(der_u.at(i, j), 2 * (i + 1.0));
			}
			else {
				EXPECT_EQ(der_u.at(i, j), 0.0);
			}
		}
	}

	// With hessian
	hessian_computation = true;
	list_der = deriv_xu(f_eval, Nx, Nu, hessian_computation);

	// Tests

	// Check dimensions
	der_x = list_der[0];
	der_u = list_der[1];
	EXPECT_EQ(der_x.nrows(), n);
	EXPECT_EQ(der_x.ncols(), Nx);
	EXPECT_EQ(der_u.nrows(), n);
	EXPECT_EQ(der_u.ncols(), Nu);
	for (size_t i = 0; i < n; i++) {
		matrixdb hess_xx_i = list_der[2 + i];
		matrixdb hess_uu_i = list_der[2 + n + i];
		matrixdb hess_ux_i = list_der[2 + 2 * n + i];
		EXPECT_EQ(hess_xx_i.nrows(), Nx);
		EXPECT_EQ(hess_xx_i.ncols(), Nx);
		EXPECT_EQ(hess_uu_i.nrows(), Nu);
		EXPECT_EQ(hess_uu_i.ncols(), Nu);
		EXPECT_EQ(hess_ux_i.nrows(), Nu);
		EXPECT_EQ(hess_ux_i.ncols(), Nx);
	}

	// Check values

	// Gradients
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < Nx; j++) {
			if (i == j) {
				EXPECT_EQ(der_x.at(i, j), 2 * (i + 1.0));
			}
			else {
				EXPECT_EQ(der_x.at(i, j), 0.0);
			}
		}
	}
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < Nu; j++) {
			if (i - Nx == j) {
				EXPECT_EQ(der_u.at(i, j), 2 * (i + 1.0));
			}
			else {
				EXPECT_EQ(der_u.at(i, j), 0.0);
			}
		}
	}

	// Hessians
	for (size_t k = 0; k < n; k++) {
		matrixdb hess_xx_k = list_der[2 + k];
		matrixdb hess_uu_k = list_der[2 + n + k];
		matrixdb hess_ux_k = list_der[2 + 2 * n + k];

		// fxx
		for (size_t i = 0; i < Nx; i++) {
			for (size_t j = 0; j < Nx; j++) {
				if (i == j && i == k) {
					EXPECT_EQ(hess_xx_k.at(i, j), 2.0);
				}
				else {
					EXPECT_EQ(hess_xx_k.at(i, j), 0.0);
				}
			}
		}

		// fuu
		for (size_t i = 0; i < Nu; i++) {
			for (size_t j = 0; j < Nu; j++) {
				if (i == j && (i + Nx) == k) {
					EXPECT_EQ(hess_uu_k.at(i, j), 2.0);
				}
				else {
					EXPECT_EQ(hess_uu_k.at(i, j), 0.0);
				}
			}
		}

		// fux
		for (size_t i = 0; i < Nu; i++) {
			for (size_t j = 0; j < Nx; j++) {
				EXPECT_EQ(hess_ux_k.at(i, j), 0.0);
			}
		}
	}
}

// Test of deriv_xu_DA
TEST(TestDerivation, derivxuDA) {
	// Init DACE
	size_t nb_variables = 5;
	DA::init(2, nb_variables);

	// Init
	size_t Nx = 3;
	size_t Nu = 2;
	size_t n = nb_variables;
	vectordb f(n);
	for (size_t i = 0; i < n; i++) {
		f[i] = 1.0 * (i + 1.0);
	}
	vectorDA f_eval = id_vector(f, 0, 0, nb_variables);
	f_eval = sqr(f_eval);
	bool hessian_computation = false;
	vector<matrixDA> list_der = deriv_xu_DA(f_eval, Nx, Nu, hessian_computation);	

	// Tests

	// Check dimensions
	matrixDA der_x = list_der[0];
	matrixDA der_u = list_der[1];
	EXPECT_EQ(der_x.nrows(), n);
	EXPECT_EQ(der_x.ncols(), Nx);
	EXPECT_EQ(der_u.nrows(), n);
	EXPECT_EQ(der_u.ncols(), Nu);

	// Check values
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < Nx; j++) {
			if (i == j) {
				EXPECT_EQ((der_x.at(i, j) - 2*(i + 1.0)).cons(), 0);
			}
			else {
				EXPECT_EQ(der_x.at(i, j).cons(), 0.0);
			}
		}
	}
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < Nu; j++) {
			if (i - Nx == j) {
				EXPECT_EQ((der_u.at(i, j) - 2 * (i + 1.0)).cons(), 0);
			}
			else {
				EXPECT_EQ(der_u.at(i, j).cons(), 0.0);
			}
		}
	}

	// With hessian
	hessian_computation = true;
	list_der = deriv_xu_DA(f_eval, Nx, Nu, hessian_computation);

	// Tests

	// Check dimensions
	der_x = list_der[0];
	der_u = list_der[1];
	EXPECT_EQ(der_x.nrows(), n);
	EXPECT_EQ(der_x.ncols(), Nx);
	EXPECT_EQ(der_u.nrows(), n);
	EXPECT_EQ(der_u.ncols(), Nu);
	for (size_t i = 0; i < n; i++) {
		matrixDA hess_xx_i = list_der[2 + i];
		matrixDA hess_uu_i = list_der[2 + n + i];
		matrixDA hess_ux_i = list_der[2 + 2 * n + i];
		EXPECT_EQ(hess_xx_i.nrows(), Nx);
		EXPECT_EQ(hess_xx_i.ncols(), Nx);
		EXPECT_EQ(hess_uu_i.nrows(), Nu);
		EXPECT_EQ(hess_uu_i.ncols(), Nu);
		EXPECT_EQ(hess_ux_i.nrows(), Nu);
		EXPECT_EQ(hess_ux_i.ncols(), Nx);
	}

	// Check values

	// Gradients
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < Nx; j++) {
			if (i == j) {
				EXPECT_EQ((der_x.at(i, j) - 2*(i + 1.0)).cons(), 0);
			}
			else {
				EXPECT_EQ(der_x.at(i, j).cons(), 0.0);
			}
		}
	}
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < Nu; j++) {
			if (i - Nx == j) {
				EXPECT_EQ((der_u.at(i, j) - 2 * (i + 1.0)).cons(), 0);
			}
			else {
				EXPECT_EQ(der_u.at(i, j).cons(), 0.0);
			}
		}
	}

	// Hessians
	for (size_t k = 0; k < n; k++) {
		matrixDA hess_xx_k = list_der[2 + k];
		matrixDA hess_uu_k = list_der[2 + n + k];
		matrixDA hess_ux_k = list_der[2 + 2 * n + k];

		// fxx
		for (size_t i = 0; i < Nx; i++) {
			for (size_t j = 0; j < Nx; j++) {
				if (i == j && i == k) {
					EXPECT_EQ((hess_xx_k.at(i, j) - 2.0).cons(), 0);
				}
				else {
					EXPECT_EQ(hess_xx_k.at(i, j).cons(), 0.0);
				}
			}
		}

		// fuu
		for (size_t i = 0; i < Nu; i++) {
			for (size_t j = 0; j < Nu; j++) {
				if (i == j && (i + Nx) == k) {
					EXPECT_EQ((hess_uu_k.at(i, j)- 2.0).cons(), 0);
				}
				else {
					EXPECT_EQ(hess_uu_k.at(i, j).cons(), 0.0);
				}
			}
		}

		// fux
		for (size_t i = 0; i < Nu; i++) {
			for (size_t j = 0; j < Nx; j++) {
				EXPECT_EQ(hess_ux_k.at(i, j).cons(), 0.0);
			}
		}
	}
}

// Test of deriv_x
TEST(TestDerivation, derivx) {
	// Init DACE
	size_t nb_variables = 5;
	DA::init(2, nb_variables);

	// Init
	size_t Nx = 3;
	size_t n = nb_variables;
	vectordb f(n);
	for (size_t i = 0; i < n; i++) {
		f[i] = 1.0 * (i + 1.0);
	}
	vectorDA f_eval = id_vector(f, 0, 0, nb_variables);
	f_eval = sqr(f_eval);
	bool hessian_computation = false;
	vector<matrixdb> list_der = deriv_x(f_eval, Nx, hessian_computation);

	// Tests

	// Check dimensions
	matrixdb der_x = list_der[0];
	EXPECT_EQ(der_x.nrows(), n);
	EXPECT_EQ(der_x.ncols(), Nx);

	// Check values
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < Nx; j++) {
			if (i == j) {
				EXPECT_EQ(der_x.at(i, j), 2 * (i + 1.0));
			}
			else {
				EXPECT_EQ(der_x.at(i, j), 0.0);
			}
		}
	}

	// With hessian
	hessian_computation = true;
	list_der = deriv_x(f_eval, Nx, hessian_computation);

	// Tests

	// Check dimensions
	der_x = list_der[0];
	EXPECT_EQ(der_x.nrows(), n);
	EXPECT_EQ(der_x.ncols(), Nx);
	for (size_t i = 0; i < n; i++) {
		matrixdb hess_xx_i = list_der[1 + i];
		EXPECT_EQ(hess_xx_i.nrows(), Nx);
		EXPECT_EQ(hess_xx_i.ncols(), Nx);
	}

	// Check values

	// Gradients
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < Nx; j++) {
			if (i == j) {
				EXPECT_EQ(der_x.at(i, j), 2 * (i + 1.0));
			}
			else {
				EXPECT_EQ(der_x.at(i, j), 0.0);
			}
		}
	}

	// Hessians
	for (size_t k = 0; k < n; k++) {
		matrixdb hess_xx_k = list_der[1 + k];

		// fxx
		for (size_t i = 0; i < Nx; i++) {
			for (size_t j = 0; j < Nx; j++) {
				if (i == j && i == k) {
					EXPECT_EQ(hess_xx_k.at(i, j), 2.0);
				}
				else {
					EXPECT_EQ(hess_xx_k.at(i, j), 0.0);
				}
			}
		}
	}
}

// Test of deriv_x_DA
TEST(TestDerivation, derivxDA) {
	// Init DACE
	size_t nb_variables = 5;
	DA::init(2, nb_variables);

	// Init
	size_t Nx = 3;
	size_t n = nb_variables;
	vectordb f(n);
	for (size_t i = 0; i < n; i++) {
		f[i] = 1.0 * (i + 1.0);
	}
	vectorDA f_eval = id_vector(f, 0, 0, nb_variables);
	f_eval = sqr(f_eval);
	bool hessian_computation = false;
	vector<matrixDA> list_der = deriv_x_DA(f_eval, Nx, hessian_computation);

	// Tests

	// Check dimensions
	matrixDA der_x = list_der[0];
	EXPECT_EQ(der_x.nrows(), n);
	EXPECT_EQ(der_x.ncols(), Nx);

	// Check values
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < Nx; j++) {
			if (i == j) {
				EXPECT_EQ((der_x.at(i, j) - 2 * (i + 1.0)).cons(), 0);
			}
			else {
				EXPECT_EQ(der_x.at(i, j).cons(), 0.0);
			}
		}
	}

	// With hessian
	hessian_computation = true;
	list_der = deriv_x_DA(f_eval, Nx, hessian_computation);

	// Tests

	// Check dimensions
	der_x = list_der[0];
	EXPECT_EQ(der_x.nrows(), n);
	EXPECT_EQ(der_x.ncols(), Nx);
	for (size_t i = 0; i < n; i++) {
		matrixDA hess_xx_i = list_der[1 + i];
		EXPECT_EQ(hess_xx_i.nrows(), Nx);
		EXPECT_EQ(hess_xx_i.ncols(), Nx);
	}

	// Check values

	// Gradients
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < Nx; j++) {
			if (i == j) {
				EXPECT_EQ((der_x.at(i, j)- 2 * (i + 1.0)).cons(), 0);
			}
			else {
				EXPECT_EQ(der_x.at(i, j).cons(), 0.0);
			}
		}
	}

	// Hessians
	for (size_t k = 0; k < n; k++) {
		matrixDA hess_xx_k = list_der[1 + k];

		// fxx
		for (size_t i = 0; i < Nx; i++) {
			for (size_t j = 0; j < Nx; j++) {
				if (i == j && i == k) {
					EXPECT_EQ((hess_xx_k.at(i, j) - 2.0).cons(), 0);
				}
				else {
					EXPECT_EQ(hess_xx_k.at(i, j).cons(), 0.0);
				}
			}
		}
	}
}