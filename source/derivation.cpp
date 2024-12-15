/**
	derivation.cpp

	Purpose: Implementation of the DA-based
	automatic differentiation methods.

	@author Thomas Caleb

	@version 2.0 02/09/2024
*/

#include "derivation.h"

using namespace DACE;
using namespace std;

double abs_cons(double const& value) {return abs(value);}
double abs_cons(DA const& value) {return abs(value.cons());}

// DACE vector conversion

// Turns a vectordb into a vectorDA without perturbation.
vectorDA db_2_DA(
	vectordb const& vector_db) {
	size_t n = vector_db.size();
	vectorDA vector_DA; vector_DA.reserve(n);
	for (size_t i = 0; i <n; i++) {
		vector_DA.push_back(vector_db[i]);
	}
	return vector_DA;
}

// Turns a vectordb into a vectorDA with identity perturbation.
vectorDA id_vector(
	vectordb const& vector_db,
	unsigned int const& begin_index,
	unsigned int const& begin_da_index,
	unsigned int const& end_da_index) {
	
	
	size_t n = vector_db.size();
	vectorDA vector_DA; vector_DA.reserve(vector_db.size());
	for (unsigned int i = 0; i < n; i++) {

		// Copy constant part
		DA coord_i = vector_db[i];

		// Add linear terms starting from begin_index
		// From variable begin_da_index to end_da_index
		if (i >= begin_index && i < begin_index + end_da_index - begin_da_index)
			coord_i += DA(begin_da_index + i - begin_index + 1);

		// Assign
		vector_DA.push_back(coord_i);
	}

	return vector_DA;
}

// Turns a statedb into a stateDA with identity perturbation.
stateDA id_vector(
	statedb const& state_db,
	unsigned int const& begin_index,
	unsigned int const& begin_da_index,
	unsigned int const& end_da_index) {
	stateDA state_DA(id_vector(state_db.nominal_state(), begin_index, begin_da_index, end_da_index));
	state_DA.set_Sigma(state_db.Sigma());
	state_DA.set_der_dynamics(state_db.der_dynamics());
	return state_DA;
}

// Turns a controldb into a controlDA with identity perturbation.
controlDA id_vector(
	controldb const& control_db,
	unsigned int const& begin_index,
	unsigned int const& begin_da_index,
	unsigned int const& end_da_index) {
	controlDA control_DA(id_vector(control_db.nominal_control(), begin_index, begin_da_index, end_da_index));
	control_DA.set_feedback_gain(control_db.feedback_gain());
	return control_DA;
}


// Differentiation functions

// Differentiates with respect to x, and u.
vector<matrixdb> deriv_xu(
	vectorDA const& f_eval,
	unsigned int const& Nx, unsigned int const& Nu,
	bool const& hessian_computation) {
	// Unpack
	size_t n = f_eval.size();

	if (!hessian_computation) {
		matrixdb fx, fu;
		if (n == 0) {
			fx = matrixdb(n,Nx);
			fu = matrixdb(n,Nu);
		} else {

			// Get global jacobian
			matrixdb jac_f = f_eval.linear();

			// First order derivatives
			fx = jac_f.submat(0, 0, n - 1, Nx - 1);
			fu = jac_f.submat(0, Nx, n - 1, Nx + Nu - 1);
		}
		return vector<matrixdb>{ fx, fu };
	}
	else {
		matrixdb fx, fu;
		vector<matrixdb> output(2 + 3 * n);
		if (n == 0) {
			fx = matrixdb(n,Nx);
			fu = matrixdb(n,Nu);
			output[0] = fx;
			output[1] = fu;
		} else {

			// Build DA jacobian
			matrixDA df_eval(n, Nx + Nu);
			for (size_t i = 0; i < Nx + Nu; i++) {
				df_eval.setcol(i, f_eval.deriv(i + 1));
			}
			matrixdb jac_f = df_eval.cons();

			// First order derivatives
			fx = jac_f.submat(0, 0, n - 1, Nx - 1);
			fu = jac_f.submat(0, Nx, n - 1, Nx + Nu - 1);

			// Assign
			output[0] = fx; output[1] = fu;

			// Second order derivatives
			for (size_t i = 0; i < n; i++) {
				matrixdb hessian = vectorDA(df_eval.getrow(i)).linear();

				// Slice hesssians
				output[2 + i] = hessian.submat(0, 0, Nx - 1, Nx - 1); // fxx
				output[2 + i + n] = hessian.submat(Nx, Nx, Nx + Nu - 1, Nx + Nu - 1); // fuu
				output[2 + i + 2*n] = hessian.submat(Nx, 0, Nx + Nu - 1, Nx - 1); // fux
			}
		}
		return output;
	}
}

// Differentiates with respect to x, and u. Returns DA matrices.
vector<matrixDA> deriv_xu_DA(
	vectorDA const& f_eval,
	unsigned int const& Nx, unsigned int const& Nu,
	bool const& hessian_computation) {
	// Unpack
	size_t n = f_eval.size();

	if (!hessian_computation) {
		if (n == 0) {
			matrixDA fx(n,Nx);
			matrixDA fu(n,Nu);
			return vector<matrixDA>{fx, fu};
		}

		// Get global jacobian
		matrixDA jac_f(n, Nx+Nu);
		for (size_t i = 0; i < Nx + Nu; i++) {
			jac_f.setcol(i, f_eval.deriv(i + 1));
		}

		// First order derivatives
		matrixDA fx = jac_f.submat(0, 0, n - 1, Nx - 1);
		matrixDA fu = jac_f.submat(0, Nx, n - 1, Nx + Nu - 1);

		return vector<matrixDA>{ fx, fu };
	}
	else {
		matrixDA fx, fu;
		vector<matrixDA> output(2 + 3 * n);
		if (n == 0) {
			fx = matrixDA(n,Nx);
			fu = matrixDA(n,Nu);
			output[0] = fx;
			output[1] = fu;
		} else {

			// Build DA jacobian
			matrixDA jac_f(n, Nx + Nu);
			for (size_t i = 0; i < Nx + Nu; i++) {
				jac_f.setcol(i, f_eval.deriv(i + 1));
			}

			// First order derivatives
			fx = jac_f.submat(0, 0, n - 1, Nx - 1);
			fu = jac_f.submat(0, Nx, n - 1, Nx + Nu - 1);

			// Assign
			output[0] = fx; output[1] = fu;

			// Second order derivatives
			for (size_t i = 0; i < n; i++) {
				// Build hessian
				matrixDA hessian(Nx + Nu, Nx + Nu);
				vectorDA jac_f_i(jac_f.getrow(i));
				for (size_t j = 0; j < Nx + Nu; j++) {
					hessian.setrow(j, jac_f_i.deriv(j + 1));
				}

				// Slice hesssians
				output[2 + i] = hessian.submat(0, 0, Nx - 1, Nx - 1); // fxx
				output[2 + i + n] = hessian.submat(Nx, Nx, Nx + Nu - 1, Nx + Nu - 1); // fuu
				output[2 + i + 2*n] = hessian.submat(Nx, 0, Nx + Nu - 1, Nx - 1); // fux
			}
		}
		return output;
	}
}

// Differentiates with respect to x.
vector<matrixdb> deriv_x(
	vectorDA const& f_eval,
	unsigned int const& Nx,
	bool const& hessian_computation) {
	// Unpack
	size_t n = f_eval.size();

	if (!hessian_computation) {

		if (n == 0) {
			matrixdb fx(n,Nx);
			return vector<matrixdb>{fx};
		}

		// Get global jacobian
		matrixdb jac_f = f_eval.linear();

		// First order derivatives
		matrixdb fx = jac_f.submat(0, 0, n - 1, Nx - 1);
		return vector<matrixdb>{fx};
	}
	else {
		// Init
		vector<matrixdb> output; 
		output.reserve(1 + n);

		// Build DA jacobian
		matrixDA df_eval(n, Nx);
		for (size_t i = 0; i < Nx; i++) {
			df_eval.setcol(i, f_eval.deriv(i + 1));
		}
		matrixdb jac_f = df_eval.cons();

		// First order derivatives
		output.emplace_back(jac_f);

		// Add hessians
		for (size_t i = 0; i < n; i++) {
			output.emplace_back(vectorDA(
				df_eval.getrow(i)
			).linear().submat(Nx - 1, Nx - 1));
		}
		return output;
	}
}

// Differentiates with respect to x. Returns DA matrices.
vector<matrixDA> deriv_x_DA(
	vectorDA const& f_eval,
	unsigned int const& Nx,
	bool const& hessian_computation) {
	// Unpack
	size_t n = f_eval.size();

	if (!hessian_computation) {
		if (n == 0) {
			return vector<matrixDA>(1, matrixDA(n,Nx));
		}

		// Get global jacobian
		matrixDA jac_f(n, Nx);
		for (size_t i = 0; i < Nx; i++) {
			jac_f.setcol(i, f_eval.deriv(i + 1));
		}
		return vector<matrixDA>{jac_f};
	}
	else {
		if (f_eval.size() == 0) {
			return vector<matrixDA>(1, matrixDA(n,Nx));
		}

		// Make output
		vector<matrixDA> output(1 + n);

		// Build DA jacobian
		matrixDA jac_f(n, Nx);
		for (size_t i = 0; i < Nx; i++) {
			jac_f.setcol(i, f_eval.deriv(i + 1));
		}

		// First order derivatives
		output[0] = jac_f;

		// Second order derivatives
		for (size_t i = 0; i < n; i++) {
			// Build and assign hessian
			matrixDA hessian(Nx, Nx);
			vectorDA jac_f_i(jac_f.getrow(i));
			for (size_t j = 0; j < Nx; j++) {
				hessian.setrow(j, jac_f_i.deriv(j + 1));
			}
			output[1 + i] = hessian;
		}
		return output;
	}
}