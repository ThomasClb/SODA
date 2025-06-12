/**
	transcription.cpp

	Purpose: Implementation of the transcription methods.

	@author Thomas Caleb

	@version 2.0 26/11/2024
*/

#include "transcription.h"

using namespace DACE;
using namespace std;


// Adaptive first order transcription for Newton method.
// DOI: 10.48550/arXiv.2502.15949
vectorDA first_order_transcription(
	vectorDA const& constraints_eval,
	vector<matrixdb> const& list_Sigma, vector<matrixdb> const& list_feedback_gain,
	double const& eta,
	SpacecraftParameters const& spacecraft_parameters, Constants const& constants,
	SolverParameters const& solver_parameters) {
	// Unpack
	size_t N = solver_parameters.N();
	size_t Nx = solver_parameters.Nx();
	size_t Nu = solver_parameters.Nu();
	size_t Nineq = solver_parameters.Nineq();
	size_t Ntineq = solver_parameters.Ntineq();
	double tol = solver_parameters.PN_tol();
	double eps(EPS);
	size_t d = N*Nineq + Ntineq;
	double quantile = solver_parameters.path_quantile();
	double beta_star = 1 - chi_2_cdf(d, sqr(quantile));

	// Init
	vectorDA output(constraints_eval);
	vectorDA norm_vector(d);
	matrixDA D_f, prod_k;
	vector<matrixDA> list_der;

	// Loop on all steps

	if (Nineq != 0) {
		for (size_t k=0; k<N; k++) {

			// Evaluate constraints derivatives
			list_der = deriv_xu_DA(
				constraints_eval.extract(
					k*(Nineq + 0), (k + 1)*Nineq - 1), Nx, Nu, false);

			// Compute the contraints covariance
			D_f = (list_der[0] + list_der[1] * list_feedback_gain[k]);
			prod_k = D_f*list_Sigma[k];

			// Assign
			for (size_t i=0; i<Nineq; i++) {
				norm_vector[k*Nineq + i] = vectorDA(prod_k.getrow(i)).dot(vectorDA(D_f.getrow(i)));
			}
		}
	}

	// Terminal constraints

	// Evaluate constraints derivatives
	if (Ntineq != 0) {
		list_der = deriv_x_DA(
					constraints_eval.extract(
						N*Nineq, N*Nineq + Ntineq - 1), Nx, false);
				
		// Compute the contraints covariance
		D_f = list_der[0];
		prod_k = D_f*list_Sigma[N];
	}

	// Assign
	for (size_t i=0; i<Ntineq; i++) {
		norm_vector[N*Nineq + i] = vectorDA(prod_k.getrow(i)).dot(vectorDA(D_f.getrow(i)));
	}

	// Build output vector
	quantile = sqrt(inv_chi_2_cdf((int)(d*eta), 1 - beta_star));
	for (size_t k=0; k<N; k++) {
		// Assign
		for (size_t i=0; i<Nineq; i++) {
			size_t index(k*Nineq + i);
			DA norm_vector_i(norm_vector[index]);
			double norm_vector_cons_i(norm_vector_i.cons());

			if (norm_vector_cons_i > 0) {
				output[index] += quantile*sqrt(norm_vector_i);
			}
			else if (norm_vector_cons_i == 0 && output[index].cons() <= 0)
				output[index] = -1e15;
		}
	}

	// Terminal constraints

	// Assign
	for (size_t i=0; i<Ntineq; i++) {	
		size_t index(N*Nineq + i);
		DA norm_vector_i(norm_vector[index]);
		double norm_vector_cons_i(norm_vector_i.cons());

		if (norm_vector_cons_i > 0) {
			output[index] += quantile*sqrt(norm_vector_i);
		}
		else if (norm_vector_cons_i == 0 && output[index].cons() <= 0)
			output[index] = -1e15;
	}

	return output;
}

// First order method for path constraints.
// DOI: 10.48550/arXiv.2502.15949
vectorDA first_order_path_transcription(
	vectorDA const& constraints_eval, stateDA const& x_DA, controlDA const& u_DA,
	SpacecraftParameters const& spacecraft_parameters, Constants const& constants,
	SolverParameters const& solver_parameters) {
	// Unpack
	size_t Nx = solver_parameters.Nx();
	size_t Nu = solver_parameters.Nu();
	vectordb rho_parameters = solver_parameters.backward_sweep_regulation_parameters();
	double rho_min = rho_parameters[1];
	size_t d = constraints_eval.size();

	// Evaluate constraints derivatives
	vector<matrixDA> list_der = deriv_xu_DA(
		constraints_eval, Nx, Nu, false);

	// Compute the contraints covariance
	matrixDA D_f = (list_der[0] + list_der[1] * u_DA.feedback_gain());
	matrixDA Sigma = D_f*x_DA.Sigma()*D_f.transpose();
	vectorDA norm_vector(d);
	for (size_t i=0; i<d; i++) {
		if (Sigma.at(i,i).cons() != 0.0)
			norm_vector[i] = sqrt(Sigma.at(i,i));
	}

	// Make final constraints 
	vectorDA transcribed_constraints = constraints_eval + solver_parameters.path_quantile()*norm_vector;
	return transcribed_constraints;
}

// First order method for terminal constraints.
// DOI: 10.48550/arXiv.2502.15949
vectorDA first_order_terminal_transcription(
	vectorDA const& constraints_eval, stateDA const& x_DA, statedb const& x_goal,
	SpacecraftParameters const& spacecraft_parameters, Constants const& constants,
	SolverParameters const& solver_parameters) {
	// Unpack
	size_t Nx = solver_parameters.Nx();
	vectordb rho_parameters = solver_parameters.backward_sweep_regulation_parameters();
	double rho_min = rho_parameters[1];
	size_t d = constraints_eval.size();

	// Evaluate constraints derivatives
	vector<matrixDA> list_der = deriv_x_DA(
		constraints_eval, Nx, false);

	// Compute the contraints covariance
	matrixDA D_f = list_der[0];
	matrixDA Sigma = D_f*x_DA.Sigma()*D_f.transpose();
	vectorDA norm_vector(d);
	for (size_t i=0; i<d; i++) {
		if (Sigma.at(i,i).cons() != 0.0) {
			norm_vector[i] = sqrt(Sigma.at(i,i));
		}
	}

	// Make final constraints 
	vectorDA transcribed_constraints = constraints_eval + solver_parameters.terminal_quantile()*norm_vector;
	return transcribed_constraints;
}

// Spectral radius method for path constraints.
// Generalisation of [Ridderhof et al. 2021]
// DOI: 10.48550/arXiv.2502.15949
vectorDA spectral_radius_path_transcription(
	vectorDA const& constraints_eval, stateDA const& x_DA, controlDA const& u_DA,
	SpacecraftParameters const& spacecraft_parameters, Constants const& constants,
	SolverParameters const& solver_parameters) {
	// Unpack
	size_t Nx = solver_parameters.Nx();
	size_t Nu = solver_parameters.Nu();
	vectordb rho_parameters = solver_parameters.backward_sweep_regulation_parameters();
	double rho_min = rho_parameters[1];
	size_t d = constraints_eval.size();
	
	// Evaluate constraints derivatives
	vector<matrixDA> list_der = deriv_xu_DA(
		constraints_eval, Nx, Nu, false);

	// Compute the contraints covariance
	matrixDA D_f = (list_der[0] + list_der[1] * u_DA.feedback_gain());
	matrixDA Sigma = D_f*x_DA.Sigma()*D_f.transpose();
	double norm_Sigma = frobenius_norm_(Sigma.cons());
	Sigma = (1.0/norm_Sigma)*Sigma;

	// Get spectral radius
	vectorDA eigenvalues = norm_Sigma*(jacobi_eigenvalue_(Sigma).first);
	DA spectral_radius = 0; double max = 0;
	for (size_t i=0; i<d; i++) {
		if (eigenvalues[i].cons() > max)
			spectral_radius = eigenvalues[i];
	}
	DA sqrt_spectral_radius(0.0);
	if (abs_cons(spectral_radius) != 0) {
		sqrt_spectral_radius = sqrt(spectral_radius);
	}

	// Make boundary
	vectorDA norm_vector(d, sqrt_spectral_radius);

	// Make final constraints 
	vectorDA transcribed_constraints = constraints_eval + solver_parameters.path_quantile()*norm_vector;
	return transcribed_constraints;
}

// Spectral radius method for terminal constraints.
// Generalisation of [Ridderhof et al. 2021]
// DOI: 10.48550/arXiv.2502.15949
vectorDA spectral_radius_terminal_transcription(
	vectorDA const& constraints_eval, stateDA const& x_DA, statedb const& x_goal,
	SpacecraftParameters const& spacecraft_parameters, Constants const& constants,
	SolverParameters const& solver_parameters) {
	// Unpack
	size_t Nx = solver_parameters.Nx();
	size_t d = constraints_eval.size();

	// Evaluate constraints derivatives
	vector<matrixDA> list_der = deriv_x_DA(
		constraints_eval, Nx, false);

	// Compute the contraints covariance
	matrixDA D_f = list_der[0];
	matrixDA Sigma = D_f*x_DA.Sigma()*D_f.transpose();
	double norm_Sigma = frobenius_norm_(Sigma.cons());
	Sigma = (1.0/norm_Sigma)*Sigma;

	// Get spectral radius
	vectorDA eigenvalues = norm_Sigma*(jacobi_eigenvalue_(Sigma).first);
	DA spectral_radius = 0; double max = 0;
	for (size_t i=0; i<d; i++) {
		if (eigenvalues[i].cons() > max)
			spectral_radius = eigenvalues[i];
	}
	DA sqrt_spectral_radius(0.0);
	if (abs_cons(spectral_radius) != 0) {
		sqrt_spectral_radius = sqrt(spectral_radius);
	}

	// Make boundary
	vectorDA norm_vector(d, sqrt_spectral_radius);

	// Make final constraints 
	vectorDA transcribed_constraints = constraints_eval + solver_parameters.terminal_quantile()*norm_vector;
	return transcribed_constraints;
}
