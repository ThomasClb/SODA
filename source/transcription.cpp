/**
	transcription.cpp

	Purpose: Implementation of the transcription methods.

	@author Thomas Caleb

	@version 2.0 26/11/2024
*/

#include "transcription.h"

using namespace DACE;
using namespace std;

// First order method.
// DOI: WIP
// Path constraints.
// DDP version.
vectorDA first_order_path_inequality_transcription(
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
		if (Sigma.at(i,i).cons() != 0.0) {
			norm_vector[i] = sqrt(Sigma.at(i,i));
		}
	}

	// Make final constraints 
	vectorDA transcribed_constraints = constraints_eval + solver_parameters.path_quantile()*norm_vector;
	return transcribed_constraints;
}

// First order method.
// DOI: WIP
// Path constraints.
// PN version.
vectorDA first_order_path_inequality_transcription(
	vectorDA const& constraints_eval,
	matrixdb const& Sigma_k, matrixdb const& feedback_gain,
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
	matrixDA D_f = (list_der[0] + list_der[1] * feedback_gain);
	matrixDA Sigma_kp1 = D_f*Sigma_k*D_f.transpose();
	vectorDA norm_vector(d);
	for (size_t i=0; i<d; i++) {
		if (Sigma_kp1.at(i,i).cons() != 0.0) {
			norm_vector[i] = sqrt(Sigma_kp1.at(i,i));
		}
	}

	// Make final constraints 
	vectorDA transcribed_constraints = constraints_eval + solver_parameters.path_quantile()*norm_vector;
	return transcribed_constraints;
}

// First order method.
// DOI: WIP
// Terminal constraints.
// DDP version.
vectorDA first_order_terminal_inequality_transcription(
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

// First order method.
// DOI: WIP
// Terminal constraints.
// PN version.
vectorDA first_order_terminal_inequality_transcription(
	vectorDA const& constraints_eval, matrixdb const& Sigma_k,
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
	matrixDA Sigma = D_f*Sigma_k*D_f.transpose();
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


// Spectral radius method.
// Generalisation of [Ozaki et al. 2020]
// DOI: https://doi.org/10.2514/1.G004363
// Path constraints.
// DDP version.
vectorDA spectral_radius_path_inequality_transcription(
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

// Spectral radius method.
// Generalisation of [Ozaki et al. 2020]
// DOI: https://doi.org/10.2514/1.G004363
// Path constraints.
// PN version.
vectorDA spectral_radius_path_inequality_transcription(
	vectorDA const& constraints_eval,
	matrixdb const& Sigma_k, matrixdb const& feedback_gain,
	SpacecraftParameters const& spacecraft_parameters, Constants const& constants,
	SolverParameters const& solver_parameters) {
	// Unpack
	size_t Nx = solver_parameters.Nx();
	size_t Nu = solver_parameters.Nu();
	size_t d = constraints_eval.size();

	// Evaluate constraints derivatives
	vector<matrixDA> list_der = deriv_xu_DA(
		constraints_eval, Nx, Nu, false);

	// Compute the contraints covariance
	matrixDA D_f = (list_der[0] + list_der[1] * feedback_gain);
	matrixDA Sigma_kp1 = D_f*Sigma_k*D_f.transpose();
	double norm_Sigma = frobenius_norm_(Sigma_kp1.cons());
	Sigma_kp1 = (1.0/norm_Sigma)*Sigma_kp1;

	// Get spectral radius
	vectorDA eigenvalues = norm_Sigma*(jacobi_eigenvalue_(Sigma_kp1).first);
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

// Spectral radius method.
// Generalisation of [Ozaki et al. 2020]
// DOI: https://doi.org/10.2514/1.G004363
// Terminal constraints.
// DDP version.
vectorDA spectral_radius_terminal_inequality_transcription(
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

// Spectral radius method.
// Generalisation of [Ozaki et al. 2020]
// DOI: https://doi.org/10.2514/1.G004363
// Terminal constraints.
// PN version.
vectorDA spectral_radius_terminal_inequality_transcription(
	vectorDA const& constraints_eval, matrixdb const& Sigma_k,
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
	matrixDA Sigma = D_f*Sigma_k*D_f.transpose();
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