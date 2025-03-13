/**
	transcription.cpp

	Purpose: Implementation of the transcription methods.

	@author Thomas Caleb

	@version 2.0 26/11/2024
*/

#include "transcription.h"

using namespace DACE;
using namespace std;


vectorDA dth_order_path_inequality_transcription( // PN version
	vectorDA const& constraints_eval,
	vector<matrixdb> const& list_Sigma, vector<matrixdb> const& list_feedback_gain,
	SpacecraftParameters const& spacecraft_parameters, Constants const& constants,
	SolverParameters const& solver_parameters) {
	// Unpack
	size_t N = solver_parameters.N();
	size_t Nx = solver_parameters.Nx();
	size_t Nu = solver_parameters.Nu();
	size_t Nineq = solver_parameters.Nineq();
	size_t Ntineq = solver_parameters.Ntineq();
	double eps(solver_parameters.PN_tol());
	size_t d = constraints_eval.size();
	double beta = 1.0 - chi_2_cdf(d, sqr(solver_parameters.path_quantile()));

	// Init
	vectorDA transcribed_constraints(constraints_eval);
	vectorDA output(constraints_eval);
	vectorDA sigma(d);
	vectorDA norm_vector(d);
	size_t index_max(0);
	vector<size_t> list_index_max(N + 1);

	// Loop on all steps
	for (size_t k=0; k<N; k++) {
		// Unpack
		matrixdb Sigma_k(list_Sigma[k]);
		matrixdb feedback_gain(list_feedback_gain[k]);

		// Evaluate constraints derivatives
		vector<matrixDA> list_der = deriv_xu_DA(
			constraints_eval.extract(k*Nineq, (k+1)*Nineq - 1), Nx, Nu, false);

		// Compute the contraints covariance
		matrixDA D_f = (list_der[0] + list_der[1] * feedback_gain);
		matrixDA Sigma_kp1 = D_f*Sigma_k*D_f.transpose(); // TO DO optimize
		for (size_t i=0; i<Nineq; i++) {
			if (Sigma_kp1.at(i,i).cons() != 0.0)
				norm_vector[k*Nineq + i] = Sigma_kp1.at(i,i);
		}
	}

	// Terminal constraints

	// Unpack
	matrixdb Sigma_k(list_Sigma[N]);

	// Evaluate constraints derivatives
	vector<matrixDA> list_der = deriv_x_DA(
		constraints_eval.extract(N*Nineq, N*Nineq + Ntineq - 1), Nx, false);

	// Compute the contraints covariance
	matrixDA D_f = list_der[0];
	matrixDA Sigma = D_f*Sigma_k*D_f.transpose(); // TO DO optimize
	double max_y(0);
	for (size_t i=0; i<Ntineq; i++) {
		if (Sigma.at(i,i).cons() != 0.0)
			norm_vector[N*Nineq + i] = Sigma.at(i,i);
	}

	// Build output vector
	for (size_t k=0; k<N; k++) {
		// Extract
		vectorDA constraints_eval_k = constraints_eval;
		vectorDA norm_vector_k = norm_vector;

		// Remove terms from other steps
		for (size_t i=0; i < N*Nineq + Ntineq; i++) {
			if (i < k*Nineq || i >= (k+1)*Nineq) {
				constraints_eval_k[i] = constraints_eval_k[i].cons();
				norm_vector_k[i] = norm_vector_k[i].cons();
			}
		}

		// Get beta_d
		DA dth_order_k = dth_order_risk_estimation(
			constraints_eval_k, norm_vector_k);
		DA dth_m1(1.0 - dth_order_k/beta), dth_p1(1.0 + dth_order_k/beta); 
		for (size_t i=k*Nineq; i<(k+1)*Nineq; i++) {
			double cons_part(output[i].cons());
			if (cons_part >= eps)
				output[i] *= 1 + dth_order_k/beta;
			else if (cons_part <= -eps)
				output[i] *= 1 - dth_order_k/beta;
			else 
				output[i] += 2*eps;
		}
	}

	// Extract
	vectorDA constraints_eval_k = constraints_eval;
	vectorDA norm_vector_k = norm_vector;

	// Remove terms from other steps
	for (size_t i=0; i<N*Nineq; i++) {
		constraints_eval_k[i] = constraints_eval_k[i].cons();
		norm_vector_k[i] = norm_vector_k[i].cons();
	}

	// Get beta_d
	DA dth_order_k = dth_order_risk_estimation(
		constraints_eval_k, norm_vector_k);
	DA dth_m1(1-dth_order_k/beta), dth_p1(1+dth_order_k/beta); 

	//cout << dth_order_k << endl;

	for (size_t i=N*Nineq; i<N*Nineq + Ntineq; i++) {
		double cons_part(output[i].cons());
		if (cons_part >= eps)
			output[i] *= 1 + dth_order_k/beta;
		else if (cons_part <= -eps)
			output[i] *= 1 - dth_order_k/beta;
		else 
			output[i] += 2*eps;
	}
	cout << output.cons() << endl;
	return output;
}

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



vectorDA first_order_path_inequality_transcription( // PN version
	vectorDA const& constraints_eval,
	vector<matrixdb> const& list_Sigma, vector<matrixdb> const& list_feedback_gain,
	SpacecraftParameters const& spacecraft_parameters, Constants const& constants,
	SolverParameters const& solver_parameters) {
	// Unpack
	size_t N = solver_parameters.N();
	size_t Nx = solver_parameters.Nx();
	size_t Nu = solver_parameters.Nu();
	size_t Nineq = solver_parameters.Nineq();
	size_t Ntineq = solver_parameters.Ntineq();
	size_t d = constraints_eval.size();

	// Init
	vectorDA transcribed_constraints(constraints_eval);
	vectorDA sigma(d);
	vectorDA norm_vector(d);

	// Loop on all steps
	for (size_t k=0; k<N; k++) {
		// Unpack
		matrixdb Sigma_k(list_Sigma[k]);
		matrixdb feedback_gain(list_feedback_gain[k]);

		// Evaluate constraints derivatives
		vector<matrixDA> list_der = deriv_xu_DA(
			constraints_eval.extract(k*Nineq, (k+1)*Nineq - 1), Nx, Nu, false);

		// Compute the contraints covariance
		matrixDA D_f = (list_der[0] + list_der[1] * feedback_gain);
		matrixDA Sigma_kp1 = D_f*Sigma_k*D_f.transpose(); // TO DO optimize
		for (size_t i=0; i<Nineq; i++) {
			if (Sigma_kp1.at(i,i).cons() != 0.0) {
				norm_vector[k*Nineq + i] = sqrt(Sigma_kp1.at(i,i));
			}
		}
	}

	// Terminal constraints

	// Unpack
	matrixdb Sigma_k(list_Sigma[N]);

	// Evaluate constraints derivatives
	vector<matrixDA> list_der = deriv_x_DA(
		constraints_eval.extract(N*Nineq, N*Nineq + Ntineq - 1), Nx, false);

	// Compute the contraints covariance
	matrixDA D_f = list_der[0];
	matrixDA Sigma = D_f*Sigma_k*D_f.transpose(); // TO DO optimize
	for (size_t i=0; i<Ntineq; i++) {
		if (Sigma.at(i,i).cons() != 0.0) {
			norm_vector[N*Nineq + i] = sqrt(Sigma.at(i,i));
		}
	}

	// Make final constraints 
	return transcribed_constraints + solver_parameters.path_quantile()*norm_vector;
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