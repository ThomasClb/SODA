/**
	parameters.h

	Purpose: Implementation of the SpacecraftParameters and SolverParameters classes.

	@author Thomas Caleb

	@version 2.0 12/12/2024
*/

#ifndef DEF_PARAMETERS
#define DEF_PARAMETERS

#pragma once

#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include <dace/dace_s.h>

#include "settings.h"
#include "constants.h"
#include "stats.h"

/*

	SPACECRAFT PARAMETERS

*/
class SpacecraftParameters {

// Attributes
protected:
	Constants constants_; // Dynamics normalisation constants 
	double initial_mass_; // Spacecraft mass [MASSU]
	double dry_mass_; // Spacecraft dry mass [MASSU]
	double thrust_; // Spacecraft thrust [THRUSTU]
	double Isp_; // Spacecraft thrust Isp [TU]

// Methods
public:
	// Constructors

	// Default constructors
	SpacecraftParameters();

	// Constructors
	SpacecraftParameters(Constants const& constants);
	SpacecraftParameters(
		Constants const& constants,
		double const& initial_mass,
		double const& dry_mass,
		double const& thrust,
		double const& Isp);

	// Copy constructor
	SpacecraftParameters(
		SpacecraftParameters const& param);

	// Loader constructor
	SpacecraftParameters(std::string const& file_name);

	// Destructors
	~SpacecraftParameters();

	// Getters
	const Constants constants() const;
	const double initial_mass() const;
	const double dry_mass() const;
	const double initial_wet_mass() const;
	const double thrust() const;
	const double Isp() const;
	const double ejection_velocity() const;
	const double mass_flow() const;


	// IO operator
	friend std::ostream& operator<<(std::ostream& os, const SpacecraftParameters& param);
	friend std::istream& operator>>(std::istream& is, SpacecraftParameters& param);
	void save(std::string const& file_name) const;
	void load(std::string const& file_name);
};


/*

	SOLVER PARAMETERS

*/
class SolverParameters {

// Attributes
protected:
	unsigned int N_; // Number of steps [-]
	unsigned int Nx_; // Size of state vector including dt [-]
	unsigned int Nu_; // Size of control vector [-]
	unsigned int Nineq_; // Size of equality constraints vector [-]
	unsigned int Ntineq_; // Size of terminal equality constraints vector [-]
	double ToF_; // Time-of-flight, initialized at solving [TU]
	bool with_J2_; // With J2 dynamics [True/False]
	double stage_cost_gain_; // sc magnitude [-]
	double terminal_cost_gain_; // tc magnitude [-]
	DACE::matrixdb terminal_cost_inv_covariance_; 
	DACE::matrixdb navigation_error_covariance_; 
	double transcription_beta_; // Target failure risk [0, 1];
	double path_quantile_; // Quantile for path constraints [-];
	double terminal_quantile_; // Quantile for terminal constraints [-];
	double mass_leak_; // Mass leak term for fuel optimal optimisation [0, <<1]
	double homotopy_coefficient_; // Homotpy coefficient between the NRJ optimal and the fuel optimal optimisation [0, 1]
	double huber_loss_coefficient_; // Coefficient of the Huber loss regularisation
	DACE::vectordb homotopy_coefficient_sequence_; // Homotpy coefficient sequence for fuel optimal optimisation
	DACE::vectordb huber_loss_coefficient_sequence_; // Coefficient of the Huber loss regularisation sequence for fuel optimal optimisation
	unsigned int DDP_type_; // Choice of the DDP method 0 = classic, 1 = DA-based
	double DDP_tol_; // DDPSolver tolerance [-]
	double AUL_tol_; // AULSolver tolerance [-]
	double PN_tol_; // PNSolver tolerance [-]
	double LOADS_tol_; // LOADS tolerance [-]
	double LOADS_max_depth_; // LOADS max depth (propabilistic) [-]
	unsigned int DDP_max_iter_; // Maximum number of iterations for DDP [-]
	unsigned int AUL_max_iter_; // Maximum number of iterations for AUL [-]
	unsigned int PN_max_iter_; // Maximum number of iterations for PN [-]
	std::vector<DACE::vectordb> list_lambda_; // List of Lagrange multiplicator
	std::vector<DACE::vectordb> list_mu_; // List of penalty factors
	DACE::vectordb line_search_parameters_; // [alpha_lb, alpha_ub, alpha_factor, max_iter]
	bool backward_sweep_regulation_; // Regulation on back 
	DACE::vectordb backward_sweep_regulation_parameters_; // [rho_initial_value, rho_lb, rho_ub, rho_factor]
	DACE::vectordb lambda_parameters_; // [lambda_initial_value, lambda_ub]
	DACE::vectordb mu_parameters_; // [mu_initial_value, mu_ub, mu_factor]
	double PN_regularisation_; // PN regulation
	double PN_active_constraint_tol_; // PN active constraint tol
	double PN_cv_rate_threshold_; // PN convergence rate threshold
	double PN_alpha_; // PN intial line search parameter
	double PN_gamma_; // PN line search reduction factor
	DACE::vectordb PN_transcription_parameters_; // [eta_initial_value, eta_lb, eta_min_step, eta_max_step]
	unsigned int verbosity_; // Quantity of display data, 0=full, 1=no DDP.
	unsigned int saving_iterations_; // Quantity of iterations saved, 0=final solution, 1=final AUL and final PN, 2=AUL and final PN, 3=DDP + AUL and final PN.

// Methods
public:
	// Constructors

	// Default constructor
	SolverParameters();

	// Constructor
	SolverParameters(
		unsigned int const& N,
		unsigned int const& Nx, unsigned int const& Nu,
		unsigned int const& Nineq, unsigned int const& Ntineq,
		bool const& with_J2,
		double const& stage_cost_gain, double const& terminal_cost_gain,
		DACE::matrixdb const& terminal_cost_inv_covariance,
		DACE::matrixdb const& navigation_error_covariance,
		double const& transcription_beta_, double const& mass_leak,
		double const& homotopy_coefficient, double const& huber_loss_coefficient,
		DACE::vectordb const& homotopy_coefficient_sequence,
		DACE::vectordb const& huber_loss_coefficient_sequence,
		unsigned int const& DDP_type,
		double const& DDP_tol, double const& AUL_tol,
		double const& PN_tol, double const& LOADS_tol,
		double const& LOADS_max_depth,
		unsigned int const& DDP_max_iter, unsigned int const& AUL_max_iter,
		unsigned int const& PN_max_iter,
		DACE::vectordb const& line_search_parameters,
		bool const& backward_sweep_regulation,
		DACE::vectordb const& backward_sweep_regulation_parameters,
		DACE::vectordb const& lambda_parameters, DACE::vectordb const& mu_parameters,
		double const& PN_regularisation, double const& PN_active_constraint_tol,
		double const& PN_cv_rate_threshold, double const& PN_alpha,
		double const& PN_gamma,
		DACE::vectordb const& PN_transcription_parameters,
		unsigned int const& verbosity, unsigned int const& saving_iterations);

	// Copy constructor
	SolverParameters(SolverParameters const& param);

	// Destructors
	~SolverParameters();

	// Getters
	const unsigned int N() const;
	const unsigned int Nx() const;
	const unsigned int Nu() const;
	const unsigned int Nineq() const;
	const unsigned int Ntineq() const;
	const double ToF() const;
	const bool with_J2() const;
	const double stage_cost_gain() const;
	const double terminal_cost_gain() const;
	const DACE::matrixdb terminal_cost_inv_covariance() const;
	const DACE::matrixdb navigation_error_covariance() const;
	const double transcription_beta() const;
	const double path_quantile() const;
	const double terminal_quantile() const;
	const double mass_leak() const;
	const double homotopy_coefficient() const;
	const double huber_loss_coefficient() const;
	const DACE::vectordb homotopy_coefficient_sequence() const;
	const DACE::vectordb huber_loss_coefficient_sequence() const;
	const unsigned int DDP_type() const;
	const double DDP_tol() const;
	const double AUL_tol() const;
	const double PN_tol() const;
	const double LOADS_tol() const;
	const double LOADS_max_depth() const;
	const unsigned int DDP_max_iter() const;
	const unsigned int AUL_max_iter() const;
	const unsigned int PN_max_iter() const;
	const std::vector<DACE::vectordb> list_lambda() const;
	const std::vector<DACE::vectordb> list_mu() const;
	const DACE::vectordb line_search_parameters() const;
	const bool backward_sweep_regulation() const;
	const DACE::vectordb backward_sweep_regulation_parameters() const;
	const DACE::vectordb lambda_parameters() const;
	const DACE::vectordb mu_parameters() const;
	const double PN_regularisation() const;
	const double PN_active_constraint_tol() const;
	const double PN_cv_rate_threshold() const;
	const double PN_alpha() const;
	const double PN_gamma() const;
	const DACE::vectordb PN_transcription_parameters() const;
	const unsigned int verbosity() const;
	const unsigned int saving_iterations() const;

	// Setters
	void set_navigation_error_covariance(DACE::matrixdb const& navigation_error_covariance);
	void set_transcription_beta(double const& transcription_beta);
	void set_path_quantile(double const& path_quantile);
	void set_terminal_quantile(double const& terminal_quantile);
	void set_homotopy_coefficient(double const& homotopy_coefficient);
	void set_huber_loss_coefficient(double const& huber_loss_coefficient);
	void set_ToF(double const& ToF);
	void set_list_lambda(std::vector<DACE::vectordb> const& list_lambda);
	void set_list_mu(std::vector<DACE::vectordb> const& list_mu);

};

#endif
