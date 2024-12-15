/**
	aul_solver.cpp

	Purpose: Implementation of the AULSolver class.

	@author Thomas Caleb

	@version 2.0 02/09/2024
*/

#include "aul_solver.h"

using namespace DACE;
using namespace std;
using namespace std::chrono;

// Constructors
AULSolver::AULSolver() : DDPsolver_(DDPSolver()),
	list_x_(vector<statedb>(0)), list_u_(vector<controldb>(0)),
	cost_(0), violation_(0),
	d_th_order_failure_risk_(1.0),
	list_ineq_(vector<vectordb>(0)),
	list_lambda_(vector<vectordb>(0)), list_mu_(vector<vectordb>(0)),
	AUL_n_iter_(0), DDP_n_iter_(0) {}

AULSolver::AULSolver(
	SolverParameters const& solver_parameters,
	SpacecraftParameters const& spacecraft_parameters,
	Dynamics const& dynamics) : DDPsolver_(
		DDPSolver(solver_parameters, spacecraft_parameters, dynamics)),
	list_x_(vector<statedb>(0)), list_u_(vector<controldb>(0)),
	cost_(0), violation_(0),
	list_ineq_(vector<vectordb>(0)),
	list_lambda_(solver_parameters.list_lambda()), list_mu_(solver_parameters.list_mu()),
	d_th_order_failure_risk_(1.0),
	AUL_n_iter_(0), DDP_n_iter_(0) {}

// Copy constructor
AULSolver::AULSolver(
	AULSolver const& solver) : DDPsolver_(solver.DDPsolver_),
	list_x_(solver.list_x_), list_u_(solver.list_u_),
	cost_(solver.cost_), violation_(solver.violation_),
	d_th_order_failure_risk_(solver.d_th_order_failure_risk_),
	list_ineq_(solver.list_ineq_),
	list_lambda_(solver.list_lambda_), list_mu_(solver.list_mu_),
	AUL_n_iter_(solver.AUL_n_iter_), DDP_n_iter_(solver.DDP_n_iter_) {}

// Destructors
AULSolver::~AULSolver() {}

// Getters
const DDPSolver AULSolver::DDPsolver() const { return DDPsolver_; }
const vector<statedb> AULSolver::list_x() const { return list_x_; }
const vector<controldb> AULSolver::list_u() const { return list_u_; }
const vector<vectordb> AULSolver::list_ineq() const { return list_ineq_; }
const vectordb AULSolver::tineq() const { return tineq_; }
const double AULSolver::cost() const { return cost_; }
const double AULSolver::violation() const { return violation_; }
const double AULSolver::d_th_order_failure_risk() const { return d_th_order_failure_risk_; }
const vector<vectordb> AULSolver::list_lambda() const { return list_lambda_; }
const vector<vectordb> AULSolver::list_mu() const { return list_mu_; }
const unsigned int AULSolver::AUL_n_iter() const { return AUL_n_iter_; }
const unsigned int AULSolver::DDP_n_iter() const { return DDP_n_iter_; }

// Setters
void AULSolver::set_ToF(double const& ToF) {
	DDPsolver_.set_ToF(ToF);
}
void AULSolver::set_homotopy_coefficient(double const& homotopy_coefficient) {
	DDPsolver_.set_homotopy_coefficient(homotopy_coefficient);
}
void AULSolver::set_huber_loss_coefficient(double const& huber_loss_coefficient) {
	DDPsolver_.set_huber_loss_coefficient(huber_loss_coefficient);
}
void AULSolver::set_path_quantile(double const& path_quantile) {
	DDPsolver_.set_path_quantile(path_quantile);
}
void AULSolver::set_terminal_quantile(double const& terminal_quantile) {
	DDPsolver_.set_terminal_quantile(terminal_quantile);
}
void AULSolver::set_navigation_error_covariance(matrixdb const& navigation_error_covariance) {
	DDPsolver_.set_navigation_error_covariance(navigation_error_covariance);
}

// Update dual state in Augmented Lagrangian formulation.
void AULSolver::update_lambda_() {
	// Unpack parameters
	SolverParameters solver_parameters = DDPsolver_.solver_parameters();
	unsigned int N = solver_parameters.N();
	unsigned int Nx = solver_parameters.Nx();
	unsigned int Nu = solver_parameters.Nu();
	unsigned int Nineq = solver_parameters.Nineq();
	unsigned int Ntineq = solver_parameters.Ntineq();
	double lambda_ub = solver_parameters.lambda_parameters()[1];
	
	// Init output
	vector<vectordb> list_lambda;
	list_lambda.reserve(list_lambda.size());

	// Loop on path constraints
	for (size_t i = 0; i < N; i++) {
		
		// Unpack
		vectordb lambda = list_lambda_[i];
		vectordb mu = list_mu_[i];
		vectordb ineq = list_ineq_[i];
		
		// Iterate on constraints
		double buff = 0.0;
		for (size_t j = 0; j < Nineq; j++) {
			buff = max(0.0, lambda[j] + mu[j] * ineq[j]);
			lambda[j] = min(buff, lambda_ub);
		}

		// Assign
		list_lambda.push_back(lambda);
	}

	// Terminal contraints
	
	// Unpack
	vectordb lambda = list_lambda_[N];
	vectordb mu = list_mu_[N];
	vectordb ineq = tineq_;

	// Iterate on constraints
	double buff = 0.0;
	for (size_t j = 0; j < Ntineq; j++) {
		buff = max(0.0, lambda[j] + mu[j] * ineq[j]);
		lambda[j] = min(buff, lambda_ub);
	}

	// Assign
	list_lambda.push_back(lambda);
	DDPsolver_.set_list_lambda(list_lambda);
}

// Update penalty in Augmented Lagrangian formulation.
void AULSolver::update_mu_() {
	// Unpack parameters
	SolverParameters solver_parameters = DDPsolver_.solver_parameters();
	unsigned int N = solver_parameters.N();
	unsigned int Nx = solver_parameters.Nx();
	unsigned int Nu = solver_parameters.Nu();
	unsigned int Nineq = solver_parameters.Nineq();
	unsigned int Ntineq = solver_parameters.Ntineq();
	vectordb mu_parameters = solver_parameters.mu_parameters();
	double mu_init = mu_parameters[0]; double mu_ub = mu_parameters[1];
	double mu_factor = mu_parameters[2];

	vector<vectordb> list_mu;
	list_mu.reserve(list_mu_.size());
	double buff = 0.0;
	for (size_t i = 0; i < N; i++) {
		// Unpack
		vectordb mu = list_mu_[i];

		// Iterate on constraints
		for (size_t j = 0; j < Nineq; j++) {
			buff = mu_factor * mu[j];
			mu[j] = max(0, min(buff, mu_ub));
		}

		// Assign
		list_mu.push_back(mu);
	}

	// Unpack
	vectordb mu = list_mu_[N];

	// Iterate on constraints
	for (size_t j = 0; j < Ntineq; j++) {
		buff = mu_factor * mu[j];
		mu[j] = max(0, min(buff, mu_ub));
	}

	// Assign
	list_mu.push_back(mu);
	DDPsolver_.set_list_mu(list_mu);
}

// Computes the d-th order risk.
// DOI: WIP
double AULSolver::evaluate_risk() {
	// Unpack parameters
	SolverParameters solver_parameters = DDPsolver_.solver_parameters();
	unsigned int N = solver_parameters.N();
	unsigned int Nx = solver_parameters.Nx();
	unsigned int Nu = solver_parameters.Nu();
	unsigned int Nineq = solver_parameters.Nineq();
	unsigned int Ntineq = solver_parameters.Ntineq();
	vector<vectorDA> list_constraints_eval(DDPsolver_.list_deterministic_constraints_eval());
	
	// Compute diagonal blocks of Sigma
	double max_beta_d(-1e15);

	// Unpack
	matrixdb Sigma_x_i = list_x_[0].Sigma();
	matrixdb feedback_gain_i = list_u_[0].feedback_gain();
	vector<matrixdb> der_constraints = deriv_xu(
		list_constraints_eval[0], Nx, Nu, false);
	matrixdb A_i(der_constraints[0]); matrixdb B_i(der_constraints[1]);

	// Get first matrices
	matrixdb Delta_i = A_i + B_i*feedback_gain_i;
	matrixdb R_i = Delta_i*Sigma_x_i*Delta_i.transpose();
	vectordb constraints_eval = list_constraints_eval[0].cons();
	double beta_d_i = dth_order_risk_estimation(constraints_eval, get_diag_vector_(R_i));
	if (beta_d_i > max_beta_d)
			max_beta_d = beta_d_i;

	// Get diag and mean
	for (size_t i=1; i<N; i++) {
		// Unpack
		Sigma_x_i = list_x_[i].Sigma();
		feedback_gain_i = list_u_[i].feedback_gain();
		der_constraints = deriv_xu(
			list_constraints_eval[i], Nx, Nu, false);
		A_i = der_constraints[0]; B_i = der_constraints[1];
		Delta_i = A_i + B_i*feedback_gain_i;
		R_i = Delta_i*Sigma_x_i*Delta_i.transpose();
		constraints_eval = list_constraints_eval[i].cons();
		beta_d_i = dth_order_risk_estimation(constraints_eval, get_diag_vector_(R_i));
		if (beta_d_i > max_beta_d)
				max_beta_d = beta_d_i;
	}

	// Unpack
	Sigma_x_i = list_x_[N].Sigma();
	der_constraints = deriv_x(
			list_constraints_eval[N], Nx, false);
	A_i = der_constraints[0];
	R_i = A_i*Sigma_x_i*A_i.transpose();
	constraints_eval = list_constraints_eval[N].cons();
	beta_d_i = dth_order_risk_estimation(constraints_eval, get_diag_vector_(R_i));
	if (beta_d_i>max_beta_d)
			max_beta_d = beta_d_i;
	
	// Return
	return max_beta_d;
}

// Performs AUL solving given a starting point,
// initial controls and a final state.
// Inspired from ALTRO (Julia).
// See: https://github.com/RoboticExplorationLab/Altro.jl
void AULSolver::solve(
	statedb const& x0,
	vector<controldb> const& list_u_init,
	statedb const& x_goal) {
	// Unpack parameters
	SolverParameters solver_parameters = DDPsolver_.solver_parameters();
	unsigned int N = solver_parameters.N();
	unsigned int Nineq = solver_parameters.Nineq();
	unsigned int Ntineq = solver_parameters.Ntineq();
	double AUL_tol = solver_parameters.AUL_tol();
	int AUL_max_iter = solver_parameters.AUL_max_iter();
	vectordb mu_parameters = solver_parameters.mu_parameters();
	vectordb lambda_parameters = solver_parameters.lambda_parameters();
	unsigned int verbosity = solver_parameters.verbosity();
	unsigned int saving_iterations = solver_parameters.saving_iterations();
	Constants constants = DDPsolver_.dynamics().constants();

	// Init DDPsolver
	DDPsolver_.set_ToF(x_goal.nominal_state()[x_goal.nominal_state().size() - 1]);
	DDPsolver_.set_recompute_dynamics(true);

	// Output
	auto start_aul = high_resolution_clock::now();
	if (verbosity < 1) {
		cout << "############################################################################" << endl;
		cout << "#                                                                          #" << endl;
		cout << "#                             START AUL SOLVING                            #" << endl;
		cout << "#                                                                          #" << endl;
		cout << "############################################################################" << endl << endl;
		cout << "Homotopy coefficient [0, 1] : " << solver_parameters.homotopy_coefficient() << endl;
		cout << "Huber-loss coefficient [0, 1] : " << solver_parameters.huber_loss_coefficient() << endl;
		cout << "Target failure risk [%] : " << 100*solver_parameters.transcription_beta() << endl;
		cout << "Intial covariance norm [-] : " << frobenius_norm_(x0.Sigma());
		cout << endl << endl << endl;
	}
	else if (verbosity < 2) {
		cout << endl;
		cout << "AUL solving - Homotopy coefficient [0, 1] : " << solver_parameters.homotopy_coefficient() << endl;
		cout << "            - Huber-loss coefficient [0, 1] : " << solver_parameters.huber_loss_coefficient() << endl;
		cout << "            - Target failure risk [%] : " << 100*solver_parameters.transcription_beta() << endl;
		cout << "            - Intial covariance norm [-] : " << frobenius_norm_(x0.Sigma()) << endl;
		cout << "	ITERATION [-], DDP ITERATIONS [-], RUNTIME [s], FINAL MASS [kg], MAX CONSTRAINT [-], Dth-ORDER RISK [%], NLI [-]" << endl;
	}

	// Init lists dual state and penalty factors lists
	double lambda_0 = lambda_parameters[0];
	double mu_0 = mu_parameters[0];
	list_lambda_ = vector<vectordb>(); 
	list_mu_ = vector<vectordb>();
	list_lambda_.reserve(N + 1); list_mu_.reserve(N + 1);
	for (size_t i = 0; i < N; i++) {
		list_lambda_.emplace_back(Nineq, lambda_0);
		list_mu_.emplace_back(Nineq, mu_0);
	}
	list_lambda_.emplace_back(Ntineq, lambda_0);
	list_mu_.emplace_back(Ntineq, mu_0);

	// Set dual state and penalty factors lists
	DDPsolver_.set_list_lambda(list_lambda_);
	DDPsolver_.set_list_mu(list_mu_);
	list_u_ = list_u_init;

	// Init loop variables
	bool loop = true;
	double cost = 1e15;
	cost_ = cost;
	d_th_order_failure_risk_ = 1.0;
	AUL_n_iter_ = 0;
	size_t counter_rejected = 0;
	violation_ = cost;
	while (loop && AUL_n_iter_ < AUL_max_iter) { // TO DO
		// Solve DDP problem
		auto start = high_resolution_clock::now();
		DDPsolver_.solve(x0, list_u_, x_goal);
		auto stop = high_resolution_clock::now();
		auto duration_mapping = duration_cast<microseconds>(stop - start);

		// Store results
		list_x_ = DDPsolver_.list_x();
		list_u_ = DDPsolver_.list_u();
		list_ineq_ = DDPsolver_.list_ineq();
		tineq_ = DDPsolver_.tineq();
		cost_ = DDPsolver_.cost();
		d_th_order_failure_risk_ = evaluate_risk();

		// Check constraints and that the solver is not stuck
		double max_constraint_new = DDPsolver_.get_max_constraint_();
		if (abs((max_constraint_new - violation_)/ violation_) < AUL_tol) {
			DDPsolver_.set_recompute_dynamics(false);

			if (max_constraint_new - violation_ ==0) {
				for (size_t j=0; j<list_u_.size(); j++) {
					controldb u_j = list_u_[j];
					u_j.set_nominal_control( // Perturbation to avoid staying stuck
						u_j.nominal_control() - 0*min(violation_, AUL_tol)); 
					list_u_[j] = u_j;
				}
			}
		}
		else 
			DDPsolver_.set_recompute_dynamics(true);
		violation_ = max_constraint_new;

		// Get NLI
		vectordb list_nli = nl_index(
	    	DDPsolver_.list_dynamic_eval(), list_x_,
	    	list_u_, solver_parameters.transcription_beta());
		double nli = 0.0;
		for (size_t i=0; i<N; i++) {
			if (list_nli[i]>nli)
				nli = list_nli[i];
		}

		// Output
		if (verbosity < 1) {
			cout << AUL_n_iter_ << " - RUNTIME [s] : "
				<< to_string(static_cast<double>(duration_mapping.count()) / 1e6) << ", "
				<< "FINAL MASS [kg] : " << DDPsolver_.list_x()[N].nominal_state()[SIZE_VECTOR] * constants.massu() << ", "
				<< "MAX CONSTRAINT [-] : " << violation_ << ", "
				<< "Dth-ORDER RISK [%] : " << 100*d_th_order_failure_risk_ << ", "
				<< "NLI [-] : " << nli << endl << endl;
		}
		else if (verbosity < 2) {
			cout << "	" << AUL_n_iter_ << ", " << DDPsolver_.n_iter()
				<< ", "	<< to_string(static_cast<double>(duration_mapping.count()) / 1e6)
				<< ", " << DDPsolver_.list_x()[N].nominal_state()[SIZE_VECTOR] * constants.massu()
				<< ", " << violation_
				<< ", " << 100*d_th_order_failure_risk_
				<< ", " << nli << endl;
		}

		// Update dual state and penalities
		update_lambda_(); update_mu_();		
		list_lambda_ = DDPsolver_.solver_parameters().list_lambda();
		list_mu_ = DDPsolver_.solver_parameters().list_mu();

		// Check constraints
		bool force_continue_loop = violation_ > AUL_tol;

		// Stopping conditions
		bool force_stop_loop = AUL_n_iter_ > AUL_max_iter;
		loop = !force_stop_loop && force_continue_loop;
		cost = DDPsolver_.cost();
		DDP_n_iter_ += DDPsolver_.n_iter(); AUL_n_iter_++;
	}

	// Output
	auto stop_aul = high_resolution_clock::now();
	auto duration_aul = duration_cast<microseconds>(stop_aul - start_aul);
	if (verbosity < 1) {
		cout << "Runtime : " + to_string(static_cast<double>(duration_aul.count()) / 1e6) + "s" << endl;
	}
	else if (verbosity < 2) {
		cout << "Runtime : " + to_string(static_cast<double>(duration_aul.count()) / 1e6) + "s" << endl;
	}
}
