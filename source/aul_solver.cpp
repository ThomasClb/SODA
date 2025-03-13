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
	trajectory_split_(),
	cost_(0), violation_(0),
	d_th_order_failure_risk_(1.0),
	list_ineq_(vector<vectordb>(0)),
	list_lambda_(vector<vectordb>(0)), list_mu_(vector<vectordb>(0)),
	AUL_n_iter_(0), DDP_n_iter_(0), list_lambda_mu_() {}

AULSolver::AULSolver(
	SolverParameters const& solver_parameters,
	SpacecraftParameters const& spacecraft_parameters,
	Dynamics const& dynamics) : DDPsolver_(
		DDPSolver(solver_parameters, spacecraft_parameters, dynamics)),
	trajectory_split_(),
	cost_(0), violation_(0),
	list_ineq_(vector<vectordb>(0)),
	list_lambda_(solver_parameters.list_lambda()), list_mu_(solver_parameters.list_mu()),
	d_th_order_failure_risk_(1.0),
	AUL_n_iter_(0), DDP_n_iter_(0), list_lambda_mu_() {}

// Copy constructor
AULSolver::AULSolver(
	AULSolver const& solver) : DDPsolver_(solver.DDPsolver_),
	trajectory_split_(solver.trajectory_split_),
	cost_(solver.cost_), violation_(solver.violation_),
	d_th_order_failure_risk_(solver.d_th_order_failure_risk_),
	list_ineq_(solver.list_ineq_),
	list_lambda_(solver.list_lambda_), list_mu_(solver.list_mu_),
	AUL_n_iter_(solver.AUL_n_iter_), DDP_n_iter_(solver.DDP_n_iter_),
	list_lambda_mu_(solver.list_lambda_mu_) {}

// Destructors
AULSolver::~AULSolver() {}

// Getters
const DDPSolver AULSolver::DDPsolver() const { return DDPsolver_; }
const TrajectorySplit AULSolver::trajectory_split() const { return trajectory_split_; }
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
	double AUL_tol = solver_parameters.AUL_tol();
	vector<vectorDA> list_constraints_eval(DDPsolver_.list_deterministic_constraints_eval());
	
	// Compute diagonal blocks of Sigma
	double max_beta_d(-1e15);

	// Get diag and mean
	matrixdb Sigma_x_i, feedback_gain_i, A_i, B_i, Delta_i, R_i;
	vector<matrixdb> der_constraints;
	vectordb constraints_eval;
	double beta_d_i;
	for (size_t i=0; i<N; i++) {
		// Unpack
		Sigma_x_i = trajectory_split_.list_x()[i].Sigma();
		feedback_gain_i = trajectory_split_.list_u()[i].feedback_gain();
		der_constraints = deriv_xu(
			list_constraints_eval[i], Nx, Nu, false);
		A_i = der_constraints[0]; B_i = der_constraints[1];
		Delta_i = A_i + B_i*feedback_gain_i;
		R_i = Delta_i*Sigma_x_i*Delta_i.transpose();
		constraints_eval = list_constraints_eval[i].cons() - AUL_tol/100;
		beta_d_i = dth_order_risk_estimation(constraints_eval, get_diag_vector_(R_i));
		if (beta_d_i > max_beta_d)
			max_beta_d = beta_d_i;
	}

	// Unpack
	Sigma_x_i = trajectory_split_.list_x()[N].Sigma();
	der_constraints = deriv_x(
			list_constraints_eval[N], Nx, false);
	A_i = der_constraints[0];
	R_i = A_i*Sigma_x_i*A_i.transpose();
	constraints_eval = list_constraints_eval[N].cons() - AUL_tol/100;
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
	deque<TrajectorySplit>* const& p_list_trajectory_split,
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
	double transcription_beta = solver_parameters.transcription_beta();
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
		cout << "Target failure risk [%] : " << 100*transcription_beta << endl;
		// cout << "Intial covariance norm [-] : " << frobenius_norm_(x0.Sigma());
		cout << endl << endl << endl;
	}
	else if (verbosity < 2) {
		cout << endl;
		cout << "AUL solving - Homotopy coefficient [0, 1] : " << solver_parameters.homotopy_coefficient() << endl;
		cout << "            - Huber-loss coefficient [0, 1] : " << solver_parameters.huber_loss_coefficient() << endl;
		cout << "            - Target failure risk [%] : " << 100*transcription_beta << endl;
		// cout << "            - Intial covariance norm [-] : " << frobenius_norm_(x0.Sigma()) << endl;
		cout << "	ITERATION [-], DDP ITERATIONS [-], RUNTIME [s], FINAL MASS [kg], MAX CONSTRAINT [-], Dth-ORDER RISK [%], NLI [-]" << endl;
	}

	// Set quantiles
	double beta_star(solver_parameters.transcription_beta());
	double sum_beta_T(0);
	this->set_path_quantile(sqrt(inv_chi_2_cdf(Nineq + 1, 1 - beta_star)));
	this->set_terminal_quantile(sqrt(inv_chi_2_cdf(Ntineq + 1, 1 - beta_star)));

	// Init Robust Trajectory
	deque<TrajectorySplit> list_trajectory_split;
	deque<pair<vector<vectordb>,vector<vectordb>>> list_lambda_mu = list_lambda_mu_;
	list_lambda_mu_.clear();

	// Iterate until the list is empty.
	bool first_round(true);
	while (!p_list_trajectory_split->empty()) {

		// Get trajectory
		trajectory_split_ = p_list_trajectory_split->front();
		p_list_trajectory_split->pop_front();

		// Ouput
		cout << "Tackling split: " << trajectory_split_.splitting_history();
		cout << "Alpha: " << trajectory_split_.splitting_history().alpha() << endl;
		cout << "Beta_star: " << 100*beta_star << endl;
		cout << "sum Beta_T: " << 100*sum_beta_T << endl;

		// Init lists dual state and penalty factors lists.
		if (list_lambda_mu.size() == 0 || first_round) { // Init first round
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
			if (list_lambda_mu.size() != 0) // Remove the useless first split
				list_lambda_mu.pop_front();
		} else { // Use the lambda and mu of the previous optimised trajectory splits
			list_lambda_ = list_lambda_mu.front().first;
			list_mu_ = list_lambda_mu.front().second;
			list_lambda_mu.pop_front();
		}
		first_round = false;

		// Set dual state and penalty factors lists.
		DDPsolver_.set_list_lambda(list_lambda_);
		DDPsolver_.set_list_mu(list_mu_);

		// Init loop variables.
		bool loop = true;
		bool merged = false;
		double cost = 1e15;
		cost_ = cost;
		AUL_n_iter_ = 0;
		violation_ = cost;
		d_th_order_failure_risk_ = 1.0;
		double LOADS_tol = 1e-3; // TO DO remove make tol LOADS attribute of solver params
		unsigned int max_depth = 3; // TO DO remove make tol LOADS attribute of solver params
		double max_depth_p = 0.5;
		// max_depth_p = 2.0;
		while (loop && AUL_n_iter_ < AUL_max_iter) {
			
			// Make the trajectory split sufficiently linear
			if (trajectory_split_.list_dynamics_eval().size() == N) {

				double max_nli; 
				bool loop_nli(true);
				bool splitted(false);
				while (loop_nli) { 

					// If an addistional split si possible
					if (trajectory_split_.splitting_history().alpha() > max_depth_p) {
						// Compute NLI
						max_nli = 0;
						vectordb list_nli = nl_index(
					    	trajectory_split_.list_dynamics_eval(), trajectory_split_.list_x(),
					    	trajectory_split_.list_u(), beta_star);
						for (size_t i=0; i<list_nli.size(); i++) {
							if (list_nli[i] > LOADS_tol) {
								if (trajectory_split_.splitting_history().size() < max_depth) 
									max_nli = list_nli[i];
								splitted = true;

								cout << "split" << endl;
								cout << "	step: " << i << endl;
								cout << "	COV0 norm: " << frobenius_norm_(trajectory_split_.list_x()[0].Sigma()) << endl;
								cout << "	" << list_nli[i] << endl;

								// Find dir split
								unsigned int dir(trajectory_split_.find_splitting_direction(i, beta_star));
								cout << "	dir: " << dir << endl;

								// Split trajecectory_split_ at dir of vector 0
								pair<TrajectorySplit, TrajectorySplit> new_traj = trajectory_split_.split(
									dir, DDPsolver_);

								// Add 2 side splits to p_list_trajectory_split
								p_list_trajectory_split->push_back(new_traj.first);
								p_list_trajectory_split->push_back(new_traj.second);
								list_lambda_mu.push_back(pair<vector<vectordb>,vector<vectordb>>(list_lambda_, list_mu_));
								list_lambda_mu.push_back(pair<vector<vectordb>,vector<vectordb>>(list_lambda_, list_mu_));
								break;
							} else if (list_nli[i] > max_nli) {
								max_nli = list_nli[i];
								splitted = false;
							}
						}
						splitted = false;
					}
					loop_nli = splitted;
				}
			}

			// Solve DDP problem
			auto start = high_resolution_clock::now();
			DDPsolver_.solve(
				trajectory_split_.list_x()[0], trajectory_split_.list_u(),
				x_goal, trajectory_split_.list_dynamics_eval());
			auto stop = high_resolution_clock::now();
			auto duration_mapping = duration_cast<microseconds>(stop - start);

			// Store results
			trajectory_split_.set_list_dynamics_eval(DDPsolver_.list_dynamics_eval());
			trajectory_split_.set_list_x(DDPsolver_.list_x());
			trajectory_split_.set_list_u(DDPsolver_.list_u());
			list_ineq_ = DDPsolver_.list_ineq();
			tineq_ = DDPsolver_.tineq();
			cost_ = DDPsolver_.cost();
			d_th_order_failure_risk_ = evaluate_risk();

			// Check constraints and that the solver is not stuck
			double max_constraint_new = DDPsolver_.get_max_constraint_();
			if (abs((max_constraint_new - violation_)/ violation_) < AUL_tol)
				DDPsolver_.set_recompute_dynamics(false);
			else 
				DDPsolver_.set_recompute_dynamics(true);
			violation_ = max_constraint_new;

			// Output
			if (verbosity < 1) {
				cout << AUL_n_iter_ << " - RUNTIME [s] : "
					<< to_string(static_cast<double>(duration_mapping.count()) / 1e6) << ", "
					<< "FINAL MASS [kg] : " << DDPsolver_.list_x()[N].nominal_state()[SIZE_VECTOR] * constants.massu() << ", "
					<< "MAX CONSTRAINT [-] : " << violation_ << ", "
					<< "Dth-ORDER RISK [%] : " << 100*d_th_order_failure_risk_ << endl << endl;
			}
			else if (verbosity < 2) {
				cout << "	" << AUL_n_iter_ << ", " << DDPsolver_.n_iter()
					<< ", "	<< to_string(static_cast<double>(duration_mapping.count()) / 1e6)
					<< ", " << DDPsolver_.list_x()[N].nominal_state()[SIZE_VECTOR] * constants.massu()
					<< ", " << violation_
					<< ", " << 100*d_th_order_failure_risk_  << endl;
			}

			// Update dual state and penalities
			update_lambda_(); update_mu_();		
			list_lambda_ = DDPsolver_.solver_parameters().list_lambda();
			list_mu_ = DDPsolver_.solver_parameters().list_mu();

			// Remerge if possible
			bool merge=false;
			SplittingHistory history(trajectory_split_.splitting_history());
			if (history.size() != 0) {
				bool check_merge(history.back().second == 0);
				while (check_merge) {
					check_merge = false;

					// Check for -1 and 1 in p_list_trajectory_split
					int index_m1(-1), index_p1(-1);
					for (size_t i=0; i<p_list_trajectory_split->size(); i++) {

						// Find the neighbouring splits
						SplittingHistory history_i(p_list_trajectory_split->at(i).splitting_history());
						if (history.can_merge(history_i)) {
							if (history_i.back().second == 1)
								index_m1 = i;
							else
								index_p1 = i;
						}

						// If the splits are found
						if (index_m1 != -1 && index_p1 != -1) {
							// Check inflated NLI (/SIGMA_TILDE^2)
							TrajectorySplit trajectory_split_merge(trajectory_split_);
							trajectory_split_merge.merge(
								history.back().first, DDPsolver_);
							vectordb list_nli = nl_index(
							    trajectory_split_merge.list_dynamics_eval(), trajectory_split_merge.list_x(),
							    trajectory_split_merge.list_u(), beta_star);
							merge = true;
							for (size_t i=0; i<list_nli.size(); i++) {
								if (list_nli[i] > LOADS_tol) {
									merge = false;
									break;
								}
							}

							// Merge
							if (merge) {
								// Assign and remove -1 and 1 from p_list_trajectory_split
								trajectory_split_ = trajectory_split_merge;
								history = trajectory_split_.splitting_history();
								p_list_trajectory_split->erase(p_list_trajectory_split->begin() + max(index_m1, index_p1));
								p_list_trajectory_split->erase(p_list_trajectory_split->begin() + min(index_m1, index_p1));
								check_merge = history.back().second == 0;

								cout << "merge" << endl;
								cout << "	COV0 norm: " << frobenius_norm_(trajectory_split_.list_x()[0].Sigma()) << endl;
							}
							break;
						}
					}
				}
			}

			// Check constraints
			bool force_continue_loop = violation_ > AUL_tol || merge;

			// Stopping conditions
			bool force_stop_loop = AUL_n_iter_ > AUL_max_iter;
			loop = !force_stop_loop && force_continue_loop;
			cost = DDPsolver_.cost();
			DDP_n_iter_ += DDPsolver_.n_iter(); AUL_n_iter_++;
		}

		// Spread solution to nearby splits.
		vectordb nominal_state(trajectory_split_.list_x()[0].nominal_state());
		for (size_t i=0; i<p_list_trajectory_split->size(); i++) {
			// Unpack
			vectordb nominal_state_i(p_list_trajectory_split->at(i).list_x()[0].nominal_state());

			// Check if history_i is child of history
			double distance((nominal_state_i-nominal_state).vnorm());

			// Check if a closer parent has already be optimised in list_trajectory_split
			bool update_child(true);
			for (size_t j=0; j<list_trajectory_split.size(); j++) {
				vectordb nominal_state_j(list_trajectory_split[j].list_x()[0].nominal_state());
				if ((nominal_state_j-nominal_state).vnorm() <= distance) {
					update_child = false;
					break;
				}
			}

			// If not, update
			if (update_child) {
				SplittingHistory history_i(p_list_trajectory_split->at(i).splitting_history());				
				p_list_trajectory_split->at(i) = trajectory_split_.get_splited_trajectory(
	    			p_list_trajectory_split->at(i).list_x()[0].nominal_state(),
	    			p_list_trajectory_split->at(i).list_x()[0].Sigma(),
	    			DDPsolver_, true);
    			p_list_trajectory_split->at(i).set_splitting_history(history_i);
    			pair<vector<vectordb>,vector<vectordb>> list_lambda_mu_i(list_lambda_, list_mu_);
  				for (size_t k=0; k<N+1; k++) {
  					list_lambda_mu_i.first[k] = list_lambda_mu_i.first[k]*0.05; // perturbation lambda
  				}
    			list_lambda_mu[i] = list_lambda_mu_i;
			}
		}
		list_trajectory_split.push_back(trajectory_split_);
		list_lambda_mu_.push_back(pair<vector<vectordb>,vector<vectordb>>(list_lambda_, list_mu_));

		// Update beta_star
		double beta_T_i(min(beta_star, d_th_order_failure_risk_));
		double delta_i(beta_star-beta_T_i);
		double alpha_i(trajectory_split_.splitting_history().alpha());
		if (p_list_trajectory_split->size() > 0) {
			
			double alpha_ip1(p_list_trajectory_split->at(0).splitting_history().alpha());
			beta_star = solver_parameters.transcription_beta() + alpha_i/alpha_ip1*delta_i;
			this->set_path_quantile(sqrt(inv_chi_2_cdf(Nineq + 1, 1 - beta_star)));
			this->set_terminal_quantile(sqrt(inv_chi_2_cdf(Ntineq + 1, 1 - beta_star)));
		}
		sum_beta_T += beta_T_i*alpha_i;
	}
	*p_list_trajectory_split = list_trajectory_split;

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
