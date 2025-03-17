/**
	pn_solver.cpp

	Purpose: Implementation of the PNSolver class.

	@author Thomas Caleb

	@version 2.0 02/09/2024
*/

#include "pn_solver.h"

using namespace DACE;
using namespace std;
using namespace std::chrono;

// Empty constructors
PNSolver::PNSolver() : AULsolver_(),
	solver_parameters_(), dynamics_(),
	spacecraft_parameters_(),
	cost_(0), d_th_order_failure_risk_(1.0),
	list_der_cost_(vector<vector<matrixdb>>(0)),
	list_dynamics_eval_(0), list_constraints_eval_(0),
	list_feedback_gain_(0), list_Sigma_(0),
	X_U_(0), INEQ_(0), der_INEQ_(0), correction_(0), violation_(1e15) {}

// Constructors
PNSolver::PNSolver(AULSolver const& AULsolver) : AULsolver_(AULsolver),
	cost_(AULsolver.cost()), violation_(AULsolver.violation()),
	d_th_order_failure_risk_(AULsolver_.d_th_order_failure_risk()),
	list_der_cost_(vector<vector<matrixdb>>(AULsolver.list_ineq().size() + 1)),
	list_dynamics_eval_(), list_constraints_eval_(),
	list_feedback_gain_(), list_Sigma_(),
	X_U_(), INEQ_(), der_INEQ_(), correction_() {
	// Unpack
	DDPSolver ddp_solver = AULsolver.DDPsolver();
	solver_parameters_ = ddp_solver.solver_parameters();
	spacecraft_parameters_ = ddp_solver.spacecraft_parameters();
	dynamics_ = ddp_solver.dynamics();
	unsigned int N = solver_parameters_.N();
	unsigned int Nx = solver_parameters_.Nx();
	unsigned int Nu = solver_parameters_.Nu();
	unsigned int Nineq = solver_parameters_.Nineq();
	unsigned int Ntineq = solver_parameters_.Ntineq();
	X_U_ = vectordb(N*(Nx + Nu));
	list_constraints_eval_ = vector<vectorDA>(N + 1);
	correction_ = vectordb(N*(Nx + Nu), 0);
	INEQ_ = vectordb(N*(Nx + Nineq + 1) + Ntineq + 1);
	der_INEQ_ = vector<matrixdb>(N*6 + 2);
}

// Copy constructor
PNSolver::PNSolver(
	PNSolver const& solver) : AULsolver_(solver.AULsolver_),
	solver_parameters_(solver.solver_parameters_), dynamics_(solver.dynamics_),
	spacecraft_parameters_(solver.spacecraft_parameters_),
	cost_(solver.cost_), d_th_order_failure_risk_(solver.d_th_order_failure_risk_),
	list_dynamics_eval_(solver.list_dynamics_eval_),
	list_constraints_eval_(solver.list_constraints_eval_),
	list_feedback_gain_(solver.list_feedback_gain_), list_Sigma_(solver.list_Sigma_),
	list_der_cost_(solver.list_der_cost_),
	X_U_(solver.X_U_), INEQ_(solver.INEQ_),
	der_INEQ_(solver.der_INEQ_), correction_(solver.correction_), violation_(solver.violation_) {}

// Destructors
PNSolver::~PNSolver() {}

// Getters
const AULSolver PNSolver::AULsolver() const { return AULsolver_; }
const DDPSolver PNSolver::DDPsolver() const { return AULsolver_.DDPsolver(); }
const double PNSolver::cost() const { return cost_; }
const double PNSolver::violation() const { return violation_; }
const double PNSolver::d_th_order_failure_risk() const { return d_th_order_failure_risk_; }
const vector<vectorDA> PNSolver::list_dynamics_eval() const { return list_dynamics_eval_; }
const size_t PNSolver::n_iter() const { return n_iter_; }

// Store the data of a TrajectorySplit in X_U_.
void PNSolver::set_list_x_u(TrajectorySplit const& trajectory_split) {
	// Unpack
	DDPSolver ddp_solver = AULsolver_.DDPsolver();
	solver_parameters_ = ddp_solver.solver_parameters();
	spacecraft_parameters_ = ddp_solver.spacecraft_parameters();
	dynamics_ = ddp_solver.dynamics();
	unsigned int N = solver_parameters_.N();
	unsigned int Nx = solver_parameters_.Nx();
	unsigned int Nu = solver_parameters_.Nu();
	X_U_ = vectordb(N*(Nx + Nu));
	list_dynamics_eval_ = trajectory_split.list_dynamics_eval();

	// First control
	list_Sigma_.push_back(trajectory_split.list_x()[0].Sigma());
	vectordb u_0 = trajectory_split.list_u()[0].nominal_control();
	list_feedback_gain_.push_back(trajectory_split.list_u()[0].feedback_gain());
	for (size_t j=0; j<Nu; j++)  {
		X_U_[j] = u_0[j];
	}

	// Loop
	for (size_t i=1; i<N; i++) {
		vectordb x_i = trajectory_split.list_x()[i].nominal_state();
		vectordb u_i = trajectory_split.list_u()[i].nominal_control();
		list_Sigma_.push_back(trajectory_split.list_x()[i].Sigma());
		list_feedback_gain_.push_back(trajectory_split.list_u()[i].feedback_gain());
		for (size_t j=0; j<Nx; j++)  {
			X_U_[(i - 1)*(Nx + Nu) + Nu + j] = x_i[j];
		}
		for (size_t j=0; j<Nu; j++)  {
			X_U_[i*(Nx + Nu) + j] = u_i[j];
		}
	}

	// Last state
	vectordb x_N = trajectory_split.list_x()[N].nominal_state();
	list_Sigma_.push_back(trajectory_split.list_x()[N].Sigma());
	for (size_t j=0; j<Nx; j++)  {
		X_U_[(N - 1)*(Nx + Nu) + Nu + j] = x_N[j];
	}
}

// Store the data of X_U_ in the correct TrajectorySplit.
void PNSolver::update_robust_trajectory(
	deque<TrajectorySplit>* const& p_list_trajectory_split,
	size_t const& k) {
	// Unpack
	DDPSolver ddp_solver = AULsolver_.DDPsolver();
	solver_parameters_ = ddp_solver.solver_parameters();
	spacecraft_parameters_ = ddp_solver.spacecraft_parameters();
	dynamics_ = ddp_solver.dynamics();
	unsigned int N = solver_parameters_.N();
	unsigned int Nx = solver_parameters_.Nx();
	unsigned int Nu = solver_parameters_.Nu();
	unsigned int Nineq = solver_parameters_.Nineq();
	unsigned int Ntineq = solver_parameters_.Ntineq();

	// Init
	vector<statedb> list_x{p_list_trajectory_split->at(k).list_x()[0]};
	vector<controldb> list_u;
	list_x.reserve(N); list_u.reserve(N);
	
	// Transfer data to the output lists
	controldb u(X_U_.extract( 0, Nu - 1));
	u.set_feedback_gain(list_feedback_gain_[0]);
	list_u.push_back(u);
	statedb x;
	for (size_t i=1; i< N; i++) {
		x = X_U_.extract(
			(i - 1) * (Nx + Nu) + Nu, (i - 1) * (Nx + Nu) + Nu + Nx - 1);
		x.set_Sigma(list_Sigma_[i]);
		x.set_der_dynamics(list_dynamics_eval_[i].linear());
		list_x.push_back(x);
		u = X_U_.extract(
			i * (Nx + Nu), i * (Nx + Nu) + Nu - 1);
		u.set_feedback_gain(list_feedback_gain_[i]);
		list_u.push_back(u);
	}
	x = X_U_.extract(
		(N - 1) * (Nx + Nu) + Nu, (N - 1) * (Nx + Nu) + Nu + Nx - 1);
	x.set_Sigma(list_Sigma_[N-1]);
	x.set_der_dynamics(list_dynamics_eval_[N-1].linear());
	list_x.push_back(x);

	// Assign
	p_list_trajectory_split->at(k).set_list_x(list_x);
	p_list_trajectory_split->at(k).set_list_u(list_u);
	p_list_trajectory_split->at(k).set_list_dynamics_eval(list_dynamics_eval_);
}


// Solves the optimisation problem with a projected Newton method
// Inspired from ALTRO (Julia).
// See: https://github.com/RoboticExplorationLab/Altro.jl
void PNSolver::solve(
	deque<TrajectorySplit>* const& p_list_trajectory_split,
	statedb const& x_goal) {
	// Unpack
	unsigned int N = solver_parameters_.N();
	unsigned int Nx = solver_parameters_.Nx();
	unsigned int Nu = solver_parameters_.Nu();
	unsigned int Nineq = solver_parameters_.Nineq();
	unsigned int Ntineq = solver_parameters_.Ntineq();
	unsigned int constraints_stchastic_dim((N*(Nineq + 1) + Ntineq + 1));
	size_t max_iter = solver_parameters_.PN_max_iter();
	double constraint_tol = solver_parameters_.PN_tol();
	double cv_rate_threshold = solver_parameters_.PN_cv_rate_threshold();
	unsigned int verbosity = solver_parameters_.verbosity();
	Constants constants = dynamics_.constants();

	// Output
	if (verbosity < 1) {
		cout << "############################################################################" << endl;
		cout << "#                                                                          #" << endl;
		cout << "#                              START PN SOLVING                            #" << endl;
		cout << "#                                                                          #" << endl;
		cout << "############################################################################" << endl << endl << endl;
	}
	else if (verbosity < 2) {
		cout << endl << "PN solving" << endl;
	}

	// Set quantiles
	double beta_star(solver_parameters_.transcription_beta());
	double sum_beta_T(0);
	AULsolver_.set_path_quantile(sqrt(inv_chi_2_cdf(constraints_stchastic_dim, 1 - beta_star)));

	// Loop on trajectory splits
	for (size_t k=0; k<p_list_trajectory_split->size(); k++) {
		set_list_x_u(p_list_trajectory_split->at(k));

		// Evaluate constraints
		update_constraints_(x_goal, false);
		violation_ = get_max_constraint_(INEQ_);
		double continuity_violation = get_max_continuity_constraint_(INEQ_);
		d_th_order_failure_risk_ = evaluate_risk();
		

		// Init loop
		double cv_rate = 1e15;
		double duration = 0.0;
		double violation_prev = 1e15;
		n_iter_ = 0;
		auto start = high_resolution_clock::now();
		auto stop = high_resolution_clock::now();
		for (size_t i = 0; i < max_iter; i++) {
			n_iter_ ++;

			// Output
			if (verbosity == 0) {
				stop = high_resolution_clock::now();
				auto duration = duration_cast<microseconds>(stop - start);
				string duration_str = to_string(static_cast<double>(duration.count()) / 1e6);
				start = high_resolution_clock::now();
				cout << i 
					<< " - " << duration_str
					<< " - " << X_U_[X_U_.size() - 2] * constants.massu()
					<< " - " << violation_ 
					<< " - " << 100*d_th_order_failure_risk_ << endl;
			} else if (verbosity == 1) {
				if (i % 5 == 0) {
					stop = high_resolution_clock::now();
					auto duration = duration_cast<microseconds>(stop - start);
					string duration_str = to_string(static_cast<double>(duration.count()) / 1e6);
					start = high_resolution_clock::now();
					cout << i 	
						<< " - " << duration_str
						<< " - " << X_U_[X_U_.size() - 2] * constants.massu()
						<< " - " << violation_ 
						<< " - " << 100*d_th_order_failure_risk_ << endl;
				}
			}

			// Check termination constraints
			if (violation_ < constraint_tol // If constraints are small
				|| violation_ == violation_prev) // If no progress is made
				break;

			// Update the constraints using DA
			else if (i != 0) {
				update_constraints_(x_goal, false);
				continuity_violation = get_max_continuity_constraint_(INEQ_);
				d_th_order_failure_risk_ = evaluate_risk();
			}

			// Reset corrections
			correction_ = vectordb(N*(Nx + Nu), 0);

			// Build active constraint vector and its gradient
			linearised_constraints constraints = get_linearised_constraints_();
			vector<matrixdb> block_D = get<1>(constraints);
			
			// Make S = D * D^t as a tridiagonal matrix
			sym_tridiag_matrixdb S = get_block_S_sq_(
				block_D, get<2>(constraints));

			// Compute block tridiagonal Cholesky factorisation of S.
			sym_tridiag_matrixdb L = cholesky_(S);

			// Line search loop
			cv_rate = 1e15; double violation_0 = violation_;
			violation_prev = violation_;
			for (size_t j = 0; j < 20; j++) {
				// Termination checks
				if (violation_ < constraint_tol || cv_rate < cv_rate_threshold)
					break;

				// Line search
				violation_ = line_search_(x_goal, L, block_D, get<0>(constraints), violation_0);

				// Update cv_rate
				cv_rate = log(violation_) / log(violation_0);
				violation_0 = violation_;
			}
		}

		// Output
		if (verbosity == 0) {
			stop = high_resolution_clock::now();
			auto duration = duration_cast<microseconds>(stop - start);
			string duration_str = to_string(static_cast<double>(duration.count()) / 1e6);
			start = high_resolution_clock::now();
			cout  << "OUT"
				<< " - " << duration_str
				<< " - " << X_U_[X_U_.size() - 2] * constants.massu()
				<< " - " << violation_ 
				<< " - " << 100*d_th_order_failure_risk_ << endl;
		} else if (verbosity == 1) {
			stop = high_resolution_clock::now();
			auto duration = duration_cast<microseconds>(stop - start);
			string duration_str = to_string(static_cast<double>(duration.count()) / 1e6);
			start = high_resolution_clock::now();
			cout << "OUT"
				<< " - " << duration_str
				<< " - " << X_U_[X_U_.size() - 2] * constants.massu()
				<< " - " << violation_ 
				<< " - " << 100*d_th_order_failure_risk_  << endl;
		}
		update_robust_trajectory(p_list_trajectory_split, k);

		// Update beta_star
		d_th_order_failure_risk_ = evaluate_risk();
		double beta_T_k(min(beta_star, d_th_order_failure_risk_));
		double delta_k(beta_star-beta_T_k);
		double alpha_k(p_list_trajectory_split->at(k).splitting_history().alpha());
		if (k+1 < p_list_trajectory_split->size()) {
			double alpha_kp1(p_list_trajectory_split->at(k + 1).splitting_history().alpha());
			beta_star = solver_parameters_.transcription_beta() + alpha_k/alpha_kp1*delta_k;
			AULsolver_.set_path_quantile(sqrt(inv_chi_2_cdf(constraints_stchastic_dim, 1 - beta_star)));
		}
		sum_beta_T += alpha_k*beta_T_k;
		cout << "	beta_star : " << 100*beta_star << endl;
		cout << "	sum beta_T : " << 100*sum_beta_T << endl;
	}
	return;
}

// Iterative line search for PN.
// Inspired from ALTRO (Julia).
// See: https://github.com/RoboticExplorationLab/Altro.jl
double PNSolver::line_search_(
	statedb const& x_goal,
	sym_tridiag_matrixdb const& L,
	vector<matrixdb> const& block_D,
	vectordb const& d_0,
	double const& violation_0) {
	// Unpack
	size_t max_iter = solver_parameters_.PN_max_iter();
	double alpha = solver_parameters_.PN_alpha();
	double gamma = solver_parameters_.PN_gamma();
	unsigned int N = solver_parameters_.N();
	unsigned int Nx = solver_parameters_.Nx();
	unsigned int Nu = solver_parameters_.Nu();
	matrixdb mat_z(d_0.size(), 1);

	// Solve S * z = d <=> L * L^t * z = d using
	// block tridiagonal Cholesky factorisation.
	vectordb z = solve_cholesky_(L, d_0);
	mat_z.setcol(0, z);
	vectordb correction(N * (Nu + Nx));

	// Compute correction
	size_t counter_z = 0; size_t counter_corr = 0;
	for (size_t j = 0; j < block_D.size(); j++) {
		vectordb z_j;
		matrixdb D_j = block_D[j];
		size_t N_j = D_j.nrows();

		if (N_j != 0)
			z_j = z.extract(counter_z, N_j + counter_z - 1);

		matrixdb mat_z_j(z_j.size(), 1); mat_z_j.setcol(0, z_j);
		vectordb corr_j = (mat_z_j.transpose() * D_j).getrow(0);

		// Assign
		for (size_t k = 0; k < corr_j.size(); k++) {
			correction[counter_corr + k] += alpha * corr_j[k];
		}
		counter_z += N_j;
		counter_corr += corr_j.size() - Nx;
	}

	// Init loop
	double violation(violation_0);
	for (size_t i = 0; i < 20; i++) {
		// Compute X_U
		vectordb X_U = X_U_ - correction;

		// Evaluate constraints
		vectordb INEQ = update_constraints_double_( 
			x_goal, X_U, correction);

		// Get the max constraint
		violation = get_max_constraint_(INEQ);

		// Check exit
		if (violation <= violation_0) {
			// Update states corrections and constraints
			X_U_ = X_U;
			correction_ -= correction;
			INEQ_ = INEQ;
			return violation;
		}
		else {
			// Restore constraints
			violation = violation_0;
			correction *= gamma;
		}
	}
	return violation_0;
}

// Evaluates the dth order risk.
// DOI: WIP.
// TO DO refine.
double PNSolver::evaluate_risk() {

	// Unpack parameters
	unsigned int N = solver_parameters_.N();
	unsigned int Nx = solver_parameters_.Nx();
	unsigned int Nu = solver_parameters_.Nu();
	
	// Get diag and mean
	vectordb mean, diag_Sigma;
	matrixdb A_i, B_i, Delta_i, Sigma;
	vector<matrixdb> der_constraints;
	vectordb constraints_eval;
	for (size_t i=0; i<N; i++) {
		// Unpack
		der_constraints = deriv_xu(
			list_constraints_eval_[i], Nx, Nu, false);
		Delta_i = der_constraints[0] + der_constraints[1]*list_feedback_gain_[i];
		Sigma = Delta_i*list_Sigma_[i]*Delta_i.transpose();
		constraints_eval = list_constraints_eval_[i].cons();
		for (size_t j=0; j<constraints_eval.size(); j++) {
			mean.push_back(constraints_eval[j]);
			diag_Sigma.push_back(Sigma.at(j,j));
		}
	}

	// Unpack
	der_constraints = deriv_x(
			list_constraints_eval_[N], Nx, false);
	Sigma = der_constraints[0]*list_Sigma_[N]*der_constraints[0].transpose();
	constraints_eval = list_constraints_eval_[N].cons();
	for (size_t j=0; j<constraints_eval.size(); j++) {
		mean.push_back(constraints_eval[j]);
		diag_Sigma.push_back(Sigma.at(j,j));
	}

	// Return
	return dth_order_risk_estimation(mean, diag_Sigma);
}

// Computes the maximum constraints given eq in ineq constraints
double PNSolver::get_max_constraint_(
	DACE::vectordb const& INEQ) {
	// Unpack parameters
	unsigned int N = solver_parameters_.N();
	unsigned int Nx = solver_parameters_.Nx();
	unsigned int Nineq = solver_parameters_.Nineq();
	unsigned int Ntineq = solver_parameters_.Ntineq();

	// Loop on all steps
	double maximum = -1e15;
	for (size_t i = 0; i < N; i++) {

		// Find maximum
		for (size_t j = 0; j < Nx; j++) {
			maximum = max(maximum, abs(INEQ[i*(Nineq + 1 + Nx) + j]));
		}
		for (size_t j = 0; j < Nineq + 1; j++) {
			maximum = max(maximum, INEQ[i*(Nineq + 1 + Nx) + Nx + j]);
		}
	}

	// Terminal constraints
	for (size_t j = 0; j < Ntineq + 1; j++) {
		maximum = max(maximum, INEQ[N*(Nineq + 1 + Nx) + j]);
	}

	return maximum;
}

// Computes the maximum constraints given eq in ineq constraints
double PNSolver::get_max_continuity_constraint_(
	DACE::vectordb const& INEQ) {
	// Unpack parameters
	unsigned int N = solver_parameters_.N();
	unsigned int Nx = solver_parameters_.Nx();
	unsigned int Nineq = solver_parameters_.Nineq();

	// Loop on all steps
	double maximum = -1e15;
	for (size_t i = 0; i < N; i++) {
		for (size_t j = 0; j < Nx; j++) {
			maximum = max(maximum, abs(INEQ[i*(Nineq + 1 + Nx) + j]));
		}
	}

	return maximum;
}

// Computes the new constraints given states and controls.
// Can recompute all DA maps and the dynamics.
// Updates the derivatives.
void PNSolver::update_constraints_(
	statedb const& x_goal,
	bool const& force_DA) {
	// Unpack
	Constants constants = dynamics_.constants();
	double tol = solver_parameters_.PN_active_constraint_tol();
	unsigned int N = solver_parameters_.N();
	unsigned int Nx = solver_parameters_.Nx();
	unsigned int Nu = solver_parameters_.Nu();
	unsigned int Nineq = solver_parameters_.Nineq();
	unsigned int Ntineq = solver_parameters_.Ntineq();

	// Make id matrix
	matrixdb id_Nx(Nx, Nx, 0.0), null_Nx(Nineq + 1, Nx, 0.0);
	for (size_t i = 0; i < Nx; i++) { id_Nx.at(i, i) = -1.0; }

	// Update path constraints

	// Loop on all steps
	vectorDA x_DA, u_DA, dx_u_DA;
	vectordb x_kp1, dx_u(Nx + Nu);
	vectorDA INEQ_stochastic(N*(Nineq + 1) + Ntineq + 1, 0);
	for (size_t i = 0; i < N; i++) {

		// Get dx_u and x_kp1
		if (!force_DA) {
			for (size_t j=0; j<Nu; j++) {dx_u[Nx + j] = correction_[i*(Nu + Nx) + j];}
			if (i == 0) 
				for (size_t j=0; j<Nx; j++) {dx_u[j] = 0.0;}
			else
				for (size_t j=0; j<Nx; j++) {dx_u[j] = correction_[(i - 1)*(Nu + Nx) + Nu + j];}
			dx_u_DA = id_vector(dx_u, 0, 0, Nx + Nu);
		}
		x_kp1 = X_U_.extract(
				Nu + i*(Nx + Nu),
				Nu + Nx - 1 + i*(Nx + Nu));

		// Get DA x, u
		if (i == 0)
			x_DA = id_vector(X_U_.extract(Nu, Nu + Nx - 1),
				0, 0, Nx - 1);
		else 
			x_DA = id_vector(
				X_U_.extract(
					Nu + (i - 1)*(Nx + Nu),
					Nu + Nx - 1 + (i - 1)*(Nx + Nu)),
				0, 0, Nx - 1);
		u_DA = id_vector(
			X_U_.extract(
				i*(Nx + Nu),
				Nu - 1 + i*(Nx + Nu)),
			0, Nx, Nu + Nx);

		// Continuity constraints
		vectorDA x_kp1_eval;
		if (force_DA) {
			x_kp1_eval = dynamics_.dynamic()(
				x_DA, u_DA,
				spacecraft_parameters_, constants, solver_parameters_);
		} else {
			// Get conv radius
			vectorDA dynamics_eval = list_dynamics_eval_[i];
			double radius = convRadius(dynamics_eval, tol);
			double norm = dx_u.vnorm();

			// Compute the next step using the previous map
			if (norm < radius) {
				x_kp1_eval = dynamics_eval.eval(dx_u_DA);	
			}

			// Compute from scratch	
			else {
				x_kp1_eval = dynamics_.dynamic()(
					x_DA, u_DA,
					spacecraft_parameters_, constants, solver_parameters_);
			}
		}
		list_dynamics_eval_[i] = x_kp1_eval;

		// Get Sigma
		matrixdb der_x = x_kp1_eval.linear();
		matrixdb mat_detla_i = der_x.submat(0, 0, Nx - 1, Nx - 1)
			+ der_x.submat(0, Nx, Nx - 1, Nx + Nu -1)*list_feedback_gain_[i];
		list_Sigma_[i+1] = mat_detla_i*list_Sigma_[i]*mat_detla_i.transpose()
			+ solver_parameters_.navigation_error_covariance();

		// Constraints evaluations
		vectorDA constraints_eval = dynamics_.inequality_constraints()(
			x_DA, u_DA, spacecraft_parameters_, constants, solver_parameters_);

		// Assign
		list_constraints_eval_[i] = constraints_eval;
		for (size_t k = 0; k < Nineq; k++) {
			INEQ_stochastic[i*(Nineq + 1) + k] = constraints_eval[k];}

		// Add continuity constraints
		vectorDA eq_eval = x_kp1_eval - x_kp1;

		// Get derivatives
		vector<matrixdb> der_eq = deriv_xu(
			eq_eval, Nx, Nu, false);
		
		// Dependency with kp1
		der_INEQ_[6*i] = der_eq[0];
		der_INEQ_[6*i + 1] = der_eq[1];
		der_INEQ_[6*i + 2] = id_Nx;
		der_INEQ_[6*i + 5] = null_Nx;

		// Assign
		for (size_t k = 0; k < Nx; k++) {
			INEQ_[i*(Nx + Nineq + 1) + k] = eq_eval[k].cons();}
	}

	// Update terminal constraints

	// Get DA x, u
	x_DA = id_vector(
		X_U_.extract(
			Nu + (N - 1)*(Nx + Nu),
			Nu + Nx - 1 + (N - 1)*(Nx + Nu)),
		0, 0, Nx - 1);

	// Constraints evaluations
	vectorDA constraints_eval = dynamics_.terminal_inequality_constraints()(
		x_DA, x_goal.nominal_state(), spacecraft_parameters_, constants, solver_parameters_);

	// Assign
	list_constraints_eval_[N] = constraints_eval;
	for (size_t k = 0; k < Ntineq; k++) {
		INEQ_stochastic[N*(Nineq  + 1) + k] = constraints_eval[k];}

	// Transcription
	if (TRANSCRIPTION_METHOD == 0) {
		// TO DO
	}
	/*
	else if (TRANSCRIPTION_METHOD == 1) {
		INEQ_stochastic = first_order_inequality_transcription(
			INEQ_stochastic,
			list_Sigma_, list_feedback_gain_,
			spacecraft_parameters_, constants,
			solver_parameters_);
		
	} */else if (TRANSCRIPTION_METHOD == 1) {
		/**/
		INEQ_stochastic = dth_order_inequality_transcription(
			INEQ_stochastic,
			list_Sigma_, list_feedback_gain_,
			spacecraft_parameters_, constants,
			solver_parameters_);
		
	}

	// Get derivatives

	// Path
	for (size_t i=0; i<N; i++) {
		vector<matrixdb> der_ineq = deriv_xu(
				INEQ_stochastic.extract(i*(Nineq + 1), (i+1)*(Nineq + 1) - 1), Nx, Nu, false);
		der_INEQ_[6*i + 3] = der_ineq[0];
		der_INEQ_[6*i + 4] = der_ineq[1];
	}

	// Terminal
	vector<matrixdb> der_tineq = deriv_xu(
		INEQ_stochastic.extract(N*(Nineq + 1), N*(Nineq  + 1) + Ntineq + 1 - 1), Nx, Nu, false);
	der_INEQ_[6*N] = der_tineq[0];
	der_INEQ_[6*N + 1] = der_tineq[1];

	// Assign
	for (size_t i=0; i<N; i++) {
		for (size_t k = 0; k < Nineq + 1; k++) {
			INEQ_[i*(Nx + Nineq + 1) + Nx + k] = INEQ_stochastic[i*(Nineq  + 1) + k].cons();}
	}
	for (size_t k = 0; k < Ntineq + 1; k++) {
		INEQ_[N*(Nx + Nineq + 1) + k] = INEQ_stochastic[N*(Nineq  + 1) + k].cons();}
}

// Computes the new constraints given states and controls.
// Uses the DA mapping of the dynamics.
// Faster than update_constraints_.
vectordb PNSolver::update_constraints_double_(
	statedb const& x_goal,
	vectordb const& X_U,
	vectordb const& correction) {
	// Unpack
	Constants constants = dynamics_.constants();
	double tol = solver_parameters_.PN_active_constraint_tol();
	unsigned int N = solver_parameters_.N();
	unsigned int Nx = solver_parameters_.Nx();
	unsigned int Nu = solver_parameters_.Nu();
	unsigned int Nineq = solver_parameters_.Nineq();
	unsigned int Ntineq = solver_parameters_.Ntineq();
	vectordb corr = correction_ - correction;

	// Init lists
	vectordb INEQ(N*(Nx + Nineq + 1) + Ntineq + 1);
	vectorDA INEQ_stochastic(N*(Nineq  + 1) + Ntineq + 1, 0);

	// Update path constraints

	// Loop on all steps
	vectordb x(Nx), xp1(Nx), u(Nu), dx_u(Nx + Nx);
	vectorDA x_DA(Nx), u_DA(Nu);
	for (size_t i = 0; i < N; i++) {

		// Get x, u, xp1, du, dx
		u = X_U.extract(i * (Nu + Nx), i * (Nu + Nx) + Nu - 1);
		for (size_t j=0; j<Nu; j++) {
			dx_u[Nx + j] = corr[i*(Nu + Nx) + j];
		}
		if (i == 0) {
			x = X_U_.extract(Nu, Nu + Nx - 1); // Never changes
			for (size_t j=0; j<Nx; j++) {
				dx_u[j] = 0.0;
			}
		}
		else {
			x = X_U.extract((i - 1) * (Nu + Nx) + Nu, (i - 1) * (Nu + Nx) + Nu + Nx - 1);
			for (size_t j=0; j<Nx; j++) {
				dx_u[j] = corr[(i - 1)*(Nu + Nx) + Nu + j];
			}
		}
		xp1 = X_U.extract(i * (Nu + Nx) + Nu, i * (Nu + Nx) + Nu + Nx - 1);

		// Continuity constraints

		// Get conv radius
		vectorDA dynamics_eval = list_dynamics_eval_[i];
		double radius = convRadius(dynamics_eval, tol);
		double norm = dx_u.vnorm();

		// Compiute the next step
		vectordb x_kp1_eval;
		if (norm < radius)
			x_kp1_eval = dynamics_eval.eval(dx_u) - xp1;
		else
			x_kp1_eval = dynamics_.dynamic_db()(
				x, u, spacecraft_parameters_, constants, solver_parameters_) - xp1;

		// Get DA x, u
		x_DA = id_vector(x, 0, 0, Nx - 1);
		u_DA = id_vector(u, 0, Nx, Nu + Nx);

		// Constraints evaluations
		vectorDA constraints_eval = dynamics_.inequality_constraints()(
			x_DA, u_DA, spacecraft_parameters_, constants, solver_parameters_);

		// Assign
		for (size_t k = 0; k < Nx; k++) {
			INEQ[i*(Nx + Nineq + 1) + k] = x_kp1_eval[k];}
		for (size_t k = 0; k < Nineq; k++) {
			INEQ_stochastic[i*(Nineq + 1) + k] = constraints_eval[k];}
	}

	// Update terminal constraints

	// Get x_f
	x = X_U.extract((N - 1) * (Nu + Nx) + Nu, (N - 1) * (Nu + Nx) + Nu + Nx - 1);
	x_DA = id_vector(x, 0, 0, Nx - 1);

	// Constraints evaluations
	vectorDA constraints_eval = dynamics_.terminal_inequality_constraints()(
		x_DA, x_goal.nominal_state(), spacecraft_parameters_, constants, solver_parameters_);

	// Assign
	for (size_t k = 0; k < Ntineq; k++) {
		INEQ_stochastic[N*(Nineq + 1) + k] = constraints_eval[k];}

	// Transcription
	if (TRANSCRIPTION_METHOD == 0) {
		// TO DO
	}
	/*
	else if (TRANSCRIPTION_METHOD == 1) {
		INEQ_stochastic = first_order_inequality_transcription(
			INEQ_stochastic,
			list_Sigma_, list_feedback_gain_,
			spacecraft_parameters_, constants,
			solver_parameters_);
		
	} */
	else if (TRANSCRIPTION_METHOD == 1) {
		/**/
		INEQ_stochastic = dth_order_inequality_transcription(
			INEQ_stochastic,
			list_Sigma_, list_feedback_gain_,
			spacecraft_parameters_, constants,
			solver_parameters_);
	}

	// Assign
	for (size_t i = 0; i < N; i++) {
		for (size_t k = 0; k < Nineq + 1; k++) {
			INEQ[i*(Nx + Nineq + 1) + Nx + k] = INEQ_stochastic[i*(Nineq + 1) + k].cons();}
	}
	for (size_t k = 0; k < Ntineq + 1; k++) {
		INEQ[N*(Nx + Nineq + 1) + k] = INEQ_stochastic[N*(Nineq + 1) + k].cons();}
	return INEQ;
}

// Returns the vector of active constraints and their gradients
// - first it the active constraints vector.
// - second is a pair with the list of gradients of constraints first.
// - third is the list of active constraints.
linearised_constraints PNSolver::get_linearised_constraints_() {
	// Unpack
	unsigned int N = solver_parameters_.N();
	unsigned int Nx = solver_parameters_.Nx();
	unsigned int Nu = solver_parameters_.Nu();
	unsigned int Nineq = solver_parameters_.Nineq();
	unsigned int Ntineq = solver_parameters_.Ntineq();
	double active_constraint_tol = solver_parameters_.PN_active_constraint_tol();

	// Init the output
	vectordb d; vector<matrixdb> block_D;
	vector<vector<size_t>> list_active_constraint_index(2*N + 1);
	d.reserve(N*(Nineq + 1 + Nx) + Ntineq + 1); block_D.reserve(N + 1);

	// Path constraints
	for (size_t i = 0; i < N; i++) {
		// Init lists
		vector<size_t> list_active_ineq_index, list_active_c_index;
		list_active_ineq_index.reserve(Nineq + 1);
		list_active_c_index.reserve(Nx); 

		// Inequality constraints
		for (size_t j = 0; j < Nineq + 1; j++) {
			double ineq_i_j = INEQ_[i*(Nx + Nineq + 1) + Nx + j];
			if (ineq_i_j > active_constraint_tol) {
				// Assign to d
				d.push_back(ineq_i_j);

				// Save index
				list_active_ineq_index.push_back(j);
			}
		}

		// Continuity constraints
		for (size_t j = 0; j < Nx; j++) {
			double eq_i_j = INEQ_[i*(Nx + Nineq + 1) + j];
			if (abs(eq_i_j) > active_constraint_tol) {
				// Assign to d
				d.push_back(eq_i_j);

				// Save index
				list_active_c_index.push_back(j);
			}
		}

		// Build gradient
		size_t Nineq_active(list_active_ineq_index.size());
		size_t Nc_active(list_active_c_index.size());
		matrixdb D_i(
			Nineq_active + Nc_active,
			(2 * Nx + Nu));
		matrixdb der_eq_x = der_INEQ_[6*i];
		matrixdb der_eq_u = der_INEQ_[6*i + 1];
		matrixdb der_eq_x_kp1 = der_INEQ_[6*i + 2];
		matrixdb der_ineq_x = der_INEQ_[6*i + 3];
		matrixdb der_ineq_u = der_INEQ_[6*i + 4];
		matrixdb der_ineq_x_kp1 = der_INEQ_[6*i + 5];

		// Inequality constraints
		for (size_t j = 0; j < list_active_ineq_index.size(); j++) {
			// Get derivatives
			size_t active_ineq_index_j = list_active_ineq_index[j];
			vectordb dx(der_ineq_x.getrow(active_ineq_index_j));
			vectordb du(der_ineq_u.getrow(active_ineq_index_j));
			vectordb dxkp1(der_ineq_x_kp1.getrow(active_ineq_index_j));

			// Assign to D_i
			size_t index_j = j;
			for (size_t k = 0; k < Nx; k++) {
				D_i.at(index_j, k) = dx[k];
				D_i.at(index_j, k + Nx + Nu) = dxkp1[k];
			}
			for (size_t k = 0; k < Nu; k++) {
				D_i.at(index_j, k + Nx) = du[k];
			}
		}

		// Continuity constraints
		for (size_t j = 0; j < list_active_c_index.size(); j++) {
			// Get derivatives
			size_t active_c_index_j = list_active_c_index[j];
			vectordb dx(der_eq_x.getrow(active_c_index_j));
			vectordb du(der_eq_u.getrow(active_c_index_j));
			vectordb dxkp1(der_eq_x_kp1.getrow(active_c_index_j));

			// Assign to D_i
			size_t index_j = j + Nineq_active;
			for (size_t k = 0; k < Nx; k++) {
				D_i.at(index_j, k) = dx[k];
				D_i.at(index_j, k + Nx + Nu) = dxkp1[k];
			}
			for (size_t k = 0; k < Nu; k++) {
				D_i.at(index_j, k + Nx) = du[k];
			}
		}

		// x0 is not a variable
		if (i == 0) {
			if (D_i.nrows() != 0)
				D_i = D_i.submat(0, Nx, D_i.nrows() - 1, 2 * Nx + Nu - 1);
			else
				D_i = matrixdb(0, Nx + Nu);
		}

		// Assign
		block_D.push_back(D_i);
		list_active_constraint_index[2*i + 0] = list_active_ineq_index;
		list_active_constraint_index[2*i + 1] = list_active_c_index;
	}

	// Terminal constraints
	vector<size_t> list_active_tineq_index;
	list_active_tineq_index.reserve(Ntineq + 1);

	// Inequality constraints
	for (size_t j = 0; j < Ntineq + 1; j++) {
		double tineq_j = INEQ_[N*(Nx + Nineq + 1) + j];
		if (tineq_j > active_constraint_tol) {
			// Assign to d
			d.push_back(tineq_j);

			// Save index
			list_active_tineq_index.push_back(j);
		}
	}

	// Build gradient
	size_t Ntineq_active(list_active_tineq_index.size());
	matrixdb D_i(
		Ntineq_active,
		Nx);

	// Inequality constraints
	matrixdb der_tineq_x = der_INEQ_[6*N];
	for (size_t j = 0; j < list_active_tineq_index.size(); j++) {
		vectordb dx(der_tineq_x.getrow(list_active_tineq_index[j]));
		D_i.setrow(j, dx);
	}

	// Assign
	block_D.push_back(D_i);
	list_active_constraint_index[2*N + 0] = list_active_tineq_index;

	// Return d & D
	return linearised_constraints(d, block_D, list_active_constraint_index);
}

// Return the matrix S = D_a * D_a^t 
// Where D_a is the gradient of the active linearized constraints.
// Using tridiagonal symetric block computation.
sym_tridiag_matrixdb PNSolver::get_block_S_sq_(
	vector<matrixdb> const& block_Delta,
	vector<vector<size_t>> const& list_active_index) {
	// Unpack
	unsigned int N = solver_parameters_.N();
	unsigned int Nx = solver_parameters_.Nx();
	unsigned int Nu = solver_parameters_.Nu();

	//  Init lists
	vector<matrixdb> list_diag, list_subdiag;
	list_diag.reserve(block_Delta.size());
	list_subdiag.reserve(block_Delta.size() - 1);

	// Path elements
	for (size_t i = 0; i < N; i++) {

		// Unpack
		matrixdb Delta_i = block_Delta[i];
		size_t NINEQ_i(list_active_index[2*i].size());
		size_t Nc_i(list_active_index[2*i + 1].size());

		// Unpack contraints gradient blocks
		matrixdb A_i(NINEQ_i, Nx), C_i(Nc_i, Nx);
		matrixdb B_i(NINEQ_i, Nu), D_i(Nc_i, Nu);

		if (NINEQ_i!= 0) {
			if (i == 0)
				B_i = Delta_i.submat(
					0, 0,
					NINEQ_i- 1, Nu - 1);
			else {
				A_i = Delta_i.submat(
					0, 0,
					NINEQ_i - 1, Nx - 1);
				B_i = Delta_i.submat(
					0, Nx,
					NINEQ_i - 1, Nx + Nu - 1);
			}
		}
		if (Nc_i != 0) {
			if (i == 0)
				D_i = Delta_i.submat(
					NINEQ_i, 0,
					NINEQ_i + Nc_i - 1, Nu - 1);
			else {
				C_i = Delta_i.submat(
					NINEQ_i, 0,
					NINEQ_i + Nc_i - 1, Nx - 1);
				D_i = Delta_i.submat(
					NINEQ_i, Nx,
					NINEQ_i + Nc_i - 1, Nx + Nu - 1);
			}
		}
		matrixdb B_i_t = B_i.transpose(); matrixdb D_i_t = D_i.transpose();
		matrixdb A_i_t = A_i.transpose(); matrixdb C_i_t = C_i.transpose();

		// Block formulation
		matrixdb R_i, S_i, T_i;
		matrixdb W_i, Y_i;
		if (i == 0) {
			// Diagonal
			R_i = B_i * B_i_t;
			S_i = D_i * D_i_t;
			T_i = D_i * B_i_t;
		}
		else {
			// Diagonal
			R_i = B_i * B_i_t + A_i * A_i_t;
			S_i = D_i * D_i_t + C_i * C_i_t;
			T_i = D_i * B_i_t + C_i * A_i_t;

			// Tridiag
			W_i = -1.0 * (A_i); Y_i = -1.0 * (C_i);
		}
		size_t size_S = S_i.ncols();
		for (size_t j = 0; j < size_S; j++) { S_i.at(j, j) += 1.0; }

		// Make diagonal matrix 
		size_t size_R = R_i.ncols();
		size_t size_diag_im1 = 0;
		matrixdb diag_i(size_R + size_S);
		for (size_t j = 0; j < size_R; j++) {

			// Assign
			for (size_t k = 0; k < size_R; k++) {
				diag_i.at(j, k) = R_i.at(j, k);
			}
			for (size_t k = 0; k < size_S; k++) {
				diag_i.at(j, k + size_R) = T_i.at(k, j);
			}
		}
		for (size_t j = 0; j < size_S; j++) {

			// Assign
			for (size_t k = 0; k < size_R; k++) {
				diag_i.at(size_R + j, k) = T_i.at(j, k);
			}
			for (size_t k = 0; k < size_S; k++) {
				diag_i.at(size_R + j, k + size_R) = S_i.at(j, k);
			}
		}

		// Subdiag
		matrixdb subdiag_i;
		if (i != 0) {
			size_diag_im1 = list_diag[i - 1].ncols();
			size_t size_diag_i = diag_i.ncols();
			vector<size_t> list_active_index_c_im1 = list_active_index[2*(i - 1) + 1];
			size_t Nc_im1 = list_active_index_c_im1.size();
			subdiag_i = matrixdb(size_diag_i, size_diag_im1);
			if (size_diag_i != 0 && Nc_im1 != 0 && size_diag_im1 != 0) {
				for (size_t j = 0; j < Nc_im1; j++) {
					// Get cols
					size_t active_index_c_im1_j = list_active_index_c_im1[j];
					vectordb col_W_j = W_i.getcol(active_index_c_im1_j);
					col_W_j.reserve(Y_i.nrows());
					for (size_t k = 0; k < Y_i.nrows(); k++) {
						col_W_j.push_back(Y_i.at(k, active_index_c_im1_j));
					}

					// Assign
					subdiag_i.setcol(size_diag_im1 - Nc_im1 + j, col_W_j);
				}
			}
		}
	
		// Assign
		list_diag.push_back(diag_i);
		if (i != 0)
			list_subdiag.push_back(subdiag_i);
	}

	// Last elements
	size_t i = N;
	matrixdb Delta_i = block_Delta[N];
	size_t Ntineq(list_active_index[2*N].size());
	matrixdb A_i(Ntineq, Nx, 0.0);
	if (Ntineq != 0) {
		A_i = Delta_i.submat(
			0, 0,
			Ntineq - 1, Nx - 1);
	}
	matrixdb A_i_t = A_i.transpose();

	// Block formulation

	// Diagonal
	matrixdb R_i = A_i * A_i_t;

	// Make diagonal matrix 
	size_t size_R = R_i.ncols();
	matrixdb diag_i(R_i);

	// Tridiag
	size_t size_diag_im1 = list_diag[N - 1].ncols();
	size_t size_diag_i = diag_i.ncols();
	vector<size_t> list_active_index_c_im1 = list_active_index[2*(N - 1) + 1];
	size_t Nc_im1 = list_active_index_c_im1.size();
	matrixdb subdiag_i(size_diag_i, size_diag_im1);
	matrixdb W_i = -1.0 * A_i;
	if (size_diag_i != 0 && Nc_im1 != 0 && size_diag_im1 != 0) {
		for (size_t j = 0; j < Nc_im1; j++) {
			// Get cols
			vectordb col_W_j = W_i.getcol(list_active_index_c_im1[j]);

			// Assign
			subdiag_i.setcol(size_diag_im1 - Nc_im1 + j, col_W_j);
		}
	}

	// Assign
	list_diag.push_back(diag_i);
	list_subdiag.push_back(subdiag_i);
	
	return sym_tridiag_matrixdb(list_diag, list_subdiag);
}
