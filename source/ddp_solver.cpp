/**
	ddp_solver.cpp

	Purpose: Implementation of the DDPSolver class.

	@author Thomas Caleb

	@version 2.0 02/09/2024
*/

#include "ddp_solver.h"

using namespace DACE;
using namespace std;
using namespace std::chrono;

// Empty constructor
DDPSolver::DDPSolver() : solver_parameters_(),
	spacecraft_parameters_(Constants()), dynamics_(),
	list_x_(vector<statedb>(0)), list_u_(vector<controldb>(0)), cost_(0),
	list_ineq_(vector<vectordb>(0)),
	rho_(0.0), d_rho_(0.0) {}

// Constructor
DDPSolver::DDPSolver(
	SolverParameters const& solver_parameters,
	SpacecraftParameters const& spacecraft_parameters,
	Dynamics const& dynamics) : solver_parameters_(solver_parameters),
	spacecraft_parameters_(spacecraft_parameters), dynamics_(dynamics),
	list_x_(vector<statedb>(0)), list_u_(vector<controldb>(0)), cost_(0),
	list_ineq_(vector<vectordb>(0)),
	rho_(0.0), d_rho_(0.0), alpha_(1.0) {}

// Copy constructor
DDPSolver::DDPSolver(
	DDPSolver const& solver) : solver_parameters_(solver.solver_parameters_),
	spacecraft_parameters_(solver.spacecraft_parameters()), dynamics_(solver.dynamics()),
	list_x_(solver.list_x_), list_u_(solver.list_u_), cost_(solver.cost_),
	list_ineq_(solver.list_ineq_),
	list_dynamics_eval_(solver.list_dynamics_eval_),
	rho_(0.0), d_rho_(0.0), alpha_(1.0) {}

// Destructors
DDPSolver::~DDPSolver() {}

// Getters
const SolverParameters DDPSolver::solver_parameters() const { return solver_parameters_; }
const SpacecraftParameters DDPSolver::spacecraft_parameters() const {
	return spacecraft_parameters_;
}
const Dynamics DDPSolver::dynamics() const { return dynamics_; }
const vector<statedb> DDPSolver::list_x() const { return list_x_; }
const vector<controldb> DDPSolver::list_u() const { return list_u_; }
const vector<vectordb> DDPSolver::list_ineq() const { return list_ineq_; }
const vector<vectorDA> DDPSolver::list_dynamics_eval() const {
	return list_dynamics_eval_; }
const vector<vectorDA> DDPSolver::list_deterministic_constraints_eval() const {
	return list_deterministic_constraints_eval_; }
const vectordb DDPSolver::tineq() const { return tineq_; }
const double DDPSolver::cost() const { return cost_; }
const unsigned int DDPSolver::n_iter() const { return n_iter_; }

// Setters
void DDPSolver::set_list_lambda(vector<vectordb> const& list_lambda) {
	solver_parameters_.set_list_lambda(list_lambda);
}
void DDPSolver::set_list_mu(vector<vectordb> const& list_mu) {
	solver_parameters_.set_list_mu(list_mu);
}
void DDPSolver::set_ToF(double const& ToF) {
	solver_parameters_.set_ToF(ToF);
}
void DDPSolver::set_homotopy_coefficient(double const& homotopy_coefficient) {
	solver_parameters_.set_homotopy_coefficient(homotopy_coefficient);
}
void DDPSolver::set_huber_loss_coefficient(double const& huber_loss_coefficient) {
	solver_parameters_.set_huber_loss_coefficient(huber_loss_coefficient);
}
void DDPSolver::set_recompute_dynamics(bool const& recompute_dynamics) {
	recompute_dynamics_ = recompute_dynamics;
}
void DDPSolver::set_path_quantile(double const& path_quantile) {
	solver_parameters_.set_path_quantile(path_quantile);
}
void DDPSolver::set_terminal_quantile(double const& terminal_quantile) {
	solver_parameters_.set_terminal_quantile(terminal_quantile);
}
void DDPSolver::set_navigation_error_covariance(matrixdb const& navigation_error_covariance) {
	solver_parameters_.set_navigation_error_covariance(navigation_error_covariance);
}


// Returns the Augmented lagrangian cost-to-go:
// AUL_sc = sc + Lambda^T * c(x) + 0.5 * c(x)^T * I_mu * c(x)
vectorDA DDPSolver::get_AUL_stage_cost(
	stateDA const& x_star_DA, controlDA const& u_star_DA, size_t const& index) {
	// Unpack parameters
	unsigned int N = solver_parameters_.N();
	unsigned int Nx = solver_parameters_.Nx();
	unsigned int Nu = solver_parameters_.Nu();
	unsigned int Nineq = solver_parameters_.Nineq();
	double AUL_tol = solver_parameters_.AUL_tol();
	vectordb lambda = solver_parameters_.list_lambda()[index];
	vectordb mu = solver_parameters_.list_mu()[index];

	// Init output
	vectorDA constraints_eval;
	constraints_eval.reserve(Nineq + 1);

	// Evaluate sc
	DA sc_eval = dynamics_.stage_cost()(
		x_star_DA.nominal_state(), u_star_DA.nominal_control(),
		spacecraft_parameters_, dynamics_.constants(), solver_parameters_);

	// Constraints evaluations
	vectorDA ineq_eval = dynamics_.inequality_constraints()(
		x_star_DA.nominal_state(), u_star_DA.nominal_control(),
		spacecraft_parameters_, dynamics_.constants(), solver_parameters_);
	
	// Assign to output.
	for (size_t i = 0; i < Nineq; i++) {
		constraints_eval.push_back(ineq_eval[i]);
	}
	list_deterministic_constraints_eval_[index] = constraints_eval;
	constraints_eval.push_back(sc_eval);

	// Transcription
	if (TRANSCRIPTION_METHOD == 0) {
		constraints_eval = spectral_radius_path_transcription(
			constraints_eval, x_star_DA, u_star_DA,
			spacecraft_parameters_, dynamics_.constants(), solver_parameters_);
	}
	else if (TRANSCRIPTION_METHOD == 1) {
		constraints_eval = first_order_path_transcription(
			constraints_eval, x_star_DA, u_star_DA,
			spacecraft_parameters_, dynamics_.constants(), solver_parameters_);
	}

	// Set mu=0 for inactive constraints.
	for (size_t i = 0; i < Nineq; i++) {
		if (constraints_eval[i].cons() < 0 && lambda[i] <= 0)
			mu[i] *= 0.0;
	}

	// AUL formulation
	// Lambda^T * c(x) + 0.5 * c(x)^T * I_mu * c(x)
	sc_eval = constraints_eval[Nineq];
	for (size_t i = 0; i < Nineq; i++) {
		DA constraints_eval_i = constraints_eval[i];
		sc_eval += constraints_eval_i*(lambda[i] + (0.5 * mu[i]) *  constraints_eval_i);
	}
	constraints_eval[Nineq] = sc_eval;

	return constraints_eval;
}

// Returns the Augmented lagrangian terminal cost:
// AUL_tc = tc + Lambda^T * c(x) + 0.5 * c(x)^T * I_mu * c(x)
vectorDA DDPSolver::get_AUL_terminal_cost(
	stateDA const& x_star_DA, statedb const& x_goal) {
	// Unpack parameters
	unsigned int N = solver_parameters_.N();
	unsigned int Nx = solver_parameters_.Nx();
	unsigned int Nu = solver_parameters_.Nu();
	unsigned int Ntineq = solver_parameters_.Ntineq();
	double DDP_tol = solver_parameters_.DDP_tol();
	vectordb lambda = solver_parameters_.list_lambda()[N];
	vectordb mu = solver_parameters_.list_mu()[N];

	// Init output
	vectorDA constraints_eval;
	constraints_eval.reserve(Ntineq + 1);

	// Evaluate terminal cost
	DA tc_eval = dynamics_.terminal_cost()(
		x_star_DA.nominal_state(), x_goal.nominal_state(),
		spacecraft_parameters_, dynamics_.constants(), solver_parameters_);

	// Constraints evaluations
	vectorDA tineq_eval = dynamics_.terminal_inequality_constraints()(
		x_star_DA.nominal_state(), x_goal.nominal_state(),
		spacecraft_parameters_, dynamics_.constants(), solver_parameters_);
		
	// Assign to output.
	for (size_t i = 0; i < Ntineq; i++) {
		constraints_eval.push_back(tineq_eval[i]);
	}
	list_deterministic_constraints_eval_[N] = constraints_eval;
	constraints_eval.push_back(tc_eval);

	// Set mu=0 for inactive constraints
	for (size_t i = 0; i < Ntineq; i++) {
		if (constraints_eval[i].cons() < 0 && lambda[i] <= 0)
			mu[i] *= 0.0;
	}

	// Transcription
	if (TRANSCRIPTION_METHOD == 0) {
		constraints_eval = spectral_radius_terminal_transcription(
			constraints_eval, x_star_DA, x_goal,
			spacecraft_parameters_, dynamics_.constants(), solver_parameters_);
	}
	else if (TRANSCRIPTION_METHOD == 1) {
		constraints_eval = first_order_terminal_transcription(
			constraints_eval, x_star_DA, x_goal,
			spacecraft_parameters_, dynamics_.constants(), solver_parameters_);
	}

	// AUL formulation
	// Lambda^T * c(x) + 0.5 * c(x)^T * I_mu * c(x)
	tc_eval = constraints_eval[Ntineq];
	for (size_t i = 0; i < Ntineq; i++) {
		DA constraints_eval_i = constraints_eval[i];
		tc_eval += constraints_eval_i * (lambda[i] + (0.5 * mu[i]) * constraints_eval_i);
	}
	constraints_eval[Ntineq] = tc_eval;
	return constraints_eval;
}

// Increases the regulation to ensure Quu + rho*I
// is symeteric positive definite.
void DDPSolver::increase_regularisation_() {
	// Unpack
	vectordb rho_parameters = solver_parameters_.backward_sweep_regulation_parameters();
	double rho_init = rho_parameters[0];
	double rho_min = rho_parameters[1];
	double rho_max = rho_parameters[2];
	double rho_factor = rho_parameters[3];

	// Update d_rho_
	d_rho_ = max(d_rho_ * rho_factor, rho_factor);
	rho_ = min(rho_max, max(rho_ * d_rho_, rho_min));
}

// Increases the regulation to ensure Quu + rho*I
// is symeteric positive definite.
// With a safe guard.
void DDPSolver::increase_regularisation_(matrixdb const& Quu) {
	// Unpack
	vectordb rho_parameters = solver_parameters_.backward_sweep_regulation_parameters();
	double rho_init = rho_parameters[0];
	double rho_min = rho_parameters[1];
	double rho_max = rho_parameters[2];
	double rho_factor = rho_parameters[3];
	
	// Update d_rho_
	d_rho_ = max(d_rho_ * rho_factor, rho_factor);

	// If the reg factor is already too large
	if (rho_ >= rho_max) {
		// The small (or largest negative in absolute value) eigenvalue
		// is smaller than the norm of the symetric matrix.
		// ie |S|_1 = |D|_1 >= |eig_max|

		vectordb diag = get_diag_vector_(Quu);
		rho_ = 0;
		for (size_t i=0; i<diag.size(); i++) {
			rho_ += abs(diag[i]);
		}		
	}
	else {
		// Set rho
		rho_ = min(rho_max, max(rho_ * d_rho_, rho_min));
	}
}

// Decreases the regulation.
void DDPSolver::decrease_regularisation_() {
	// Unpack
	vectordb rho_parameters = solver_parameters_.backward_sweep_regulation_parameters();
	double rho_init = rho_parameters[0];
	double rho_min = rho_parameters[1];
	double rho_max = rho_parameters[2];
	double rho_factor = rho_parameters[3];

	// Set rho and d_rho
	d_rho_ = min(d_rho_ / rho_factor, 1/rho_factor);
	rho_ = max(min(rho_ * d_rho_, rho_max), rho_min);
}

// Computes the expected cost after backward sweep for linesearch.
double DDPSolver::expected_cost_(
	double const& alpha) {
	// Unpack parameters
	unsigned int N = solver_parameters_.N();

	// Init
	double dv_1 = 0.0;
	double dv_2 = 0.0;

	// Add expected costs
	for (size_t i = 0; i < N; i++) {
		matrixdb k_i = list_k_[N - 1 - i];
		matrixdb Qu_i = list_Qu_[N - 1 - i];
		dv_1 += (Qu_i * k_i).at(0, 0);
		dv_2 += 0.5*(Qu_i * k_i).at(0, 0);
	}
	return -1.0 * alpha * (dv_1 + alpha * dv_2);
}

// Evaluates the convergence of DDP optimisation
bool DDPSolver::evaluate_convergence_(double const& d_cost) {
	// Unpack parameters
	double tol = solver_parameters_.DDP_tol();
	unsigned int max_iter = solver_parameters_.DDP_max_iter();

	// Converged
	if ((d_cost < tol && d_cost >= 0.0 && (n_iter_ >= 3 || d_cost==0))
		|| n_iter_ >= max_iter)
		return true;

	return false;
}

// Compute the value of the max constraint.
double DDPSolver::get_max_constraint_() {
	// Unpack parameters
	SolverParameters solver_parameters = solver_parameters_;
	unsigned int N = solver_parameters.N();
	unsigned int Nineq = solver_parameters.Nineq();
	unsigned int Ntineq = solver_parameters.Ntineq();

	// Loop on all steps
	double max = -1e15;
	for (size_t i = 0; i < N; i++) {

		// Unpack
		vectordb ineq_i = list_ineq_[i];

		// Find max
		for (size_t j = 0; j < Nineq; j++) {
			double ineq_j = ineq_i[j];
			if (ineq_j > max)
				max = ineq_j;
		}
	}

	// Terminal constraints
	for (size_t j = 0; j < Ntineq; j++) {
		double tineq_j = tineq_[j];
		if (tineq_j > max)
			max = tineq_j;
	}
	return max;
}

// Returns a completed state with STM, derivatives, and covariance.
statedb DDPSolver::make_state(
	unsigned int const& Nx, unsigned int const& Nu,
	vectorDA const& x_k_DA, statedb const& x_km1, controldb const& u_km1) const {
	statedb x_k(x_k_DA.cons());
	
	// Get STM
	x_k.set_der_dynamics(x_k_DA.linear());

	// Sigma
	matrixdb der_dynamics_x_i = x_k.der_dynamics().submat(0, 0, Nx - 1, Nx - 1);
	matrixdb der_dynamics_u_i = x_k.der_dynamics().submat(0, Nx, Nx - 1, Nx + Nu -1);
	matrixdb mat_detla_i = der_dynamics_x_i + der_dynamics_u_i*u_km1.feedback_gain();
	matrixdb Sigma_x_i = mat_detla_i*x_km1.Sigma()*mat_detla_i.transpose();

	// Addition of navigation errors
	x_k.set_Sigma(Sigma_x_i + solver_parameters_.navigation_error_covariance());

	return x_k;
}

// Performs the iLQR backward sweep, that consists in the computation
// of the gains corrections.
// Inspired from ALTRO (Julia).
// See: https://github.com/RoboticExplorationLab/Altro.jl
void DDPSolver::backward_sweep_iLQR_() {
	// Unpack parameters
	unsigned int N = solver_parameters_.N();
	unsigned int Nx = solver_parameters_.Nx();
	unsigned int Nu = solver_parameters_.Nu();
	unsigned int Nineq = solver_parameters_.Nineq();
	unsigned int Ntineq = solver_parameters_.Ntineq();
	vectordb rho_parameters = solver_parameters_.backward_sweep_regulation_parameters();
	double rho_max = rho_parameters[2];

	// Evaluate terminal cost derivatives
	vector<matrixdb> der_terminal_cost = deriv_x(
		vectorDA{list_constraints_eval_[N][Ntineq] }, Nx, true);

	// Init Vx0, Vxx0
	matrixdb Vx0 = der_terminal_cost[0].transpose();
	matrixdb Vxx0 = der_terminal_cost[1];

	// Init vectors
	vectorDA vect_sc(1);

	// Backward loop
	bool success = false;
	while (!success) {
		// Init Vx, Vxx
		matrixdb Vx = Vx0;
		matrixdb Vxx = Vxx0;

		// Init gains
		list_Qu_ = vector<matrixdb>(0); list_Qu_.reserve(N);
		list_K_ = vector<matrixdb>(0); list_K_.reserve(N);
		list_k_ = vector<matrixdb>(0); list_k_.reserve(N);
		success = true;
		for (int j = N - 1; j >= 0; j--) {

			// Get derivatives
			vect_sc[0] = list_constraints_eval_[j][Nineq];
			vector<matrixdb> der_dynamic = deriv_xu(
				list_dynamics_eval_[j], Nx, Nu, false);
			vector<matrixdb> der_stage_cost = deriv_xu(
				vect_sc, Nx, Nu, true);

			// Unpack derivatives and transpose if needed
			matrixdb der_dynamic_0(der_dynamic[0]);
			matrixdb der_dynamic_1(der_dynamic[1]);
			matrixdb der_dynamic_1_t(der_dynamic_1.transpose());
			matrixdb Vx_t = Vx.transpose();

			// Get Q derivatives
			matrixdb Qx(der_stage_cost[0] + Vx_t * der_dynamic_0);
			matrixdb Qu(der_stage_cost[1] + Vx_t * der_dynamic_1);
			matrixdb Qxx(der_stage_cost[2] + der_dynamic_0.transpose() * Vxx * der_dynamic_0);
			matrixdb Quu(der_stage_cost[3] + der_dynamic_1_t * Vxx * der_dynamic_1);
			matrixdb Qux(der_stage_cost[4] + der_dynamic_1_t * Vxx * der_dynamic_0);

			// Regularisation
			for (size_t k = 0; k < Nu; k++) {
				Quu.at(k, k) += rho_;
			}

			// Check that Quu is definite positive (hence invertible)
			if (!is_def_pos_(Quu)) {
				if (rho_ == rho_max)
					increase_regularisation_(Quu); // safe version
				else
					increase_regularisation_();
				success = false; break;
			}

			// Compute Cholesky Factorisation to ease solving
			matrixdb Luu = cholesky_(Quu);

			// Compute gains
			matrixdb Qux_t = Qux.transpose();
			matrixdb K = -1.0 * solve_cholesky_(Luu, Qux);
			matrixdb k = -1.0 * solve_cholesky_(Luu, Qu.transpose());

			// Store gains and Qu
			list_Qu_.push_back(Qu);
			list_K_.push_back(K);
			list_k_.push_back(k);

			// Get Vx and Vxx for next step
			Vx = Qx.transpose() + Qux_t * k;
			Vxx = Qxx + Qux_t * K;
		}
	}
	decrease_regularisation_();
}

// Performs the DDP backward sweep, that consists in the computation
// of the gains corrections, with hessians.
void DDPSolver::backward_sweep_DDP_() {
	// Unpack parameters
	unsigned int N = solver_parameters_.N();
	unsigned int Nx = solver_parameters_.Nx();
	unsigned int Nu = solver_parameters_.Nu();
	unsigned int Nineq = solver_parameters_.Nineq();
	unsigned int Ntineq = solver_parameters_.Ntineq();
	vectordb rho_parameters = solver_parameters_.backward_sweep_regulation_parameters();
	double rho_max = rho_parameters[2];

	// Evaluate terminal cost derivatives
	vector<matrixdb> der_terminal_cost = deriv_x(
		vectorDA{list_constraints_eval_[N][Ntineq] }, Nx, true);

	// Init Vx0, Vxx0
	matrixdb Vx0 = der_terminal_cost[0].transpose();
	matrixdb Vxx0 = der_terminal_cost[1];

	// Init vectors
	vectorDA vect_sc(1);

	// Backward loop
	bool success = false;
	while (!success) {
		// Init Vx, Vxx
		matrixdb Vx = Vx0;
		matrixdb Vxx = Vxx0;

		// Init gains
		list_Qu_ = vector<matrixdb>(0); list_Qu_.reserve(N);
		list_K_ = vector<matrixdb>(0); list_K_.reserve(N);
		list_k_ = vector<matrixdb>(0); list_k_.reserve(N);
		success = true;
		for (int j = N - 1; j >= 0; j--) {

			// Get derivatives
			vect_sc[0] = list_constraints_eval_[j][Nineq];
			vector<matrixdb> der_dynamic = deriv_xu(
				list_dynamics_eval_[j], Nx, Nu, true);
			vector<matrixdb> der_cost_to_go = deriv_xu(
				vect_sc, Nx, Nu, true);

			// Unpack derivatives and transpose if needed
			matrixdb der_dynamic_0(der_dynamic[0]);
			matrixdb der_dynamic_1(der_dynamic[1]);
			matrixdb der_dynamic_1_t(der_dynamic_1.transpose());
			matrixdb Vx_t = Vx.transpose();

			// Get Q derivatives
			matrixdb Qx(der_cost_to_go[0] + Vx_t * der_dynamic_0);
			matrixdb Qu(der_cost_to_go[1] + Vx_t * der_dynamic_1);
			matrixdb Qxx(der_cost_to_go[2] + der_dynamic_0.transpose() * Vxx * der_dynamic_0);
			matrixdb Quu(der_cost_to_go[3] + der_dynamic_1_t * Vxx * der_dynamic_1);
			matrixdb Qux(der_cost_to_go[4] + der_dynamic_1_t * Vxx * der_dynamic_0);

			// Add hessians
			for (size_t k=0; k < Nx; k++) {
				matrixdb fxx_k = der_dynamic[2 + k];
				matrixdb fuu_k = der_dynamic[2 + k + Nx];
				matrixdb fux_k = der_dynamic[2 + k + 2*Nx];
				double Vx_k = Vx.at(k, 0);
				Qxx = Qxx + Vx_k * fxx_k;
				Quu = Quu + Vx_k * fuu_k;
				Qux = Qux + Vx_k * fux_k;
			}
			
			// Regularisation
			for (size_t k = 0; k < Nu; k++) {
				Quu.at(k, k) += rho_;
			}

			// Check that Quu is definite positive (hence invertible)
			if (!is_def_pos_(Quu)) {
				if (rho_ == rho_max)
					increase_regularisation_(Quu); // safe version
				else
					increase_regularisation_();
				success = false; break;
			}

			// Compute Cholesky Factorisation to ease solving
			matrixdb Luu = cholesky_(Quu);

			// Compute gains
			matrixdb Qux_t = Qux.transpose();
			matrixdb K = -1.0 * solve_cholesky_(Luu, Qux);
			matrixdb k = -1.0 * solve_cholesky_(Luu, Qu.transpose());

			// Store gains and Qu
			list_Qu_.push_back(Qu);
			list_K_.push_back(K);
			list_k_.push_back(k);

			// Get Vx and Vxx for next step
			Vx = Qx.transpose() + Qux_t * k;
			Vxx = Qxx + Qux_t * K;
		}
	}
	decrease_regularisation_();
}

// Performs the DDP forward pass, that consists in the computation
// of the new states and control after correction.
// Inspired from ALTRO (Julia).
// DA only for automatic differentiation.
// See: https://github.com/RoboticExplorationLab/Altro.jl
// DOI: 10.48550/arXiv.2502.00398 
void DDPSolver::forward_pass_(
	vector<statedb> const& list_x, vector<controldb> const& list_u, statedb const& x_goal) {
	// Unpack parameters
	double tol = solver_parameters_.DDP_tol();
	double tol_DA = solver_parameters_.AUL_tol();
	unsigned int N = solver_parameters_.N();
	unsigned int Nx = solver_parameters_.Nx();
	unsigned int Nu = solver_parameters_.Nu();
	unsigned int Nineq = solver_parameters_.Nineq();
	unsigned int Ntineq = solver_parameters_.Ntineq();
	vectordb alpha_parameters = solver_parameters_.line_search_parameters();
	double z_min = alpha_parameters[0];
	double z_max = alpha_parameters[1];
	double alpha_factor = alpha_parameters[2];
	int max_iter = alpha_parameters[3];
	matrixdb error_mat(Nx, 1);

	// Init loop variables
	bool success = false; size_t counter = 0;
	double alpha = 1.0;
	while (!success) {

		// Init lists
		vector<statedb> list_x_star(list_x_);
		vector<controldb> list_u_star(list_u_);
		vector<vectordb> list_ineq_eval;
		vectordb tineq_eval;
		vector<vectorDA> list_dynamics_eval;
		vector<vectorDA> list_constraints_eval; 
		vector<vectorDA> list_deterministic_constraints_eval_(N+1);
		vector<matrixdb> list_Sigma(N+1);
		list_Sigma[0] = list_x_[0].Sigma();
		vector<matrixdb> list_feedback_gain(N);

		// Reserve space
		list_x_star.reserve(N);
		list_u_star.reserve(N);
		list_ineq_eval.reserve(N);
		list_dynamics_eval.reserve(N);
		list_constraints_eval.reserve(N+1);

		// Rollout
		double cost = 0; double z = 1e15;
		for (size_t i = 0; i < N; i++) {

			// Get state error
			vectordb error = list_x_star[i].nominal_state() - list_x[i].nominal_state();
			error_mat.setcol(0, error);
			double error_norm = error.vnorm();

			// Retrieve gains
			matrixdb k = list_k_[N - 1 - i];
			matrixdb K = list_K_[N - 1 - i];

			// Get control correction
			vectordb correction = vectordb((alpha * k + K * error_mat).getcol(0));
			controldb u_star = list_u[i];
			u_star.set_nominal_control(u_star.nominal_control() + correction);
			u_star.set_feedback_gain(K);
			list_u_star.push_back(u_star);
			double correction_norm = correction.vnorm();

			// Concatenate two vectors to get [dx, du]
			error.reserve(Nu);
			for (size_t k = 0; k < Nu; k++) { error.push_back(correction[k]); }
			double deviation_norm = error.vnorm();

			// Get dynamics convergence radius
			vectorDA dynamic_eval_prev = list_dynamics_eval_[i];
			double dynamics_radius = convRadius(dynamic_eval_prev, tol_DA);

			// Declare variables
			stateDA x_star_DA = id_vector(list_x_star[i], 0, 0, Nx - 1);
			controlDA u_star_DA = id_vector(u_star, 0, Nx, Nx + Nu);
			vectorDA x_kp1_DA;

			// Checks if [dx, du] belongs to the convergence radius of the dynamics
			// Thus, if they can be approximated
			if (deviation_norm < dynamics_radius) {
				// Evaluate previous evaluation at [dx, du]
				vectorDA dx_star_DA = id_vector(error, 0, 0, Nx + Nu);
				x_kp1_DA = dynamic_eval_prev.eval(dx_star_DA);
			}
			else {
				// Get x_k+1 from scratchs
				x_kp1_DA = dynamics_.dynamic()(
					x_star_DA.nominal_state(), u_star_DA.nominal_control(),
					spacecraft_parameters_, dynamics_.constants(),
					solver_parameters_);
			}
			list_dynamics_eval.emplace_back(x_kp1_DA);

			// Make new state
			statedb x_kp1 = make_state(
				Nx, Nu, x_kp1_DA, list_x_star[i], u_star);
			list_x_star.emplace_back(x_kp1);
			list_Sigma[i+1] = x_kp1.Sigma();
			list_feedback_gain[i] = K;

			// Evaluate contraints
			vectorDA constraints = get_AUL_stage_cost(
					x_star_DA, u_star_DA, i);
			list_constraints_eval.emplace_back(constraints);

			// Evaluate Store constraints
			if (Nineq == 0)
				list_ineq_eval.emplace_back();
			else
				list_ineq_eval.emplace_back(constraints.extract(0, Nineq - 1).cons());

			// Update cost and append list_x_star
			cost += constraints[Nineq].cons();
		}

		// AUL terminal cost and store constraints
		statedb x_star = list_x_star[N];
		stateDA x_star_DA = id_vector(x_star, 0, 0, Nx - 1);
		vectorDA constraints = get_AUL_terminal_cost(
			x_star_DA, x_goal);
		list_constraints_eval.emplace_back(constraints);
		if (Ntineq == 0)
			tineq_eval = vectorDA(0).cons();
		else
			tineq_eval = constraints.extract(0, Ntineq - 1).cons();

		// Update cost 
		cost += constraints[Ntineq].cons();

		// Get expected cost
		double expected_cost = expected_cost_(alpha);

		// If step too small
		if (expected_cost < tol*tol) {
			list_x_ = list_x;
			list_u_ = list_u;
			increase_regularisation_();
			break;
		}
		else {
			z = (cost_ - cost) / expected_cost;

			// Check z \in interval
			if (z > z_min && z < z_max) {
				// Save all double and DA lists
				list_x_ = list_x_star;
				list_u_ = list_u_star;
				list_ineq_ = list_ineq_eval;
				tineq_ = tineq_eval;
				cost_ = cost;
				list_dynamics_eval_ = list_dynamics_eval;
				list_constraints_eval_ = list_constraints_eval;
				break;
			}
			else {
				// Decrease line search parameter
				alpha *= alpha_factor;

				// Check iteration number
				if (counter > max_iter) {
					list_x_ = list_x;
					list_u_ = list_u;
					break;
				}
			}
		}
		counter++;
	}
}

// Performs DDP solving given a starting point,
// initial controls and a final state.
// Inspired from ALTRO (Julia).
// See: https://github.com/RoboticExplorationLab/Altro.jl
void DDPSolver::solve(
	statedb const& x0,
	vector<controldb> const& list_u_init,
	statedb const& x_goal,
	vector<vectorDA> const& list_dynamics_eval) {
	// Unpack parameters
	double tol = solver_parameters_.DDP_tol();
	unsigned int max_iter = solver_parameters_.DDP_max_iter();
	unsigned int N = solver_parameters_.N();
	unsigned int Nx = solver_parameters_.Nx();
	unsigned int Nu = solver_parameters_.Nu();
	unsigned int Nineq = solver_parameters_.Nineq();
	unsigned int Ntineq = solver_parameters_.Ntineq();
	bool bs_reg = solver_parameters_.backward_sweep_regulation();
	unsigned int verbosity = solver_parameters_.verbosity();

	// Init regularisation
	rho_ = solver_parameters_.backward_sweep_regulation_parameters()[0];
	d_rho_ = 0.0;

	// Init lists
	list_constraints_eval_ = vector<vectorDA>(); 
	list_deterministic_constraints_eval_ = vector<vectorDA>(N+1); 
	list_ineq_ = vector<vectordb>();

	// Reserve
	list_ineq_.reserve(N);
	list_constraints_eval_.reserve(N+1);

	// Separte dynamic
	auto start = high_resolution_clock::now();
	if (recompute_dynamics_) {
		list_x_ = vector<statedb>();
		list_x_.reserve(N + 1);
		list_x_.push_back(x0);
		list_dynamics_eval_ = vector<vectorDA>();
		list_dynamics_eval_.reserve(N);
	} else if (list_dynamics_eval.size() == 0)
		list_dynamics_eval_ = list_dynamics_eval;

	// Make first trajectory evaluation.

	// Init cost
	cost_ = 0; n_iter_ = 0.0;
	for (size_t i = 0; i < N; i++) {

		// Make DA maps
		stateDA x_DA = id_vector(list_x_[i], 0, 0, Nx - 1);
		controlDA u_DA = id_vector(list_u_init[i], 0, Nx, Nx + Nu);
		
		// Get x_k+1
		if (recompute_dynamics_) {
			vectorDA x_kp1_DA = dynamics_.dynamic()(x_DA.nominal_state(), u_DA.nominal_control(),
				spacecraft_parameters_, dynamics_.constants(),
				solver_parameters_);
			list_dynamics_eval_.push_back(x_kp1_DA);

			// Compute the next step
			statedb x_kp1 = make_state(
				Nx, Nu, x_kp1_DA, list_x_[i], list_u_init[i]);
			list_x_.emplace_back(x_kp1);
		} else 
			list_x_[0] = x0;
		
		// Get cost to go

		// Evaluate AUL sc
		vectorDA constraints = get_AUL_stage_cost(
			x_DA, u_DA, i);
		list_constraints_eval_.emplace_back(constraints);
		if (Nineq == 0)
			list_ineq_.emplace_back(0);
		else
			list_ineq_.emplace_back(constraints.extract(0, Nineq - 1).cons());

		// Update cost
		cost_ += constraints[Nineq].cons();
	}
	list_u_ = list_u_init;

	// Evaluate AUL tc
	vectorDA constraints = get_AUL_terminal_cost(
		id_vector(list_x_[N], 0, 0, Nx - 1), x_goal);
	list_constraints_eval_.emplace_back(constraints);
	if (Ntineq == 0)
		tineq_ = vectordb(0);
	else
		tineq_ = constraints.extract(0, Ntineq - 1).cons();

	// Update cost
	cost_ += constraints[Ntineq].cons();

	// Compute costs
	double tc = constraints[Ntineq].cons();
	double sc = (cost_ - tc);

	// Output
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
	if (verbosity < 1)
		cout << "    " << 0 << " - " << cost_ << ", " << sc << ", " << tc 
		<< ", "<< to_string(static_cast<double>(duration.count()) / 1e6) << endl;
	
	// DDP solving
	vector<statedb> list_x(list_x_); vector<controldb> list_u(list_u_init);
	bool loop = true; double cost_last = 0.0;
	double d_cost(1e15);
	while (loop) {
		// Update looping variables
		n_iter_++;
		cost_last = cost_;

		// Get times
		auto start = high_resolution_clock::now();

		if (d_cost > 10*tol)
			backward_sweep_iLQR_();
		else
			backward_sweep_DDP_();

		// Forward pass init
		list_x_ = vector<statedb>(1, x0);
		list_u_ = vector<controldb>();

		// Forward pass
		forward_pass_(list_x, list_u, x_goal);

		// Compute costs
		tc = list_constraints_eval_[N][Ntineq].cons();
		sc = (cost_ - tc);

		// Evaluate convergence 
		d_cost = (cost_last - cost_)/max(abs(cost_last), abs(cost_));
		loop = !evaluate_convergence_(d_cost);

		// Update states and control
		list_x = list_x_; list_u = list_u_;

		// Output
		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(stop - start);
		if (verbosity < 1)
			cout << "    " << n_iter_ 
				<< " - " << cost_ << ", " << sc
				<< ", " << tc << ", "<< to_string(static_cast<double>(duration.count()) / 1e6) << endl;
	}
}

// Evaluates the convergence radius of a DA vector.
double convRadius(vectorDA const& x, double const& tol) {
	size_t n = x.size(); double radius = 1e15;
	for (size_t i = 0; i < n; i++) {
		double value = convRadius(x[i], tol);
		if (value < radius)
			radius = value;
	}
	return radius;
}
