/**
	soda.cpp

	Purpose: Implementation of the SODA solver class.

	@author Thomas Caleb

	@version 1.0 08/03/2024
*/

#include "soda.h"

using namespace DACE;
using namespace std;
using namespace std::chrono;

// Constructors
SODA::SODA() : 
		solver_parameters_(), dynamics_(), spacecraft_parameters_(),
		list_x_(), list_u_(),
		list_nli_(),
		AULsolver_(), PNsolver_() {}
SODA::SODA(
	SolverParameters const& solver_parameters,
	SpacecraftParameters const& spacecraft_parameters,
	Dynamics const& dynamics) :
		solver_parameters_(solver_parameters), dynamics_(dynamics), spacecraft_parameters_(spacecraft_parameters),
		list_x_(), list_u_(), list_nli_(),
		PNsolver_(), AULsolver_(solver_parameters, spacecraft_parameters, dynamics) {}

// Copy constructor
SODA::SODA(
	SODA const& solver) : 
		solver_parameters_(solver.solver_parameters_), dynamics_(solver.dynamics_),
		spacecraft_parameters_(solver.spacecraft_parameters_),
		list_x_(solver.list_x_), list_u_(solver.list_u_), list_nli_(solver.list_nli_),
		AULsolver_(solver.AULsolver_), PNsolver_(solver.PNsolver_) {}

// Destructors
SODA::~SODA() {}

// Getters
const AULSolver SODA::AULsolver() const { return AULsolver_; }
const PNSolver SODA::PNsolver() const { return PNsolver_; }
const vector<statedb> SODA::list_x() const { return list_x_; }
const vector<controldb> SODA::list_u() const { return list_u_; }
const vectordb SODA::list_nli() const { return list_nli_; }
const double SODA::nli() const { return nli_; }
const double SODA::cost() const { return PNsolver_.cost(); }
const double SODA::violation() const { return PNsolver_.violation(); }
const double SODA::d_th_order_failure_risk() const { return PNsolver_.d_th_order_failure_risk(); }
const double SODA::runtime() const { return runtime_; }
const double SODA::PN_runtime() const { return PN_runtime_; }
const double SODA::AUL_runtime() const { return AUL_runtime_; }
const size_t SODA::PN_n_iter() const { return PN_n_iter_; }
const size_t SODA::AUL_n_iter() const { return AUL_n_iter_; }
const size_t SODA::DDP_n_iter() const { return DDP_n_iter_; }


// Performs solving given a starting point,
// initial controls and a final state.
// Inspired from ALTRO (Julia).
// See: https://github.com/RoboticExplorationLab/Altro.jl
void SODA::solve(
	statedb const& x0,
	vector<controldb> const& list_u_init,
	statedb const& x_goal,
	bool const& robust_solving,
	bool const& fuel_optimal,
	bool const& pn_solving) {
	// Unpack
	unsigned int verbosity = solver_parameters_.verbosity();
	unsigned int N = solver_parameters_.N();
	unsigned int Nx = solver_parameters_.Nx();
	unsigned int Nu = solver_parameters_.Nu();
	unsigned int Nineq = solver_parameters_.Nineq();
	unsigned int Ntineq = solver_parameters_.Ntineq();
	double transcription_beta = solver_parameters_.transcription_beta();
	double AUL_tol = solver_parameters_.AUL_tol();
	double DDP_tol = solver_parameters_.DDP_tol();

	// Deterministic case
	statedb x0_det = x0;
	x0_det.set_Sigma(0*x0_det.Sigma());
	matrixdb navigation_error_covariance = solver_parameters_.navigation_error_covariance();
	AULsolver_.set_navigation_error_covariance(0*navigation_error_covariance);

	// Set quantiles
	AULsolver_.set_path_quantile(sqrt(inv_chi_2_cdf(Nineq + 1, 1 - transcription_beta)));
	AULsolver_.set_terminal_quantile(sqrt(inv_chi_2_cdf(Ntineq + 1, 1 - transcription_beta)));

	// Run DDP
	auto start = high_resolution_clock::now();
	vectordb homotopy_sequence = solver_parameters_.homotopy_coefficient_sequence();
	vectordb huber_loss_coefficient_sequence = solver_parameters_.huber_loss_coefficient_sequence();
	AUL_n_iter_ = 0; DDP_n_iter_ = 0;
	for (size_t i = 0; i < homotopy_sequence.size(); i++) {
		AULsolver_.set_homotopy_coefficient(homotopy_sequence[i]);
		AULsolver_.set_huber_loss_coefficient(huber_loss_coefficient_sequence[i]);

		// First iteration is always deterministic
		if (i == 0)
			AULsolver_.solve(x0_det, list_u_init, x_goal);
		else if (i != 0 && robust_solving) { // Fully robust case
			AULsolver_.set_navigation_error_covariance(navigation_error_covariance);
			AULsolver_.solve(x0, AULsolver_.list_u(), x_goal);
		}
		else // Fully deterministic case
			AULsolver_.solve(x0_det, AULsolver_.list_u(), x_goal);
		DDP_n_iter_ += AULsolver_.DDP_n_iter();
		AUL_n_iter_ += AULsolver_.AUL_n_iter();

		// For robust NRJ optimal solving
		if (!fuel_optimal && robust_solving) {
			AULsolver_.solve(x0, AULsolver_.list_u(), x_goal);
			DDP_n_iter_ += AULsolver_.DDP_n_iter();
			AUL_n_iter_ += AULsolver_.AUL_n_iter();
			break;
		} else if (!fuel_optimal)
			break;
	}
	vector<vectorDA> list_dynamic_eval = AULsolver_.DDPsolver().list_dynamic_eval();
	list_x_ = AULsolver_.list_x();
	list_u_ = AULsolver_.list_u();

	// PN
	auto start_inter = high_resolution_clock::now();
	PNsolver_ = PNSolver(AULsolver_);
	PN_n_iter_ = 0;
	if (pn_solving && AULsolver_.d_th_order_failure_risk() >= transcription_beta) {
		// Set quantiles
		AULsolver_.set_path_quantile(sqrt(inv_chi_2_cdf(Nineq + 1, 1 - transcription_beta)));
		AULsolver_.set_terminal_quantile(sqrt(inv_chi_2_cdf(Ntineq + 1, 1 - transcription_beta)));

		// Solve
		PNsolver_.solve(x_goal);
		list_x_ = PNsolver_.list_x();
		list_u_ = PNsolver_.list_u();
		list_dynamic_eval = PNsolver_.list_dynamic_eval();
		PN_n_iter_ = PNsolver_.n_iter();
	}
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
	auto duration_AUL = duration_cast<microseconds>(start_inter - start);
	auto duration_PN = duration_cast<microseconds>(stop - start_inter);
	PN_runtime_ = static_cast<double>(duration_PN.count()) / 1e6;
	AUL_runtime_ = static_cast<double>(duration_AUL.count()) / 1e6;
	runtime_ = static_cast<double>(duration.count()) / 1e6;
	

	// If det, recompute Sigmas
	if (!robust_solving) {
		list_x_[0] = x0;
		for (size_t i=1; i<N+1; i++) {
			matrixdb der_x_i = list_x_[i].der_dynamics();
			matrixdb mat_detla_i = der_x_i.submat(0, 0, Nx - 1, Nx - 1)
				+ der_x_i.submat(0, Nx, Nx - 1, Nx + Nu -1)*list_u_[i-1].feedback_gain();
			list_x_[i].set_Sigma(
				mat_detla_i*list_x_[i-1].Sigma()*mat_detla_i.transpose()
				+ solver_parameters_.navigation_error_covariance());
		}
	}

	// Get NLI
	list_nli_ = nl_index(
    	list_dynamic_eval, list_x_,
    	list_u_, solver_parameters_.transcription_beta());
	nli_ = 0.0;
	for (size_t i=0; i<N; i++) {
		if (list_nli_[i]>nli_)
			nli_ = list_nli_[i];
	}

	// Output
	if (verbosity <= 2) {
		cout << endl;
		cout << "Optimised" << endl;
		cout << "	Total runtime : " + to_string(runtime_) + "s" << endl;
		cout << "	AUL solver runtime : " + to_string(AUL_runtime_) + "s" << endl;
		cout << "	PN solver runtime : " + to_string(PN_runtime_) + "s" << endl;
		cout << "	FINAL COST [-] : " << PNsolver_.cost() << endl;
		cout << "	ERROR [-] : " << PNsolver_.violation() << endl;
		cout << "	Dth-ORDER RISK [%] : " << 100*PNsolver_.d_th_order_failure_risk() << endl;
		cout << "	Nonlinearity index [-] : " << nli_ << endl;
	}
}
