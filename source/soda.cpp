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
		list_trajectory_split_(),
		AULsolver_(), PNsolver_() {}
SODA::SODA(
	SolverParameters const& solver_parameters,
	SpacecraftParameters const& spacecraft_parameters,
	Dynamics const& dynamics) :
		solver_parameters_(solver_parameters), dynamics_(dynamics), spacecraft_parameters_(spacecraft_parameters),
		list_trajectory_split_(),
		PNsolver_(), AULsolver_(solver_parameters, spacecraft_parameters, dynamics) {}

// Copy constructor
SODA::SODA(
	SODA const& solver) : 
		solver_parameters_(solver.solver_parameters_), dynamics_(solver.dynamics_),
		spacecraft_parameters_(solver.spacecraft_parameters_),  list_trajectory_split_(solver.list_trajectory_split_),
		AULsolver_(solver.AULsolver_), PNsolver_(solver.PNsolver_) {}

// Destructors
SODA::~SODA() {}

// Getters
const AULSolver SODA::AULsolver() const { return AULsolver_; }
const PNSolver SODA::PNsolver() const { return PNsolver_; }
const deque<TrajectorySplit> SODA::list_trajectory_split() const { return list_trajectory_split_; }
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

	// Init robust_trajectory
	TrajectorySplit trajectory_split_init(vector<statedb>(1, x0), list_u_init, SplittingHistory());

	// Robust Trajectory
	list_trajectory_split_ = deque<TrajectorySplit>{trajectory_split_init};

	// Run DDP
	auto start = high_resolution_clock::now();
	vectordb homotopy_sequence = solver_parameters_.homotopy_coefficient_sequence();
	vectordb huber_loss_coefficient_sequence = solver_parameters_.huber_loss_coefficient_sequence();
	AUL_n_iter_ = 0; DDP_n_iter_ = 0;
	bool loop = true;
	size_t counter = 0;
	while (loop) {
		// Prepare list_trajectory_split_
		if (fuel_optimal) {
			AULsolver_.set_homotopy_coefficient(homotopy_sequence[counter]);
			AULsolver_.set_huber_loss_coefficient(huber_loss_coefficient_sequence[counter]);
		}
		if (counter == 0) {
			TrajectorySplit trajectory_split_init(vector<statedb>(1, x0_det), list_u_init, SplittingHistory());
			list_trajectory_split_ = deque<TrajectorySplit>{trajectory_split_init};
		}
		else if (counter == 1 && robust_solving) { // Fully robust case
			AULsolver_.set_navigation_error_covariance(navigation_error_covariance);
			TrajectorySplit trajectory_split_init(
				vector<statedb>(1, x0),
				list_trajectory_split_.front().list_u(),
				SplittingHistory());
			list_trajectory_split_ = deque<TrajectorySplit>{trajectory_split_init};
		} // Else list_trajectory_split_ is already configured

		// Solve
		AULsolver_.solve(&list_trajectory_split_, x_goal);
		counter ++;
		DDP_n_iter_ += AULsolver_.DDP_n_iter();
		AUL_n_iter_ += AULsolver_.AUL_n_iter();
		loop = (counter < homotopy_sequence.size() && fuel_optimal) || (counter < 2 && robust_solving);
	}

	// PN
	auto start_inter = high_resolution_clock::now();
	PNsolver_ = PNSolver(AULsolver_);
	PN_n_iter_ = 0;
	if (pn_solving) {
		// Solve
		PNsolver_.solve(&list_trajectory_split_, x_goal);
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
		vector<statedb> list_x(list_trajectory_split_.front().list_x());
		list_x[0] = x0;
		for (size_t i=1; i<N+1; i++) {
			matrixdb der_x_i = list_x[i].der_dynamics();
			matrixdb mat_detla_i = der_x_i.submat(0, 0, Nx - 1, Nx - 1)
				+ der_x_i.submat(0, Nx, Nx - 1, Nx + Nu -1)*list_trajectory_split_.front().list_u()[i-1].feedback_gain();
			list_x[i].set_Sigma(
				mat_detla_i*list_x[i-1].Sigma()*mat_detla_i.transpose()
				+ solver_parameters_.navigation_error_covariance());
		}
		list_trajectory_split_.front().set_list_x(list_x);
	}

	// TO DO make Robust trajectory
	// Mahalanobis distance

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
	}
}
