/**
	trajectory_split.cpp

	Purpose: Implementation of the TrajectorySplit class.

	@author Thomas Caleb

	@version 1.0 23/01/2025
*/

#include "trajectory_split.h"

using namespace DACE;
using namespace std;

// Default constructor
TrajectorySplit::TrajectorySplit()  = default;

// Constructor with parameters
TrajectorySplit::TrajectorySplit(
	vector<statedb> const& list_x,
    vector<controldb> const& list_u,
    SplittingHistory const& splitting_history)
    : list_x_(list_x), list_u_(list_u),
    splitting_history_(splitting_history), list_dynamics_eval_() {}

// Constructor with parameters and list dynamics eval
TrajectorySplit::TrajectorySplit(
    vector<statedb> const& list_x,
    vector<controldb> const& list_u,
    SplittingHistory const& splitting_history,
    vector<vectorDA> const& list_dynamics_eval)
    : list_x_(list_x), list_u_(list_u),
    splitting_history_(splitting_history), list_dynamics_eval_(list_dynamics_eval) {}

// Copy constructor
TrajectorySplit::TrajectorySplit(
	TrajectorySplit const& other) = default;

// Copy assignment operator
TrajectorySplit& TrajectorySplit::operator=(TrajectorySplit const& other) = default;

// Destructor
TrajectorySplit::~TrajectorySplit() = default;

// Getters
const vector<statedb> TrajectorySplit::list_x() const {return list_x_;}
const vector<controldb> TrajectorySplit::list_u() const {return list_u_;}
const SplittingHistory TrajectorySplit::splitting_history() const {return splitting_history_;}
const vector<vectorDA> TrajectorySplit::list_dynamics_eval() const {return list_dynamics_eval_;}

// Setters
void TrajectorySplit::set_list_x(vector<statedb> const& list_x) {list_x_ = list_x;}
void TrajectorySplit::set_list_u(vector<controldb> const& list_u) {list_u_ = list_u;}
void TrajectorySplit::set_splitting_history(SplittingHistory const& splitting_history) {splitting_history_ = splitting_history;}
void TrajectorySplit::set_list_dynamics_eval(vector<vectorDA> const& list_dynamics_eval) {list_dynamics_eval_ = list_dynamics_eval;}

// Find the splitting direction of highest nonlinearity.
// UNTESTED.
unsigned int TrajectorySplit::find_splitting_direction(
    size_t const& index, double const& transcription_beta) {
    // Get directional NLI
    vectordb NLI_dir = nl_index_dir(
        index, list_dynamics_eval_, list_x_,
        list_u_, transcription_beta);

    unsigned int dir = 0;
    double max_NLI_dir = 0;
    for (unsigned int i=0; i<NLI_dir.size(); i++) {
        if (NLI_dir[i] > max_NLI_dir) {
            dir = i;
            max_NLI_dir = NLI_dir[i];
        }
    }

    return dir;
}

TrajectorySplit TrajectorySplit::get_splited_trajectory(
    vectordb const& modified_x0, matrixdb const& modified_Sigma,
    DDPSolver const& DDPsolver, bool const& perturbation) {
    // Init
    size_t N = list_u_.size();
    size_t Nx(list_x_[0].nominal_state().size());
    size_t Nu(list_u_[0].nominal_control().size());
    TrajectorySplit output(*this);
    vectordb dx(modified_x0 - list_x_[0].nominal_state());
    double factor_perturbation(1.0);
    if (perturbation)
        factor_perturbation += 1e-4;

    // Change lists
    output.list_x_[0].set_nominal_state(modified_x0);
    output.list_x_[0].set_Sigma(modified_Sigma);
    vectorDA identity(vectorDA::identity(Nx + Nu));
    for (size_t i=0; i<N; i++) {
        // Fix u_i
        vectordb nominal_control(list_u_[i].nominal_control());
        matrixdb feedback_gain(list_u_[i].feedback_gain());
        vectordb du((feedback_gain*dx).extract(0, Nu - 1));
        output.list_u_[i].set_nominal_control(nominal_control*factor_perturbation + du);

        // Change x_ip1 and dynamics eval
        vectorDA delta(identity);
        for (size_t k=0; k<Nx; k++)
            delta[k] += dx[k];
        for (size_t k=0; k<Nu; k++)
            delta[k + Nx] += du[k];
        output.list_dynamics_eval_[i] = output.list_dynamics_eval_[i].eval(delta);

        // der dynamics + Sigma
        output.list_x_[i + 1] = DDPsolver.make_state(Nx, Nu, output.list_dynamics_eval_[i] , output.list_x_[i], output.list_u_[i]);
        dx = output.list_x_[i + 1].nominal_state() - list_x_[i + 1].nominal_state();
    }
    return output;
}

// Split a TrajectorySplit into 3 along a given direction.
pair<TrajectorySplit, TrajectorySplit> TrajectorySplit::split(
    unsigned int const& dir, DDPSolver const& DDPsolver) {
    // Init
    size_t N = list_u_.size();
    size_t Nx(list_x_[0].nominal_state().size());
    size_t Nu(list_u_[0].nominal_control().size());
    matrixdb navigation_error_covariance(DDPsolver.solver_parameters().navigation_error_covariance());
    vector<matrixdb> list_Sigma;

    // Get new mean and covariance
    vector<tuple<double, vectordb, matrixdb>> gmm_output = split_gmm(
    list_x_[0].nominal_state(), list_x_[0].Sigma(), dir);

    // Central split

    // Get covariance
    matrixdb Sigma_tilde(get<2>(gmm_output[0]));
    list_x_[0].set_Sigma(Sigma_tilde);
    for (size_t i=1; i<N+1; i++) {
        matrixdb der_x_i = list_x_[i].der_dynamics();
        matrixdb mat_detla_i = der_x_i.submat(0, 0, Nx - 1, Nx - 1)
            + der_x_i.submat(0, Nx, Nx - 1, Nx + Nu -1)*list_u_[i-1].feedback_gain();
        list_x_[i].set_Sigma(
            mat_detla_i*list_x_[i-1].Sigma()*mat_detla_i.transpose() + navigation_error_covariance);
    }

    // Splitting history
    splitting_history_.push_back(pair<unsigned int, int>(dir, 0));

    // Side splits

    // Init
    TrajectorySplit trajectory_split_m1(this->get_splited_trajectory(
        get<1>(gmm_output[0]), get<2>(gmm_output[0]), DDPsolver, false));
    TrajectorySplit trajectory_split_p1(this->get_splited_trajectory(
        get<1>(gmm_output[2]), get<2>(gmm_output[2]), DDPsolver, false));

    // History
    trajectory_split_m1.splitting_history_.back().second = -1;
    trajectory_split_p1.splitting_history_.back().second = 1;

   
    return pair<TrajectorySplit,TrajectorySplit>(trajectory_split_m1, trajectory_split_p1);
}

// Merge 3 TrajectorySplit into 1 along a given direction.
// The merged split is the central one.
void TrajectorySplit::merge(
    unsigned int const& dir, DDPSolver const& DDPsolver) {
    // Init
    size_t N = list_u_.size();
    size_t Nx(list_x_[0].nominal_state().size());
    size_t Nu(list_u_[0].nominal_control().size());
    matrixdb navigation_error_covariance(DDPsolver.solver_parameters().navigation_error_covariance());
    vector<matrixdb> list_Sigma;

    // Get new mean and covariance
    pair<vectordb, matrixdb> gmm_output = merge_gmm(
    list_x_[0].nominal_state(), list_x_[0].Sigma(), dir);

    // Central split

    // Get covariance
    matrixdb Sigma_tilde(gmm_output.second);
    list_x_[0].set_Sigma(Sigma_tilde);
    for (size_t i=1; i<N+1; i++) {
        matrixdb der_x_i = list_x_[i].der_dynamics();
        matrixdb mat_detla_i = der_x_i.submat(0, 0, Nx - 1, Nx - 1)
            + der_x_i.submat(0, Nx, Nx - 1, Nx + Nu -1)*list_u_[i-1].feedback_gain();
        list_x_[i].set_Sigma(
            mat_detla_i*list_x_[i-1].Sigma()*mat_detla_i.transpose() + navigation_error_covariance);
    }

    // Splitting history
    splitting_history_.pop_back();
    return;
}