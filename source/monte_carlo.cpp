/**b
	monte_carlo.cpp

	Purpose: Implementation of the Monte-Carlo validation methods.

	@author Thomas Caleb

	@version 1.0 26/09/2024
*/

#include "monte_carlo.h"

using namespace DACE;
using namespace std;
using namespace std::chrono;

// Returns a sample of a multivariate normal distribution.
matrixdb generate_normal_sample(
	size_t const& size_sample, vectordb const& mean, matrixdb const& covariance) {
	// Make generator
	default_random_engine generator;
	// generator.seed(123456789); // For debug purpose.
  	normal_distribution<double> distribution(0.0);

  	// Make sample

  	// Generate the normalised sample
  	size_t size_vector(mean.size());
  	matrixdb output(size_vector, size_sample);
  	for (size_t i=0; i<size_vector; i++) {
  		for (size_t j=0; j<size_sample; j++) {
  			output.at(i, j) = distribution(generator);
  		}
  	}

  	// Scaling
  	matrixdb cholesky_Sigma;
  	if (diagonal_error_(covariance) == 0)
		cholesky_Sigma = make_diag_matrix_(sqrt(get_diag_vector_(covariance)));
  	else
  		cholesky_Sigma = cholesky_(covariance);
  	output = cholesky_Sigma*output; // Scale

  	// Offset
  	for (size_t i=0; i<size_vector; i++) {
  		double mean_i = mean[i];
  		for (size_t j=0; j<size_sample; j++) {
  			output.at(i, j) += mean_i; // Offset
  		}
  	}

  	return output;
}

// Propagates a trajectory and returns
// the remaining fuel and the max constraint.
vectordb propagate_trajectory(
	vectordb const& x_0_sample,
	statedb const& x_goal,
	vector<statedb> const& list_x,
	vector<controldb> const& list_u,
	matrixdb const& navigation_error_sample,
	SolverParameters const& solver_parameters,
	SpacecraftParameters const& spacecraft_parameters,
	Dynamics const& dynamics, int const& test_case_id,
	matrixdb* const& p_mat_state,
	matrixdb* const& p_mat_control,
	matrixdb* const& p_mat_path_constraints,
	matrixdb* const& p_mat_terminal_constraints,
	bool const& save_trajectory) {
	// Unpack parameters
	unsigned int N = solver_parameters.N();
	unsigned int Nx = solver_parameters.Nx();
	unsigned int Nu = solver_parameters.Nu();
	unsigned int Nineq = solver_parameters.Nineq();
	unsigned int Ntineq = solver_parameters.Ntineq();
	Constants constants(dynamics.constants());
	double lu = constants.lu();
	double massu = constants.massu();
	double tu = constants.tu();
	double thrustu = constants.thrustu();
	double vu = constants.vu();

	// Init vector lists
	vector<vectordb> list_x_sample{x_0_sample};
	vector<vectordb> list_constraints;
	vectordb list_cost;
	double cost(0.0);

	// Save value if needed
	if (save_trajectory) {
		p_mat_state->setrow(0, x_0_sample);
	}

	// Loop on all steps
	for (size_t i = 0; i < N; i++) {

		// Get x, u
		vectordb x_sample = list_x_sample[i];
		controldb u = list_u[i];
		statedb x = list_x[i];
		vectordb dx = x_sample - x.nominal_state();
		vectordb du = (u.feedback_gain()*dx).extract(0, Nu-1);
		vectordb u_sample = u.nominal_control() + du;

		// Propagate dynamics
		vectordb x_kp1_sample = dynamics.dynamic_db()(
			x_sample, u_sample,
			spacecraft_parameters, constants, solver_parameters);

		// Add nav error
		x_kp1_sample += vectordb(navigation_error_sample.getcol(i));
		list_x_sample.push_back(x_kp1_sample);

		// Evaluate sc
		double sc = dynamics.stage_cost_db()(
			x_sample, u_sample,
			spacecraft_parameters, constants, solver_parameters);
		list_cost.push_back(sc);
		cost += sc;

		// Constraints evaluations
		vectordb constraints = dynamics.inequality_constraints_db()(
			x_sample, u_sample,
			spacecraft_parameters, constants, solver_parameters);
		list_constraints.push_back(constraints);

		// Save value if needed
		if (save_trajectory) {
			p_mat_control->setrow(i, u_sample);
			p_mat_state->setrow(i+1, x_kp1_sample);
			p_mat_path_constraints->setrow(i, constraints);
		}
	}

	// Get x
	vectordb x_sample = list_x_sample[N];

	// Evaluate terminal cost
	double tc = dynamics.terminal_cost_db()(
		x_sample, x_goal.nominal_state(),
		spacecraft_parameters, constants, solver_parameters);
	cost += tc;

	// Constraints evaluations
	vectordb constraints = dynamics.terminal_inequality_constraints_db()(
		x_sample, x_goal.nominal_state(),
		spacecraft_parameters, constants, solver_parameters);
	list_constraints.push_back(constraints);

	// Save value if needed
	if (save_trajectory) {
		p_mat_terminal_constraints->setrow(0, constraints);
	}
	
	// Check constraints
	double max_constraint = -1e15;
	for (size_t i = 0; i < N; i++) {
		constraints = list_constraints[i];
		for (size_t j=0; j<Nineq; j++) {
			if (constraints[j]>max_constraint)
				max_constraint = constraints[j];
		}
	}
	constraints = list_constraints[N];
	for (size_t j=0; j<Ntineq; j++) {
		if (constraints[j]>max_constraint)
			max_constraint = constraints[j];
	}

	// Output
	vectordb output(2);
	if (test_case_id != 0) {
		output[0] = list_x_sample[N][SIZE_VECTOR]*massu;
		output[1] = max_constraint;
	}
	else {
		output[0] = cost;
		output[1] = 0;
	}

	return output;
}

// Tests a robust command law with MC.
vector<vector<matrixdb>>  test_trajectory(
	vector<statedb> const& list_x, vector<controldb> const& list_u,
	statedb const& x_goal, size_t const& size_sample,
	SODA const& solver, int const& test_case_id,
	bool const& robust_optimisation,
	SolverParameters const& solver_parameters,
	SpacecraftParameters const& spacecraft_parameters,
	Dynamics const& dynamics, size_t const& size_sample_saved_) {
	// Unpack
	unsigned int N = solver_parameters.N();
	unsigned int Nx = solver_parameters.Nx();
	unsigned int Nu = solver_parameters.Nu();
	unsigned int Nineq = solver_parameters.Nineq();
	unsigned int Ntineq = solver_parameters.Ntineq();
	Constants constants(dynamics.constants());
	double lu = constants.lu();
	double massu = constants.massu();
	double tu = constants.tu();
	double thrustu = constants.thrustu();
	double vu = constants.vu();
	int verbosity = solver_parameters.verbosity();

	// Init output variables
	size_t size_sample_saved = size_sample_saved_;
	if (size_sample_saved_ > size_sample) {
		size_sample_saved = size_sample;
	}
	vector<matrixdb> list_mat_state, list_mat_control, list_mat_path_constraints, list_mat_terminal_constraints;
	matrixdb mat_state_(N+1, Nx), mat_control_(N, Nu), mat_path_constraints_(N, Nineq), mat_terminal_constraints_(1, Ntineq);
	if (size_sample_saved > 0) {
		list_mat_state.reserve(size_sample_saved);
		list_mat_control.reserve(size_sample_saved);
		list_mat_path_constraints.reserve(size_sample_saved);
		list_mat_terminal_constraints.reserve(size_sample_saved);
	}

	// First output
	if (verbosity < 3) {
		cout << endl;
		cout << "Monte-Carlo" << endl;
		cout << "	Sample size: " << size_sample << endl;
		cout << "	Est Run time [s]: " << size_sample*0.01 << endl;
	}

	// Get nominal trajectory
	auto start = high_resolution_clock::now();
	vectordb nominal_results = propagate_trajectory(
			list_x[0].nominal_state(), x_goal,
			list_x, list_u,
			matrixdb(list_x[0].nominal_state().size(), N, 0.0),
			solver_parameters, spacecraft_parameters, dynamics,
			test_case_id,
			&mat_state_, &mat_control_, 
			&mat_path_constraints_, &mat_terminal_constraints_,
			false);
	
	// Get initial conditions sample
	matrixdb ic_sample = generate_normal_sample(
		size_sample, list_x[0].nominal_state(), list_x[0].Sigma());

	// Get navigation error sample
	matrixdb nav_sample = generate_normal_sample(
		size_sample*N, 0.0*list_x[0].nominal_state(), solver_parameters.navigation_error_covariance());

	// Compute mean and store results
	vectordb mean_results(nominal_results*0);
	vector<vectordb> list_results;
	vectordb list_success(size_sample, 1.0);
	for (size_t i=0; i<size_sample; i++) {
		vectordb results = propagate_trajectory(
			ic_sample.getcol(i), x_goal,
			list_x, list_u,
			nav_sample.submat(
				0, i*N,
				list_x[0].nominal_state().size() - 1, (i+1)*N - 1),
			solver_parameters, spacecraft_parameters, dynamics,
			test_case_id,
			&mat_state_, &mat_control_,
			&mat_path_constraints_, &mat_terminal_constraints_,
			i<size_sample_saved);
		mean_results += results;
		list_results.push_back(results);

		// Check if the trajectory is valid
		if (results[1]>0)
			list_success[i] = 0.0;

		// Store trajectory
		if (i<size_sample_saved) {
			list_mat_state.push_back(mat_state_);
			list_mat_control.push_back(mat_control_);
			list_mat_path_constraints.push_back(mat_path_constraints_);
			list_mat_terminal_constraints.push_back(mat_terminal_constraints_);
		}
	}
	mean_results /= size_sample;

	// Get standard deviation
	vectordb std_results(nominal_results*0);
	for (size_t i=0; i<size_sample; i++) {
		std_results += sqr(list_results[i]-mean_results);
	}
	std_results = sqrt(std_results/size_sample);
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	// Get failure stats
	double failed_samples = size_sample - sqr(list_success.vnorm());
	double beta_r = max(0, failed_samples/(1.0*size_sample));
	double var_beta_r = beta_r*(1-beta_r);

	// Confidence intervals

	// Classic MC
	double probability_ci = 5e-2;
	double quantile_bound_traj = sqrt(inv_chi_2_cdf(1, 1 - solver_parameters.transcription_beta()));
	double quantile_bound_beta = sqrt(inv_chi_2_cdf(1, 1 - probability_ci));
	double beta_ci_size = quantile_bound_beta*sqrt(var_beta_r)/sqrt(1.0*size_sample);

	// From [Hanley Lippman-Hand 1983]
	double hanley_lippman_bound = -log(probability_ci);
	if (beta_r == 0 || beta_r == 1)
		beta_ci_size = hanley_lippman_bound/size_sample;
	vectordb ci_beta_r{max(beta_r - beta_ci_size, 0), min(beta_r + beta_ci_size, 1)};

	// Output
	if (verbosity < 3) {
		if (test_case_id != 0) {
			cout << "	Run time [s]: " << to_string(static_cast<double>(duration.count()) / 1e6) << endl;
			cout << "	Value - Final wet mass [kg] - Contraints violation [-]" << endl;
			cout << "	Nominal - "
				<< nominal_results[0] - spacecraft_parameters.dry_mass()*massu
				<< " - " << nominal_results[1] << endl;
			cout << "	MC - "
				<< mean_results[0] - spacecraft_parameters.dry_mass()*massu 
				<< " +/- " << std_results[0] 
				<< " - " << mean_results[1] 
				<< " +/- " << std_results[1] << endl;
			cout << "	Expected Beta quantile - "
				<< mean_results[0] - spacecraft_parameters.dry_mass()*massu - quantile_bound_traj*std_results[0] 
				<< " - " << nominal_results[1] + quantile_bound_traj*std_results[1] << endl;

			cout << "	Failed samples [-]: " << failed_samples << endl;
			cout << "	Failure rate [%] " << 100*beta_r 
				<<  ", within: [" << 100*ci_beta_r[0] << ", " << 100*ci_beta_r[1] << "] %" << endl;
			cout << "	Conservatism [-] " << conservatism(solver_parameters.transcription_beta(), beta_r) 
				<<  ", within: ["
				<< conservatism(solver_parameters.transcription_beta(), ci_beta_r[1])
				<< ", " << conservatism(solver_parameters.transcription_beta(), ci_beta_r[0]) << "]" << endl;
		} else {
			cout << "	Run time [s]: " << to_string(static_cast<double>(duration.count()) / 1e6) << endl;
			cout << "	Value - Final wet mass [kg] - Contraints violation [-]" << endl;
			cout << "	Nominal - "
				<< nominal_results[0] << endl;
			cout << "	MC - "
				<< mean_results[0] << " +/- " << std_results[0] << endl;
			cout << "	Expected Beta quantile - "
				<< mean_results[0] + quantile_bound_traj*std_results[0] << endl;

			cout << "	No constraints" << endl;
		}
	}
	else if (verbosity == 3) { // Benchmark mode
		// ID
		cout << test_case_id << ", ";
		cout << spacecraft_parameters.thrust()*thrustu/(spacecraft_parameters.initial_mass()*massu) << ", ";
		cout << N*list_x[0].nominal_state()[7]*tu*SEC2DAYS << ", ";
		cout << sqrt(list_x[0].Sigma().at(0, 0)) << ", ";
		cout << sqrt(list_x[0].Sigma().at(3, 3)) << ", ";
		cout << N << ", ";
		if (robust_optimisation)
			cout << solver_parameters.transcription_beta() << ", ";
		else
			cout << 1.0 << ", ";

		// Data

		// Results
		if (test_case_id != 0) {
			cout << nominal_results[0] - spacecraft_parameters.dry_mass()*massu << ", ";
			cout << nominal_results[1] << ", ";
		}
		else {
			cout << nominal_results[0]<< ", ";
			cout << "-, ";
		}
		cout << solver.nli() << ", ";

		// Convergence metrics
		cout << solver.AUL_runtime() << ", ";
		cout << solver.PN_runtime() << ", ";
		cout << solver.runtime() << ", ";
		cout << solver.DDP_n_iter() << ", ";
		cout << solver.AUL_n_iter() << ", ";
		cout << solver.PN_n_iter() << ", ";

		// MC properties
		cout << size_sample << ", ";
		cout << to_string(static_cast<double>(duration.count()) / 1e6) << ", ";

		// Mass
		if (test_case_id != 0) {
			cout << mean_results[0] - spacecraft_parameters.dry_mass()*massu << ", ";
			cout << std_results[0] << ", ";
			cout << mean_results[0] - spacecraft_parameters.dry_mass()*massu - quantile_bound_traj*std_results[0] << ", ";

			// Error
			cout << mean_results[1] << ", ";
			cout << std_results[1] << ", ";
			cout << mean_results[1] + quantile_bound_traj*std_results[1] << ", ";	

			// Beta r
			cout << beta_r << ", ";
			cout << beta_ci_size << ", ";
			cout << ci_beta_r[0] << ", ";
			cout << ci_beta_r[1] << ", ";	

			// Conservatism
			if (robust_optimisation) {
				cout << conservatism(solver_parameters.transcription_beta(), beta_r) << ", ";
				cout << conservatism(solver_parameters.transcription_beta(), ci_beta_r[1]) << ", ";
				cout << conservatism(solver_parameters.transcription_beta(), ci_beta_r[0]) << ", ";
			} else {
				cout << "-, -, -, ";
			}
		}
		else {
			cout << mean_results[0] << ", ";
			cout << std_results[0] << ", ";
			cout << mean_results[0] + quantile_bound_traj*std_results[0] << ", ";

			// Error
			cout << "-, -, -, ";

			// Beta r
			cout << "- , -, -, -, ";

			// Conservatism
			cout << "-, -, -, ";
		}
		cout <<  endl;
	}

	return vector<vector<matrixdb>>{
		list_mat_state, list_mat_control,
		list_mat_path_constraints, list_mat_terminal_constraints};
}
