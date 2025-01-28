/**
	double_integrator.cpp

	Purpose: Double_integrator execution script.

	@author Thomas Caleb

	@version 1.0 26/04/2023
*/

#include "test_cases.h"

using namespace DACE;
using namespace std::chrono;
using namespace std;

SolverParameters get_SolverParameters_double_integrator(
	unsigned int const& N, unsigned int const& DDP_type,
	matrixdb const& navigation_error_covariance,
	double const& transcription_beta,
	unsigned int verbosity) {
	// Solver parameters
	unsigned int Nx = 7;
	unsigned int Nu = 3;
	unsigned int Nineq = 0;
	unsigned int Ntineq = 0;
	bool with_J2 = false;
	double cost_to_go_gain = 1e-2;
	double terminal_cost_gain = 1e4;
	matrixdb terminal_cost_inv_covariance = inv(make_diag_matrix_(vectordb(1e-4, SIZE_VECTOR)));
	double mass_leak = 1e-8;
	double homotopy_coefficient = 0.0;
	double huber_loss_coefficient = 5e-3;
	vectordb homotopy_sequence{0};
	vectordb huber_loss_coefficient_sequence{1e-2};
	double DDP_tol = 1e-6;
	double AUL_tol = 1e-8; 
	double PN_tol = 1e-12;
	double PN_active_constraint_tol = 1e-13;
	unsigned int max_iter = 10000;
	unsigned int DDP_max_iter = 100;
	unsigned int AUL_max_iter = max_iter / DDP_max_iter;
	unsigned int PN_max_iter = 50;
	vectordb lambda_parameters{0.0, 1e8};
	vectordb mu_parameters{1, 1e8, 10};
	vectordb line_search_parameters{1e-8, 10.0, 0.5, 20};
	bool backward_sweep_regulation = true;
	vectordb backward_sweep_regulation_parameters{0, 1e-8, 1e8, 1.6};
	double PN_regularisation(1e-8);
	double PN_cv_rate_threshold(1.1);
	double PN_alpha(1.0); double PN_gamma(0.5);
	unsigned int saving_iterations = 0;

	return SolverParameters(
		N, Nx, Nu,
		Nineq, Ntineq, with_J2,
		cost_to_go_gain, terminal_cost_gain,
		terminal_cost_inv_covariance,
		navigation_error_covariance,
		transcription_beta,
		mass_leak,
		homotopy_coefficient, huber_loss_coefficient,
		homotopy_sequence,
		huber_loss_coefficient_sequence,
		DDP_type,
		DDP_tol, AUL_tol, PN_tol,
		DDP_max_iter, AUL_max_iter, PN_max_iter,
		line_search_parameters,
		backward_sweep_regulation,
		backward_sweep_regulation_parameters,
		lambda_parameters, mu_parameters,
		PN_regularisation, PN_active_constraint_tol,
		PN_cv_rate_threshold, PN_alpha, PN_gamma,
		verbosity, saving_iterations);
}

void double_integrator(int argc, char** argv) {
	// Input check
	if (argc < 11) {
		cout << "Wrong number of arguments." << endl;
		cout << "Requested number : 9" << endl;
		cout << "0 - Test case number." << endl;
		cout << "1 - SpacecraftParameter adress." << endl;
		cout << "2 - DDP type [0/1]." << endl;
		cout << "3 - Number of nodes [-]." << endl;
		cout << "4 - Time of flight [days]." << endl;
		cout << "5 - Perform robust optimisation [0/1]." << endl;
		cout << "6 - Target failure risk [0, 1]." << endl;
		cout << "7 - Perform fuel-optimal optimisation [0/1]." << endl;
		cout << "8 - Perform projected Newton solving [0/1]." << endl;
		cout << "9 - Save results [0/1]." << endl;
		cout << "10 - Verbosity [0-2]." << endl;
		return;
	}

	// Unpack inputs
	string spacecraft_parameters_file = argv[2];
	unsigned int DDP_type = atoi(argv[3]);
	unsigned int N = atoi(argv[4]);
	double ToF = atof(argv[5]);
	bool robust_solving = false;
	double transcription_beta = atof(argv[7]);
	bool fuel_optimal = false;
	bool pn_solving = false;
	bool save_results = false;
	int verbosity = atoi(argv[11]);
	if (atoi(argv[6]) == 1) { robust_solving = true; }
	if (atoi(argv[8]) == 1) { fuel_optimal = true; }
	if (atoi(argv[9]) == 1) { pn_solving = true; }
	if (atoi(argv[10]) == 1) { save_results = true; }

	// Set double precision
	typedef std::numeric_limits<double> dbl;
	cout.precision(5);

	// Set dynamics
	Dynamics dynamics = get_double_integrator_dynamics();
	
	// Normalisation constants
	Constants constants(dynamics.constants());

	// Spacecraft parameters
	SpacecraftParameters spacecraft_parameters(spacecraft_parameters_file);

	// Uncertainties
	double position_error = 1e-2; double velocity_error = 1e-2;
	vectordb init_convariance_diag{
		position_error, position_error, position_error,
		velocity_error, velocity_error, velocity_error,
		0.0};

	// Init solver parameters
	SolverParameters solver_parameters = get_SolverParameters_double_integrator(
		N, DDP_type, make_diag_matrix_(sqr(init_convariance_diag/100)),
		transcription_beta, verbosity);

	// Solver parameters
	unsigned int Nx = solver_parameters.Nx();
	unsigned int Nu = solver_parameters.Nu();
	N = solver_parameters.N();

	// Initial conditions [3*LU, 3*VU, MASSU, TU]
	vectordb x_departure{ 1.0, 1.0, 1.0, 1, 1, 1, 0};
	vectordb x_arrival{1.0, -1.0, 0.0, 0, 0, 0, 0};
	statedb x_goal(x_arrival); x_goal[Nx - 1] = ToF; // ToF

	statedb x0 = make_initial_state(
		x_departure, init_convariance_diag.extract(0, SIZE_VECTOR - 1), Nu);
	x0[Nx - 1] = 0; // Time step

	// First guess command
	controldb u_init = make_first_guess(vectordb(Nu, 2e-5), 1e-15, Nx);
	vector<controldb> list_u_init(N, u_init);

	// Set double precision
	typedef std::numeric_limits<double> dbl;
	cout.precision(5);

	// Make solver
	SODA solver(solver_parameters, spacecraft_parameters, dynamics);

	// Compute or load trajectory
	string file_name = "./data/robust_trajectory/double_integrator";
	string system_name = "DOUBLE INTEGRATOR";
	pair<vector<statedb>, vector<controldb>> robust_trajectory;
	bool load_trajectory=false;
	if (load_trajectory) {
		robust_trajectory = load_robust_trajectory(
			file_name, ToF, robust_solving,
			dynamics, spacecraft_parameters, constants, solver_parameters);
	} else {
		solver.solve(x0, list_u_init, x_goal, robust_solving, fuel_optimal, pn_solving);
		robust_trajectory = pair<vector<statedb>, vector<controldb>>(solver.list_x(), solver.list_u());
	}

	// Unpack
	vector<statedb> list_x(robust_trajectory.first);
	vector<controldb> list_u(robust_trajectory.second);
	
	// Print datasets
	if (save_results && !load_trajectory) {
		string file_name = "./data/robust_trajectory/double_integrator";
		string system_name = "DOUBLE INTEGRATOR";
		print_robust_trajectory_dataset(
			file_name, system_name,
			solver.list_trajectory_split(),
			x_departure, x_arrival, ToF, robust_solving,
			dynamics, spacecraft_parameters, constants, solver_parameters);
	}

	// Monte-Carlo validation
	bool monte_carlo_validaiton = true;
	if (monte_carlo_validaiton) {
		size_t size_sample=100000;
		vector<vector<matrixdb>> sample = test_trajectory(
			list_x, list_u, x_goal, size_sample,
			solver,	atoi(argv[1]), robust_solving, solver_parameters,
			solver.AULsolver().DDPsolver().spacecraft_parameters(),
			solver.AULsolver().DDPsolver().dynamics(),
			500);

		// Print datasets
		if (save_results) {
			string file_name = "./data/sample_trajectory/double_integrator";
			string system_name = "DOUBLE INTEGRATOR";
			print_sample_trajectory_dataset(
				file_name, system_name,
				sample[0], sample[1], sample[2], sample[3],
				ToF, robust_solving,
				spacecraft_parameters, constants, solver_parameters);
		}
	}
}
