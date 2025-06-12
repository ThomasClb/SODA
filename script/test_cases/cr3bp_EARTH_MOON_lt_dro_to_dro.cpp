/**
	cr3bp_EARTH_MOON_lt_dro_to_dro.cpp

	Purpose: Low-thrust DRO to DRO transfer execution script.
	In the Earth-Moon system.

	@author Thomas Caleb

	@version 1.0 20/11/2024
*/

#include "test_cases.h"

using namespace DACE;
using namespace std::chrono;
using namespace std;

SolverParameters get_SolverParameters_cr3bp_EARTH_MOON_lt_dro_to_dro(
	unsigned int const& N, 
	double const& position_error_sqr, double const& velocity_error_sqr, 
	matrixdb const& navigation_error_covariance,
	double const& transcription_beta, double const& LOADS_max_depth,
	bool const& robust_solving,
	unsigned int verbosity) {
	// Solver parameters
	unsigned int Nx = (SIZE_VECTOR + 1) + 1;
	unsigned int Nu = SIZE_VECTOR / 2;
	unsigned int Nineq = 2;
	unsigned int Ntineq = 1;
	double cost_to_go_gain = 1e-1;
	double terminal_cost_gain = 1e5;
	matrixdb terminal_cost_inv_covariance = make_diag_matrix_(
		vectordb{
			1/position_error_sqr, 1/position_error_sqr, 1/position_error_sqr,
			1/velocity_error_sqr, 1/velocity_error_sqr, 1/velocity_error_sqr});
	double mass_leak = 5e-6;
	double homotopy_coefficient = 0.0;
	double huber_loss_coefficient = 5e-3;
	vectordb homotopy_sequence, huber_loss_coefficient_sequence;
	double AUL_tol = 1e-6;
	double DDP_tol = 1e-4;
	double PN_tol = 1e-12;
	double LOADS_tol = 1e-3;
	double AUL_transcription_parameter = 1.0;
	double AUL_magnitude_perturbation = AUL_tol;
	double PN_active_constraint_tol = 1e-13;
	unsigned int max_iter = 10000;
	unsigned int DDP_max_iter = 100;
	unsigned int AUL_max_iter = 100;
	unsigned int PN_max_iter = 3000;
	vectordb mu_parameters{1, 1e8, 10};
	vectordb lambda_parameters{0.0, 1e8};
	vectordb line_search_parameters{1e-10, 10.0, 0.5, 20};
	bool backward_sweep_regulation = true;
	vectordb backward_sweep_regulation_parameters{0, 1e-8, 1e15, 1.5};
	double PN_regularisation(1e-8);
	double PN_cv_rate_threshold(1.1);
	double PN_alpha(1.0); double PN_gamma(0.5);
	vectordb PN_transcription_parameters{1.0, 1e-6, 1e-3, 0.2};

	if (!robust_solving){
		homotopy_sequence = vectordb{0, 0.5, 0.9, 0.995}; 
		huber_loss_coefficient_sequence = vectordb{1e-2, 2e-3, 1e-3, 5e-4};
	} else if (
		(transcription_beta == 0.05 && LOADS_max_depth == 0.5) ) { 
		homotopy_sequence = vectordb{0, 0.5, 0.9, 0.995}; 
		huber_loss_coefficient_sequence = vectordb{1e-2, 2e-3, 1e-3, 5e-4};
		AUL_transcription_parameter = 2.0;
	} else if (
		(transcription_beta == 0.05 && LOADS_max_depth == 0.05)) {
		AUL_magnitude_perturbation = 5*AUL_tol; // 5
		homotopy_sequence = vectordb{0, 0.85, 0.995};
		huber_loss_coefficient_sequence = vectordb{1e-2, 2e-3, 5e-4}; 
		mu_parameters[2] = 5;
	}

	return SolverParameters(
		N, Nx, Nu,
		Nineq, Ntineq,
		cost_to_go_gain, terminal_cost_gain,
		terminal_cost_inv_covariance,
		navigation_error_covariance,
		transcription_beta, mass_leak,
		homotopy_coefficient, huber_loss_coefficient,
		homotopy_sequence,
		huber_loss_coefficient_sequence,
		DDP_tol, AUL_tol, PN_tol,
		LOADS_tol, LOADS_max_depth,
		AUL_transcription_parameter,
		AUL_magnitude_perturbation,
		DDP_max_iter, AUL_max_iter, PN_max_iter,
		line_search_parameters,
		backward_sweep_regulation,
		backward_sweep_regulation_parameters,
		lambda_parameters, mu_parameters,
		PN_regularisation, PN_active_constraint_tol,
		PN_cv_rate_threshold, PN_alpha, PN_gamma,
		PN_transcription_parameters,
		verbosity);
}

void cr3bp_EARTH_MOON_lt_dro_to_dro(int argc, char** argv) {
	// Input check
	if (argc < 13) {
		cout << "Wrong number of arguments." << endl;
		cout << "Requested number : 12" << endl;
		cout << "0 - SpacecraftParameter adress." << endl;
		cout << "1 - DDP type [0/1]." << endl;
		cout << "2 - Number of nodes [-]." << endl;
		cout << "3 - Time of flight [days]." << endl;
		cout << "4 - Perform robust optimisation [0/1]." << endl;
		cout << "5 - LOADS max depth [0, 1]." << endl;
		cout << "6 - Target failure risk [0, 1]." << endl;
		cout << "7 - Perform fuel-optimal optimisation [0/1]." << endl;
		cout << "8 - Perform projected Newton solving [0/1]." << endl;
		cout << "9 - Save results [0/1]." << endl;
		cout << "10 - Load trajectory [0/1]." << endl;
		cout << "11 - MC sample  size [-]." << endl;
		cout << "12 - Verbosity [0-2]." << endl;
		return;
	}

	// Unpack inputs
	string spacecraft_parameters_file = argv[2];
	unsigned int N = atoi(argv[3]);
	double ToF = atof(argv[4]);
	bool robust_solving = false;
	double LOADS_max_depth = atof(argv[6]);
	double transcription_beta = atof(argv[7]);
	bool fuel_optimal = false;
	bool pn_solving = false;
	bool save_results = false;
	bool load_trajectory = false;
	if (atoi(argv[5]) == 1) { robust_solving = true; }
	if (atoi(argv[8]) == 1) { fuel_optimal = true; }
	if (atoi(argv[9]) == 1) { pn_solving = true; }
	if (atoi(argv[10]) == 1) { save_results = true; }
	if (atoi(argv[11]) == 1) { load_trajectory = true; }
	unsigned int size_sample = atoi(argv[12]);
	int verbosity = atoi(argv[13]);

	// Get dynamics
	Dynamics dynamics = get_cr3bp_EARTH_MOON_lt_dynamics();

	// Normalisation cosntants
	Constants constants(dynamics.constants());
	double lu = constants.lu();
	double massu = constants.massu();
	double tu = constants.tu();
	double thrustu = constants.thrustu();
	double vu = constants.vu();

	// Spacecraft parameters
	SpacecraftParameters spacecraft_parameters(spacecraft_parameters_file);

	// Uncertainties
	double position_error = 5e-6; double velocity_error = 5e-5; 
	vectordb init_convariance_diag{
		position_error, position_error, position_error/100,
		velocity_error, velocity_error, velocity_error/100,
		0.0, 0.0};

	// Init solver parameters
	double terminal_position_error_sqr = sqr(position_error/10);
	double terminal_velocity_error_sqr = sqr(velocity_error/10);
	SolverParameters solver_parameters = get_SolverParameters_cr3bp_EARTH_MOON_lt_dro_to_dro(
		N, terminal_position_error_sqr, terminal_velocity_error_sqr,
		make_diag_matrix_(sqr(init_convariance_diag/100)),
		transcription_beta, LOADS_max_depth, robust_solving, verbosity);

	// Solver parameters
	unsigned int Nx = solver_parameters.Nx();
	unsigned int Nu = solver_parameters.Nu();
	N = solver_parameters.N();

	// Initial conditions [3*LU, 3*VU, MASSU, TU]
	ToF /=  (SEC2DAYS * tu); // [TU]
	double dt = ToF / N; // [TU]
	vectordb x_departure{ 
		1.171359, 0, 0.0,
		0, -0.489458, 0.0,
		spacecraft_parameters.initial_mass(),
		13.4 / SEC2DAYS / tu };
	vectordb x_arrival{
		1.301844, 0, 0.0,
		0, -0.642177, 0.0,
		spacecraft_parameters.dry_mass(),
		21.6 / SEC2DAYS / tu };
	statedb x_goal = x_arrival; x_goal[Nx - 1] = ToF; // ToF
	statedb x0 = make_initial_state(
		x_departure, init_convariance_diag.extract(0, SIZE_VECTOR - 1), Nu);
	x0[Nx - 1] = dt; // Time step

	// First guess command
	controldb u_init = make_first_guess(vectordb(Nu, 1e-6 / thrustu), 1e-15 / thrustu, Nx);
	vector<controldb> list_u_init(N, u_init);

	// Set double precision
	typedef std::numeric_limits<double> dbl;
	cout.precision(5);

	// Make solver
	SODA solver(solver_parameters, spacecraft_parameters, dynamics);

	// Compute or load trajectory
	string file_name = "./data/robust_trajectory/cr3bp_EARTH_MOON_lt_dro_to_dro";
	string system_name = "CR3BP EARTH-MOON CARTESIAN LT";
	RobustTrajectory robust_trajectory;
	if (load_trajectory) {
		robust_trajectory = load_robust_trajectory(
			file_name, ToF, robust_solving,
			dynamics, spacecraft_parameters, constants, solver_parameters);
		
	} else {
		solver.solve(x0, list_u_init, x_goal, robust_solving, fuel_optimal, pn_solving);
		robust_trajectory = RobustTrajectory(solver.list_trajectory_split());
	}
	
	// Print datasets
	if (save_results && !load_trajectory) {
		print_robust_trajectory_dataset(
			file_name, system_name,
			robust_trajectory,
			x_departure, x_arrival, ToF, robust_solving,
			dynamics, spacecraft_parameters, constants, solver_parameters);
	}	

	// Monte-Carlo validation
	if (size_sample > 0) {
		vector<vector<matrixdb>> sample = test_trajectory(
			robust_trajectory,
			x0, x_goal, size_sample,
			solver,	atoi(argv[1]), robust_solving, solver_parameters,
			solver.AULsolver().DDPsolver().spacecraft_parameters(),
			solver.AULsolver().DDPsolver().dynamics(),
			500);

		// Print datasets
		if (save_results) {
			string file_name = "./data/sample_trajectory/cr3bp_EARTH_MOON_lt_dro_to_dro";
			print_sample_trajectory_dataset(
				file_name, system_name,
				sample[0], sample[1], sample[2], sample[3],
				ToF, robust_solving,
				spacecraft_parameters, constants, solver_parameters);
		}
	}
}
