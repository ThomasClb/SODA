/**
	cr3bp_EARTH_MOON_lt_haloL2_to_haloL1.cpp

	Purpose: Low-thrust Halo L2 to Halo L1 transfer execution script.
	In the Earth-Moon system.

	@author Thomas Caleb

	@version 1.0 20/11/2024
*/

#include "test_cases.h"

using namespace DACE;
using namespace std::chrono;
using namespace std;

SolverParameters get_SolverParameters_cr3bp_EARTH_MOON_lt_haloL2_to_haloL1(
	unsigned int const& N, unsigned int const& DDP_type,
	double const& position_error_sqr, double const& velocity_error_sqr, 
	matrixdb const& navigation_error_covariance,
	bool const& robust_solving,
	double const& transcription_beta, double const& LOADS_max_depth,
	unsigned int verbosity) {
	// Solver parameters
	unsigned int Nx = (SIZE_VECTOR + 1) + 1;
	unsigned int Nu = SIZE_VECTOR / 2;
	unsigned int Nineq = 2;
	unsigned int Ntineq = 1;
	bool with_J2 = false;
	double cost_to_go_gain = 1e-1;
	double terminal_cost_gain = 1e8;
	matrixdb terminal_cost_inv_covariance = make_diag_matrix_(
		vectordb{
			1/position_error_sqr, 1/position_error_sqr, 1/position_error_sqr,
			1/velocity_error_sqr, 1/velocity_error_sqr, 1/velocity_error_sqr});
	double mass_leak = 1e-5;
	double homotopy_coefficient = 0.0;
	double huber_loss_coefficient = 5e-3;
	vectordb homotopy_sequence{0, 0.5, 0.9, 0.995};
	vectordb huber_loss_coefficient_sequence{1e-2, 1e-2, 1e-2, 5e-4}; 
	double AUL_transcription_parameter = 2.5;
	vectordb PN_transcription_parameters{1.0, 1e-6, 1e-3, 0.5};
	if (
		(transcription_beta == 0.05 && LOADS_max_depth == 0.05 && robust_solving)
		) {
		vectordb homotopy_sequence{0, 0.5, 0.9, 0.995};
		vectordb huber_loss_coefficient_sequence{1e-2, 1e-2, 1e-2, 1e-3}; 
		AUL_transcription_parameter = 1;
	} else if (
		(transcription_beta == 0.05 && LOADS_max_depth == 0.5 && robust_solving)
		) {
		homotopy_sequence = vectordb{0, 0.5, 0.9, 0.99}; 
		huber_loss_coefficient_sequence = vectordb{1e-2, 1e-2, 1e-2, 1e-3}; 
		AUL_transcription_parameter = 2.5;
	}

	double DDP_tol = 1e-4;
	double AUL_tol = 1e-7;
	double PN_tol = 1e-12;
	double LOADS_tol = 1e-3;
	double PN_active_constraint_tol = 1e-13;
	unsigned int max_iter = 10000;
	unsigned int DDP_max_iter = 100;
	unsigned int AUL_max_iter = 100;
	unsigned int PN_max_iter = 5000;
	vectordb lambda_parameters{0.0, 1e8};
	vectordb mu_parameters{1, 1e8, 5};
	vectordb line_search_parameters{1e-10, 10.0, 0.5, 20};
	bool backward_sweep_regulation = true;
	vectordb backward_sweep_regulation_parameters{0, 1e-8, 1e15, 1.6};
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
		transcription_beta, mass_leak,
		homotopy_coefficient, huber_loss_coefficient,
		homotopy_sequence,
		huber_loss_coefficient_sequence,
		DDP_type,
		DDP_tol, AUL_tol, PN_tol,
		LOADS_tol, LOADS_max_depth,
		AUL_transcription_parameter,
		DDP_max_iter, AUL_max_iter, PN_max_iter,
		line_search_parameters,
		backward_sweep_regulation,
		backward_sweep_regulation_parameters,
		lambda_parameters, mu_parameters,
		PN_regularisation, PN_active_constraint_tol,
		PN_cv_rate_threshold, PN_alpha, PN_gamma,
		PN_transcription_parameters,
		verbosity, saving_iterations);
}

void cr3bp_EARTH_MOON_lt_haloL2_to_haloL1(int argc, char** argv) {
	// Input check
	if (argc < 14) {
		cout << "Wrong number of arguments." << endl;
		cout << "Requested number : 9" << endl;
		cout << "0 - Test case number." << endl;
		cout << "1 - SpacecraftParameter adress." << endl;
		cout << "2 - DDP type [0/1]." << endl;
		cout << "3 - Number of nodes [-]." << endl;
		cout << "4 - Time of flight [days]." << endl;
		cout << "5 - Perform robust optimisation [0/1]." << endl;
		cout << "6 - LOADS max depth [0, 1]." << endl;
		cout << "7 - Target failure risk [0, 1]." << endl;
		cout << "8 - Perform fuel-optimal optimisation [0/1]." << endl;
		cout << "9 - Perform projected Newton solving [0/1]." << endl;
		cout << "10 - Save results [0/1]." << endl;
		cout << "11 - Load trajectory [0/1]." << endl;
		cout << "12 - MC sample  size [-]." << endl;
		cout << "13 - Verbosity [0-2]." << endl;
		return;
	}

	// Unpack inputs
	string spacecraft_parameters_file = argv[2];
	unsigned int DDP_type = atoi(argv[3]);
	unsigned int N = atoi(argv[4]);
	double ToF = atof(argv[5]);
	bool robust_solving = false;
	double LOADS_max_depth = atof(argv[7]);
	double transcription_beta = atof(argv[8]);
	bool fuel_optimal = false;
	bool pn_solving = false;
	bool save_results = false;
	bool load_trajectory = false;
	unsigned int size_sample = atoi(argv[13]);
	int verbosity = atoi(argv[14]);
	if (atoi(argv[6]) == 1) { robust_solving = true; }
	if (atoi(argv[9]) == 1) { fuel_optimal = true; }
	if (atoi(argv[10]) == 1) { pn_solving = true; }
	if (atoi(argv[11]) == 1) { save_results = true; }
	if (atoi(argv[12]) == 1) { load_trajectory = true; }

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
	double position_error = 1e-7; double velocity_error = 1e-5;
	vectordb init_convariance_diag{
		position_error, position_error, position_error,
		velocity_error, velocity_error, velocity_error,
		0.0, 0.0};

	// Init solver parameters
	double terminal_position_error_sqr = sqr(position_error/10); double terminal_velocity_error_sqr = sqr(velocity_error/10);
	SolverParameters solver_parameters = get_SolverParameters_cr3bp_EARTH_MOON_lt_haloL2_to_haloL1(
		N, DDP_type, terminal_position_error_sqr, terminal_velocity_error_sqr,
		make_diag_matrix_(sqr(init_convariance_diag/100)), robust_solving,
		transcription_beta, LOADS_max_depth, verbosity);

	// Solver parameters
	unsigned int Nx = solver_parameters.Nx();
	unsigned int Nu = solver_parameters.Nu();
	N = solver_parameters.N();

	// Initial conditions [3*LU, 3*VU, MASSU, TU]
	ToF /=  (SEC2DAYS * tu); // [TU]
	double dt = ToF / N; // [TU]
	vectordb x_departure{
		1.1607973110000016, 0,
		-0.12269696820337475, 0,
		-0.20768326513738075, 0,
		spacecraft_parameters.initial_mass(),
		3.2746644337639852 };
	vectordb x_arrival{
		0.84871015300008812, 0,
		0.17388998538319206, 0,
		0.26350093896218163, 0,
		spacecraft_parameters.dry_mass(),
		2.5748200748171399 };
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
	string file_name = "./data/robust_trajectory/cr3bp_EARTH_MOON_lt_haloL2_to_haloL1";
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
		string file_name = "./data/sample_trajectory/cr3bp_EARTH_MOON_lt_haloL2_to_haloL1";
		print_sample_trajectory_dataset(
				file_name, system_name,
				sample[0], sample[1], sample[2], sample[3],
				ToF, robust_solving,
				spacecraft_parameters, constants, solver_parameters);
		}
	}
}
