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
	double const& transcription_beta,
	unsigned int verbosity) {
	// Solver parameters
	unsigned int Nx = (SIZE_VECTOR + 1) + 1;
	unsigned int Nu = SIZE_VECTOR / 2;
	unsigned int Nineq = 4;
	unsigned int Ntineq = 1;
	bool with_J2 = false;
	double cost_to_go_gain = 1e-1;
	double terminal_cost_gain = 1e8;
	matrixdb terminal_cost_inv_covariance = inv(make_diag_matrix_(
		vectordb{
			position_error_sqr, position_error_sqr, position_error_sqr,
			velocity_error_sqr, velocity_error_sqr, velocity_error_sqr}));
	double mass_leak = 2e-6;
	double homotopy_coefficient = 0.0;
	double huber_loss_coefficient = 5e-3;
	/* NO NAV 
	vectordb homotopy_sequence{0, 0.75, 0.9, 0.99};
	vectordb huber_loss_coefficient_sequence{1e-2, 1e-2, 2e-3, 1e-3};
	*/
	/* NAV */
	vectordb homotopy_sequence{0, 0.75, 0.9, 0.99};
	vectordb huber_loss_coefficient_sequence{1e-2, 1e-2, 2e-3, 2e-3};
	double DDP_tol = 1e-4;
	double AUL_tol = 1e-6;
	double PN_tol = 1e-12;
	double PN_active_constraint_tol = 1e-13;
	unsigned int max_iter = 10000;
	unsigned int DDP_max_iter = 100;
	unsigned int AUL_max_iter = 20;
	unsigned int PN_max_iter = 500;
	vectordb lambda_parameters{0.0, 1e8};
	vectordb mu_parameters{1, 1e8, 10};
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
		DDP_max_iter, AUL_max_iter, PN_max_iter,
		line_search_parameters,
		backward_sweep_regulation,
		backward_sweep_regulation_parameters,
		lambda_parameters, mu_parameters,
		PN_regularisation, PN_active_constraint_tol,
		PN_cv_rate_threshold, PN_alpha, PN_gamma,
		verbosity, saving_iterations);
}

void cr3bp_EARTH_MOON_lt_haloL2_to_haloL1(int argc, char** argv) {
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
	cout.precision(7);

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
	double position_error = 1e-8; double velocity_error = 5e-9;
	vectordb init_convariance_diag{
		position_error, position_error, position_error,
		velocity_error, velocity_error, velocity_error,
		0.0, 0.0};

	// Init solver parameters
	double terminal_position_error_sqr = sqr(1e-5); double terminal_velocity_error_sqr = sqr(1e-5);
	SolverParameters solver_parameters = get_SolverParameters_cr3bp_EARTH_MOON_lt_haloL2_to_haloL1(
		N, DDP_type, terminal_position_error_sqr, terminal_velocity_error_sqr,
		make_diag_matrix_(sqr(init_convariance_diag/100)), transcription_beta, verbosity);

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
		string file_name = "./data/robust_trajectory/cr3bp_EARTH_MOON_lt_haloL2_to_haloL1";
		string system_name = "CR3BP EARTH-MOON CARTESIAN LT";
		print_robust_trajectory_dataset(
			file_name, system_name,
			list_x, list_u, solver.list_nli(),
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
			string file_name = "./data/sample_trajectory/cr3bp_EARTH_MOON_lt_haloL2_to_haloL1";
			string system_name = "CR3BP EARTH-MOON CARTESIAN LT";
			print_sample_trajectory_dataset(
				file_name, system_name,
				sample[0], sample[1], sample[2], sample[3],
				ToF, robust_solving,
				spacecraft_parameters, constants, solver_parameters);
		}
	}
}
