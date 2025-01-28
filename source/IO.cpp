/**
	IO.cpp

	Purpose: Implementation of the inputs and outputs of data.

	@author Thomas Caleb

	@version 2.0 13/12/2024
*/

#include "IO.h"

using namespace DACE;
using namespace std;

// Splits a string into substring.
vector<string> split(string s, string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    string token;
    vector<string> res;
    while ((pos_end = s.find(delimiter, pos_start)) != string::npos) {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back (token);
    }
    res.push_back (s.substr (pos_start));
    return res;
}

// Function to print a dataset at a given name in order to
// produce python visuals.
void print_dataset(
	string const& file_name,
	string const& system_name,
	SpacecraftParameters const& spacecraft_parameters,
	vector<vector<string>> const& list_title,
	vector<matrixdb> const& list_data) {

	// Open file
	ofstream ofs(file_name);

	// Print header
	ofs << file_name << endl;
	ofs << system_name << endl;
	ofs << spacecraft_parameters;

	// Print lists
	
	size_t nb_lists = list_data.size();
	ofs << nb_lists << endl; // Number of lists
	
	for (size_t i = 0; i < nb_lists; i++) {
		// Unpack
		vector<string> list_title_i = list_title[i];
		matrixdb data_i = list_data[i];
		size_t nb_rows = data_i.nrows();
		size_t nb_colomns = data_i.ncols();

		// Print

		// Headers
		ofs << list_title_i[0] << endl; // Series title
		ofs << nb_rows << ", " << nb_colomns << endl; // Dataset size
		
		for (size_t l = 0; l < nb_colomns; l++) {
			ofs << list_title_i[1 + l] << ", ";
		}
		ofs << endl;
		
		// Data
		for (size_t k = 0; k < nb_rows; k++) {
			for (size_t l = 0; l < nb_colomns; l++) {
				ofs << data_i.at(k,l) << ", ";
			}
			ofs << endl;
		}
	}
	ofs.close();
	return;
}

// Function to propagate a vector without control.
matrixdb get_mat_reference_trajectory(
	vectordb const& x_0,
	Dynamics const& dynamics,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters,
	int const& nb_point) {
	// Unpack
	int Nu = solver_parameters.Nu();
	int Nx = solver_parameters.Nx();
	double period = x_0[Nx - 1];

	// Init
	matrixdb output(nb_point, Nx);
	vectordb x_i = x_0;
	double dt = period / (nb_point - 1.0);
	x_i[Nx - 1] = dt;
	vectordb null_control(Nu, 0.0);
	
	// Loop
	for (size_t i = 0; i < nb_point; i++) {
		output.setrow(i, x_i);
		x_i = dynamics.dynamic_db()(
			x_i, null_control,
			spacecraft_parameters, constants, solver_parameters);
	}
	output.setrow(nb_point-1, x_i);
	return output;
}

// Prints a transfer dataset on a standardised format
void print_sample_trajectory_dataset(
	string const& file_name, string const& system_name,
	vector<matrixdb> const& list_mat_state,
	vector<matrixdb> const& list_mat_control,
	vector<matrixdb> const& list_mat_path_constraints,
	vector<matrixdb> const& list_mat_terminal_constraints,
	double const& ToF, bool const& robust_solving,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters) {
	// Unpack
	unsigned int N = solver_parameters.N();
	unsigned int Nx = solver_parameters.Nx();
	unsigned int Nu = solver_parameters.Nu();
	unsigned int Nineq = solver_parameters.Nineq();
	unsigned int Ntineq = solver_parameters.Ntineq();
	size_t size_sample = list_mat_state.size();

	// Init
	vector<matrixdb> list_data;
	list_data.reserve(size_sample*4);

	// Make list title
	vector<vector<string>> list_title;
	list_title.reserve(size_sample*4);
	vector<string> title_state{
		" ",
		"x [LU]", "y [LU]", "z [LU]",
		"vx [VU]", "vy [VU]", "vz [VU]",
		"mass [MASSU]", "dt [TU]" };
	vector<string> title_control{
		" ",
		"ux [THRUSTU]", "uy [THRUSTU]", "uz [THRUSTU]"};
	vector<string> title_path_constraints(Nineq + 1, "TBD [TBD]");
	vector<string> title_terminal_constraints(Ntineq + 1, "TBD [TBD]");

	// Loop on all sample
	for (size_t i=0; i<size_sample; i++) {
		// Title
		title_state[0] = "Sample state " + to_string(i);
		title_control[0] = "Sample control " + to_string(i);
		title_path_constraints[0] = "Sample path constraints " + to_string(i);
		title_terminal_constraints[0] = "Sample terminal constraints " + to_string(i);
		list_title.push_back(title_state);
		list_title.push_back(title_control);
		list_title.push_back(title_path_constraints);
		list_title.push_back(title_terminal_constraints);

		// Data
		list_data.push_back(list_mat_state[i]);
		list_data.push_back(list_mat_control[i]);
		list_data.push_back(list_mat_path_constraints[i]);
		list_data.push_back(list_mat_terminal_constraints[i]);
	}

	// Make file name
	int power = 9;
	double thrust_mass = ((spacecraft_parameters.thrust() * spacecraft_parameters.constants().thrustu()
		/ (spacecraft_parameters.initial_mass() * spacecraft_parameters.constants().massu())
		) * pow(10.0, power));
	int inv_beta = 1/solver_parameters.transcription_beta();
	if (!robust_solving)
		inv_beta = 0;
	int exposant = log10(thrust_mass);
	int mantisse = static_cast<int>(thrust_mass) / static_cast<int>(pow(10, exposant));
	string str_T2m = to_string(mantisse) + "e" + to_string(exposant - power);
	string file_name_ = file_name + "_"
		+ to_string(inv_beta) + "_"
		+ str_T2m + "_"
		+ to_string((int)(ToF*spacecraft_parameters.constants().tu()*SEC2DAYS)) + "_"
		+ to_string(solver_parameters.DDP_type())
		+ ".dat";

	print_dataset(
		file_name_, system_name,
		spacecraft_parameters,
		list_title, list_data);
}

// Prints a transfer dataset on a standardised format
// with reference orbits.
void print_robust_trajectory_dataset(
	string const& file_name, string const& system_name,
	deque<TrajectorySplit> const& robust_trajectory,
	vectordb const& x_0, vectordb const& x_f,
	double const& ToF, bool const& robust_solving,
	Dynamics const& dynamics,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters) {
	// Unpack
	int Nu = solver_parameters.Nu();
	int Nx = solver_parameters.Nx();
	int N = solver_parameters.N();

	// Get reference trajectories
	matrixdb mat_departure = get_mat_reference_trajectory(
		x_0, dynamics,
		spacecraft_parameters, constants, solver_parameters,
		2*N);
	matrixdb mat_arrival = get_mat_reference_trajectory(
		x_f, dynamics,
		spacecraft_parameters, constants, solver_parameters,
		2*N);

	// Make lists
	vector<string> title_history{"Splitting history", "Splitting direction", "Splitting side"};
	vector<string> title_state{
		"Nominal state",
		"x [LU]", "y [LU]", "z [LU]",
		"vx [VU]", "vy [VU]", "vz [VU]",
		"mass [MASSU]", "dt [TU]" };
	vector<string> title_control{
		"Nominal control",
		"ux [THRUSTU]", "uy [THRUSTU]", "uz [THRUSTU]"};
	vector<string> title_der_dynamics{"Der dynamics"};
	vector<string> title_feedback_gain{"Feedback gain"};
	vector<string> title_Sigma{"Sigma"};
	for (size_t i=0; i<Nx; i++) {
		string A(title_state[1+i]), B;
		for (size_t j=0; j<Nx; j++) {
			B = title_state[1+j];
			title_Sigma.push_back(A + " * " + B);
			title_der_dynamics.push_back(A + " * " + B);
		}
		for (size_t j=0; j<Nu; j++) {
			B = title_control[1+j];
			title_der_dynamics.push_back(A + " * " + B);
		}
	}
	for (size_t i=0; i<Nu; i++) {
		string A(title_control[1+i]), B;
		for (size_t j=0; j<Nx; j++) {
			B = title_state[1+j];
			title_feedback_gain.push_back(A + " * " + B);
		}
	}
	vector<string> title_departure(title_state), title_arrival(title_state);
	title_departure[0] = "Departure orbit";
	title_arrival[0] = "Arrival orbit";
	vector<vector<string>> list_title{title_departure, title_arrival};
	list_title.reserve(robust_trajectory.size()*6);
	vector<matrixdb> list_data{mat_departure, mat_arrival};
	list_data.reserve(robust_trajectory.size()*6);

	// Loop on all splits
	matrixdb mat_nominal_state(N+1, Nx), mat_der_dynamics(N+1, Nx*(Nx+Nu)), mat_Sigma(N+1, Nx*Nx);
	matrixdb mat_nominal_control(N, Nu), mat_feeedback_gain(N, Nx*Nu);
	for (size_t k=0; k<robust_trajectory.size(); k++) {
		TrajectorySplit trajectory_split(robust_trajectory[k]);
		SplittingHistory history(trajectory_split.splitting_history());
		matrixdb mat_history(history.size(), 2);
		for (size_t i=0; i<history.size(); i++) {
			mat_history.at(i, 0) = history[i].first;
			mat_history.at(i, 1) = history[i].second;
		}

		// Make matrices
		for (size_t i=0; i<N; i++) {
			statedb x_i(trajectory_split.list_x()[i]);
			controldb u_i(trajectory_split.list_u()[i]);
			mat_nominal_state.setrow(i, x_i.nominal_state());
			mat_der_dynamics.setrow(i, matrix_to_vector(x_i.der_dynamics())); 
			mat_Sigma.setrow(i, matrix_to_vector(x_i.Sigma()));
			mat_nominal_control.setrow(i, u_i.nominal_control());
			mat_feeedback_gain.setrow(i, matrix_to_vector(u_i.feedback_gain()));
		}
		statedb x_i(trajectory_split.list_x()[N]);
		mat_nominal_state.setrow(N, x_i.nominal_state());
		mat_der_dynamics.setrow(N, matrix_to_vector(x_i.der_dynamics())); 
		mat_Sigma.setrow(N, matrix_to_vector(x_i.Sigma()));

		// Add to lists

		// Titles
		vector<string> title_history_k(title_history);
		vector<string> title_state_k(title_state);
		vector<string> title_control_k(title_control);
		vector<string> title_der_dynamics_k(title_der_dynamics);
		vector<string> title_Sigma_k(title_Sigma);
		vector<string> title_feedback_gain_k(title_feedback_gain);
		title_history_k[0] = title_history_k[0] + " GMM " + to_string(k);
		title_state_k[0] = title_state_k[0] + " GMM " + to_string(k);
		title_control_k[0] = title_control_k[0] + " GMM " + to_string(k);
		title_der_dynamics_k[0] = title_der_dynamics_k[0] + " GMM " + to_string(k);
		title_Sigma_k[0] = title_Sigma_k[0] + " GMM " + to_string(k);
		title_feedback_gain_k[0] = title_feedback_gain_k[0] + " GMM " + to_string(k);
		list_title.push_back(title_history_k);
		list_title.push_back(title_state_k);
		list_title.push_back(title_der_dynamics_k);
		list_title.push_back(title_Sigma_k);
		list_title.push_back(title_control_k);
		list_title.push_back(title_feedback_gain_k);

		// Data
		list_data.push_back(mat_history);
		list_data.push_back(mat_nominal_state);
		list_data.push_back(mat_der_dynamics);
		list_data.push_back(mat_Sigma);
		list_data.push_back(mat_nominal_control);
		list_data.push_back(mat_feeedback_gain);
	}

	// Make file name
	int power = 9;
	double thrust_mass = ((spacecraft_parameters.thrust() * spacecraft_parameters.constants().thrustu()
		/ (spacecraft_parameters.initial_mass() * spacecraft_parameters.constants().massu())
		) * pow(10.0, power));
	int inv_beta = 1/solver_parameters.transcription_beta();
	if (!robust_solving)
		inv_beta = 0;
	int exposant = log10(thrust_mass);
	int mantisse = static_cast<int>(thrust_mass) / static_cast<int>(pow(10, exposant));
	string str_T2m = to_string(mantisse) + "e" + to_string(exposant - power);
	string file_name_ = file_name + "_"
		+ to_string(inv_beta) + "_"
		+ str_T2m + "_"
		+ to_string((int)(ToF*spacecraft_parameters.constants().tu()*SEC2DAYS)) + "_"
		+ to_string(solver_parameters.DDP_type())
		+ ".dat";

	print_dataset(
		file_name_, system_name,
		spacecraft_parameters,
		list_title, list_data);
}

// Loads a printed robust trajectory from a printed file.
pair<vector<statedb>, vector<controldb>> load_robust_trajectory(
	string const& file_name,
	double const& ToF, bool const& robust_solving,
	Dynamics const& dynamics,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters) {

	// Make file name
	int power = 9;
	double thrust_mass = ((spacecraft_parameters.thrust() * spacecraft_parameters.constants().thrustu()
		/ (spacecraft_parameters.initial_mass() * spacecraft_parameters.constants().massu())
		) * pow(10.0, power));
	int inv_beta = 1/solver_parameters.transcription_beta();
	if (!robust_solving)
		inv_beta = 0;
	int exposant = log10(thrust_mass);
	int mantisse = static_cast<int>(thrust_mass) / static_cast<int>(pow(10, exposant));
	string str_T2m = to_string(mantisse) + "e" + to_string(exposant - power);
	string file_name_ = file_name + "_"
		+ to_string(inv_beta) + "_"
		+ str_T2m + "_"
		+ to_string((int)(ToF*spacecraft_parameters.constants().tu()*SEC2DAYS)) + "_"
		+ to_string(solver_parameters.DDP_type())
		+ ".dat";

	// Open file
	ifstream ifs(file_name_);
	string delimiter(", ");

	// Skip preamble
	string buffer_str;
	for (size_t i=0; i<14; i++) {getline(ifs, buffer_str);}

	// Get N, Nx
	getline(ifs, buffer_str);
	vector<string> list_str = split(buffer_str, delimiter);
	unsigned int N(stoi(list_str[0])-1), Nx(stoi(list_str[1]));

	// Init
	vector<statedb> list_x(N+1);
	vector<controldb> list_u(N);

	// Skip header
	getline(ifs, buffer_str);

	// Get nominal states
	statedb x;
	vectordb nominal_state(Nx);
	for (size_t i=0; i<N+1; i++) {
		// Get string
		getline(ifs, buffer_str);
		list_str = split(buffer_str, delimiter);

		// Assign
		for (size_t j=0; j<Nx; j++) {
			nominal_state[j] = stod(list_str[j]);
		}
		x.set_nominal_state(nominal_state);
		list_x[i] = x;
	}

	// Skip header
	getline(ifs, buffer_str);

	// Get Nu
	getline(ifs, buffer_str);
	list_str = split(buffer_str, delimiter);
	unsigned int Nu(stoi(list_str[1])/Nx - Nx);

	// Skip header
	getline(ifs, buffer_str);

	// Get der dynamics
	matrixdb der_dynamics(Nx, Nx + Nu), Sigma(Nx, Nx);
	for (size_t i=0; i<N+1; i++) {
		// Get string
		getline(ifs, buffer_str);
		list_str = split(buffer_str, delimiter);

		// Assign
		for (size_t j=0; j<Nx*(Nx+Nu); j++) {
			der_dynamics.at(j/(Nx+Nu),j%(Nx+Nu)) = stod(list_str[j]);
		}
		list_x[i].set_der_dynamics(der_dynamics);
	}

	// Skip 3*header
	getline(ifs, buffer_str); getline(ifs, buffer_str); getline(ifs, buffer_str);

	// Get Sigma
	for (size_t i=0; i<N+1; i++) {
		// Get string
		getline(ifs, buffer_str);
		list_str = split(buffer_str, delimiter);

		// Assign
		for (size_t j=0; j<Nx*Nx; j++) {
			Sigma.at(j/Nx,j%Nx) = stod(list_str[j]);
		}
		list_x[i].set_Sigma(Sigma);
	}

	// Skip 3*header
	getline(ifs, buffer_str);
	getline(ifs, buffer_str);
	getline(ifs, buffer_str);

	// Get nominal controls
	controldb u;
	vectordb nominal_control(Nu);
	matrixdb feedback_gain(Nu, Nx);
	for (size_t i=0; i<N; i++) {
		// Get string
		getline(ifs, buffer_str);
		list_str = split(buffer_str, delimiter);

		// Assign
		for (size_t j=0; j<Nu; j++) {
			nominal_control[j] = stod(list_str[j]);
		}
		u.set_nominal_control(nominal_control);
		list_u[i] = u;
	}

	// Skip 3*header
	getline(ifs, buffer_str); getline(ifs, buffer_str); getline(ifs, buffer_str);

	// Get feedback gains
	for (size_t i=0; i<N; i++) {
		// Get string
		getline(ifs, buffer_str);
		list_str = split(buffer_str, delimiter);

		// Assign
		for (size_t j=0; j<Nx*Nu; j++) {
			feedback_gain.at(j/Nx,j%Nx) = stod(list_str[j]);
		}
		list_u[i].set_feedback_gain(feedback_gain);
	}

	// Close
	ifs.close();
	return pair<vector<statedb>, vector<controldb>>(list_x, list_u);
}