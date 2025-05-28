/**
	test_cases.cpp

	Purpose: Executes a test case.

	@author Thomas Caleb

	@version 1.0 07/12/2023
*/

#include "test_cases.h"

using namespace DACE;
using namespace std::chrono;
using namespace std;

// Misc
void save_control(
	string const& file_name,
	vector<vectordb> const& list_u) {
	// Open file
	string file_name_(file_name);
	ofstream ofs(file_name_);

	// Set double precision
	typedef std::numeric_limits<double> dbl;
	ofs.precision(dbl::max_digits10);


	// Store the object to file
	size_t N = list_u.size();
	ofs << N << endl;

	for (size_t i=0; i<N; i++) {
		ofs << list_u[i];
	}

	ofs.close();
	return;
}
vector<vectordb> load_control(std::string const& file_name) {

	// Open file
	string file_name_(file_name);
	ifstream ifs(file_name_);

	// Load data
	string N_str; size_t N = 1;
	getline(ifs, N_str); 
	istringstream(N_str) >> N;

	vector<vectordb> list_u(N);
	for (size_t i=0; i<N; i++) {
		ifs >> list_u[i];
	}

	ifs.close();
	return list_u;
}

// Runs the selected test case.
void run_test_cases(int argc, char** argv) {

	// Input check
	if (argc < 2) {
		cout << "Wrong number of arguments." << endl;
		cout << "Requested number : 1" << endl;
		cout << "0 - Test case number from 0 to 4." << endl;
		return;
	}

	// Unpack inputs
	unsigned int test_case = atoi(argv[1]);

	// Excute correct test case

	// TBP
	
	// SUN centered

	// Earth-mars transfer
	if (test_case == 1)
		tbp_SUN_lt_earth_to_mars(argc, argv);

	// CR3BP

	// Earth-Moon
	
	// Halo L2 to Halo L1
	else if (test_case == 2)
		cr3bp_EARTH_MOON_lt_haloL2_to_haloL1(argc, argv);

	// Halo L2 (Gateway) to DRO
	else if (test_case == 3) {
		cr3bp_EARTH_MOON_lt_nrho_to_dro(argc, argv);
	}
	// Manyrev DRO to DRO
	else if (test_case == 4) {
		cr3bp_EARTH_MOON_lt_dro_to_dro(argc, argv);
	}
	// Lyapunov to Lyapunov
	else if (test_case == 5) {
		cr3bp_EARTH_MOON_lt_lyapunovL1_to_lyapunovL2(argc, argv);
	}
	
}
