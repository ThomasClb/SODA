/**
	parameters.cpp

	Purpose: Implementation of the SpacecraftParameters 
	and SolverParameters classes.

	@author Thomas Caleb

	@version 2.0 12/12/2024
*/

#include "parameters.h"

using namespace DACE;
using namespace std;


/*

	SPACECRAFT PARAMETERS

*/

// Constructors

// Default constructors
// Values from [Lantoine and Russell 2012]
// DOI : https://doi.org/10.1007/s10957-012-0039-0
SpacecraftParameters::SpacecraftParameters() {
	Constants constants;
	constants_ = constants;
	initial_mass_ = 1000 / constants.massu();
	dry_mass_ = 500 / constants.massu();
	thrust_ = 0.5 / constants.thrustu();
	Isp_ = 2000 / constants.tu();
}

// Constructors

// Values from [Lantoine and Russell 2012]
// DOI : https://doi.org/10.1007/s10957-012-0039-0
SpacecraftParameters::SpacecraftParameters(Constants const& constants) :
	constants_(constants),
	initial_mass_(1000 / constants.massu()), dry_mass_(500 / constants.massu()),
	thrust_(0.5 / constants.thrustu()), Isp_(2000 / constants.tu()) {}

SpacecraftParameters::SpacecraftParameters(
	Constants const& constants,
	double const& initial_mass,
	double const& dry_mass,
	double const& thrust,
	double const& Isp) :
	constants_(constants),
	initial_mass_(initial_mass), dry_mass_(dry_mass),
	thrust_(thrust), Isp_(Isp) {}

// Copy constructor
SpacecraftParameters::SpacecraftParameters(SpacecraftParameters const& param) :
	constants_(param.constants_),
	initial_mass_(param.initial_mass_), dry_mass_(param.dry_mass_),
	thrust_(param.thrust_), Isp_(param.Isp_) {}

// Loader constructor
SpacecraftParameters::SpacecraftParameters(string const& file_name) {
	// Open file
	string file_name_(file_name);
	ifstream ifs(file_name_);

	// Load data
	ifs >> *this;
	ifs.close();
	return;
}

// Destructor
SpacecraftParameters::~SpacecraftParameters() {}

// Getters
const Constants SpacecraftParameters::constants() const { return constants_; }
const double SpacecraftParameters::initial_mass() const {return initial_mass_; }
const double SpacecraftParameters::dry_mass() const { return dry_mass_; }
const double SpacecraftParameters::initial_wet_mass() const {
	return initial_mass_ - dry_mass_;
}
const double SpacecraftParameters::thrust() const { return thrust_; }
const double SpacecraftParameters::Isp() const { return Isp_; }
const double SpacecraftParameters::ejection_velocity() const {
	return (G_0 * constants_.tu() / (1000 * constants_.vu()))*Isp_;
}
const double SpacecraftParameters::mass_flow() const {
	return thrust_/ejection_velocity();
}

// IO operator
ostream& operator<<(ostream& os, const SpacecraftParameters& param) {
	// Set double precision
	typedef std::numeric_limits<double> dbl;
	os.precision(dbl::max_digits10);

	string legend(
		"SPACECRAFT PARAMETERS : INITIAL MASS [MASSU], DRY MASS [MASSU], THRUST [THRUSTU], ISP [TU].");

	// Write attributes
	os << param.constants();
	os << legend << endl;
	os << param.initial_mass() << endl;
	os << param.dry_mass() << endl;
	os << param.thrust() << endl;
	os << param.Isp() << endl;

	return os;
}
istream& operator>>(istream& is, SpacecraftParameters& param) {

	// Get constants
	is >> param.constants_;

	// Reading simple property from a line
	string initial_mass_str, dry_mass_str, thrust_str, Isp_str;
	getline(is, initial_mass_str); // Skip the legend
	getline(is, initial_mass_str);
	getline(is, dry_mass_str);
	getline(is, thrust_str);
	getline(is, Isp_str);
	istringstream(initial_mass_str) >> param.initial_mass_;
	istringstream(dry_mass_str) >> param.dry_mass_;
	istringstream(thrust_str) >> param.thrust_;
	istringstream(Isp_str) >> param.Isp_;

	return is;
}
void SpacecraftParameters::save(string const& file_name) const {
	// Open file
	string file_name_(file_name);
	ofstream ofs(file_name_);

	// Store the object to file
	ofs << *this;

	ofs.close();
}
void SpacecraftParameters::load(string const& file_name) {
	// Open file
	string file_name_(file_name);
	ifstream ifs(file_name_);

	// Load data
	ifs >> *this;

	ifs.close();
	return;
}

/*

	SOLVER PARAMETERS

*/

// Constructors

// Default constructor
SolverParameters::SolverParameters() :
	N_(40), Nx_((SIZE_VECTOR + 1) + 1), Nu_(SIZE_VECTOR / 2),
	Nineq_(2), Ntineq_(0),
	ToF_(0.0), with_J2_(false),
	stage_cost_gain_(1e-2), terminal_cost_gain_(1e4),
	terminal_cost_inv_covariance_(make_diag_matrix_(sqr(vectordb(SIZE_VECTOR + 1, 1e4)))),
	navigation_error_covariance_(make_diag_matrix_(sqr(vectordb(SIZE_VECTOR + 1, 1e-4)))),
	transcription_beta_(5e-2), mass_leak_(1e-8),
	homotopy_coefficient_(0), huber_loss_coefficient_(5e-3),
	homotopy_coefficient_sequence_(1, 0),
	huber_loss_coefficient_sequence_(1, 1e-2),
	DDP_type_(0),
	DDP_tol_(1e-4), AUL_tol_(1e-4), PN_tol_(1e-12),
	LOADS_tol_(1e-2), LOADS_max_depth_(1),
	DDP_max_iter_(100), AUL_max_iter_(100), PN_max_iter_(60),
	line_search_parameters_(vectordb{ 1e-8 , 10.0, 0.5, 20}),
	backward_sweep_regulation_(true),
	backward_sweep_regulation_parameters_(vectordb{ 0, 1e-8,  1e8, 1.6 }),
	lambda_parameters_(vectordb{ 0.0, 1e8 }),
	mu_parameters_(vectordb{ 1.0,  1e8, 10 }),
	PN_regularisation_(1e-8), PN_active_constraint_tol_(1e-13),
	PN_cv_rate_threshold_(1.1), PN_alpha_(1.0), PN_gamma_(0.5),
	list_lambda_(), list_mu_(), verbosity_(0), saving_iterations_(0) {
	// Unpack
	double lambda(lambda_parameters_[0]);
	double mu(mu_parameters_[0]);
	size_t size_path(Nineq_);
	size_t size_terminal(Ntineq_);

	// Init DACE
	DA::init(3, Nx_ + Nu_);
	DA::setEps(1e-90);

	// Quantiles
	path_quantile_ = sqrt(inv_chi_2_cdf(Nineq_ + 1, 1 - transcription_beta_));
	terminal_quantile_ = sqrt(inv_chi_2_cdf(Ntineq_ + 1, 1 - transcription_beta_));

	// Reserve space
	list_lambda_.reserve(N_ + 1);
	list_mu_.reserve(N_ + 1);
	
	// Path constraints dual state and penalties
	for (size_t i = 0; i < N_; i++) {
		list_lambda_.emplace_back(size_path, lambda);
		list_mu_.emplace_back(size_path, mu);
	}

	// Terminal constraints dual state and penalties
	list_lambda_.emplace_back(size_terminal, lambda);
	list_mu_.emplace_back(size_terminal, mu);
}

// Constructor
SolverParameters::SolverParameters(
	unsigned int const& N,
	unsigned int const& Nx, unsigned int const& Nu,
	unsigned int const& Nineq, unsigned int const& Ntineq,
	bool const& with_J2,
	double const& stage_cost_gain, double const& terminal_cost_gain,
	matrixdb const& terminal_cost_inv_covariance,
	matrixdb const& navigation_error_covariance,
	double const& transcription_beta, double const& mass_leak,
	double const& homotopy_coefficient, double const& huber_loss_coefficient,
	vectordb const& homotopy_coefficient_sequence,
	vectordb const& huber_loss_coefficient_sequence,
	unsigned int const& DDP_type,
	double const& DDP_tol, double const& AUL_tol, double const& PN_tol,
	double const& LOADS_tol, double const& LOADS_max_depth,
	unsigned int const& DDP_max_iter, unsigned int const& AUL_max_iter,
	unsigned int const& PN_max_iter,
	vectordb const& line_search_parameters, bool const& backward_sweep_regulation, 
	vectordb const& backward_sweep_regulation_parameters,
	vectordb const& lambda_parameters, vectordb const& mu_parameters,
	double const& PN_regularisation, double const& PN_active_constraint_tol,
	double const& PN_cv_rate_threshold, double const& PN_alpha,
	double const& PN_gamma, unsigned int const& verbosity,
	unsigned int const& saving_iterations) :
	N_(N), Nx_(Nx), Nu_(Nu),
	Nineq_(Nineq), Ntineq_(Ntineq), ToF_(0.0),
	with_J2_(with_J2),
	stage_cost_gain_(stage_cost_gain),
	terminal_cost_gain_(terminal_cost_gain),
	terminal_cost_inv_covariance_(terminal_cost_inv_covariance),
	navigation_error_covariance_(navigation_error_covariance),
	transcription_beta_(transcription_beta),
	mass_leak_(mass_leak),
	homotopy_coefficient_(homotopy_coefficient),
	huber_loss_coefficient_(huber_loss_coefficient),
	homotopy_coefficient_sequence_(homotopy_coefficient_sequence),
	huber_loss_coefficient_sequence_(huber_loss_coefficient_sequence),
	DDP_type_(DDP_type),
	DDP_tol_(DDP_tol), AUL_tol_(AUL_tol), PN_tol_(PN_tol),
	LOADS_tol_(LOADS_tol), LOADS_max_depth_(LOADS_max_depth),
	DDP_max_iter_(DDP_max_iter), AUL_max_iter_(AUL_max_iter), PN_max_iter_(PN_max_iter),
	line_search_parameters_(line_search_parameters), backward_sweep_regulation_(backward_sweep_regulation),
	backward_sweep_regulation_parameters_(backward_sweep_regulation_parameters),
	lambda_parameters_(lambda_parameters), mu_parameters_(mu_parameters),
	PN_regularisation_(PN_regularisation), PN_active_constraint_tol_(PN_active_constraint_tol),
	PN_cv_rate_threshold_(PN_cv_rate_threshold), PN_alpha_(PN_alpha),
	PN_gamma_(PN_gamma), list_lambda_(), list_mu_(),
	verbosity_(verbosity), saving_iterations_(saving_iterations) {
	vector<vectordb> list_lambda(N + 1, vectordb(Nineq, lambda_parameters[0]));
	vector<vectordb> list_mu(N + 1, vectordb(Nineq, mu_parameters[0]));
	// Unpack
	double lambda(lambda_parameters_[0]);
	double mu(mu_parameters_[0]);
	size_t size_path(Nineq_);
	size_t size_terminal( Ntineq_);

	// Init DACE
	DA::init(3, Nx_ + Nu_);
	DA::setEps(1e-90);

	// Quantiles
	path_quantile_ = sqrt(inv_chi_2_cdf(Nineq_ + 1, 1 - transcription_beta_));
	terminal_quantile_ = sqrt(inv_chi_2_cdf(Ntineq_ + 1, 1 - transcription_beta_));

	// Reserve space
	list_lambda_.reserve(N_ + 1);
	list_mu_.reserve(N_ + 1);

	// Path constraints dual state and penalties
	for (size_t i = 0; i < N_; i++) {
		list_lambda_.emplace_back(size_path, lambda);
		list_mu_.emplace_back(size_path, mu);
	}

	// Terminal constraints dual state and penalties
	list_lambda_.emplace_back(size_terminal, lambda);
	list_mu_.emplace_back(size_terminal, mu);
}

// Copy constructor
SolverParameters::SolverParameters(SolverParameters const& param) :
	N_(param.N_), Nx_(param.Nx_), Nu_(param.Nu_),
	Nineq_(param.Nineq_), Ntineq_(param.Ntineq_),
	ToF_(param.ToF_), with_J2_(param.with_J2_),
	stage_cost_gain_(param.stage_cost_gain_),
	terminal_cost_gain_(param.terminal_cost_gain_),
	terminal_cost_inv_covariance_(param.terminal_cost_inv_covariance_),
	navigation_error_covariance_(param.navigation_error_covariance_),
	transcription_beta_(param.transcription_beta_),
	path_quantile_(param.path_quantile_),
	terminal_quantile_(param.terminal_quantile_),
	mass_leak_(param.mass_leak_),
	homotopy_coefficient_(param.homotopy_coefficient_),
	huber_loss_coefficient_(param.huber_loss_coefficient_),
	homotopy_coefficient_sequence_(param.homotopy_coefficient_sequence_),
	huber_loss_coefficient_sequence_(param.huber_loss_coefficient_sequence_),
	DDP_type_(param.DDP_type_),
	DDP_tol_(param.DDP_tol_), AUL_tol_(param.AUL_tol_), PN_tol_(param.PN_tol_),
	LOADS_tol_(param.LOADS_tol_), LOADS_max_depth_(param.LOADS_max_depth_),
	DDP_max_iter_(param.DDP_max_iter_), AUL_max_iter_(param.AUL_max_iter_),
	PN_max_iter_(param.PN_max_iter_),
	list_lambda_(param.list_lambda_), list_mu_(param.list_mu_),
	line_search_parameters_(param.line_search_parameters_),
	backward_sweep_regulation_(param.backward_sweep_regulation_),
	backward_sweep_regulation_parameters_(param.backward_sweep_regulation_parameters_),
	lambda_parameters_(param.lambda_parameters_), mu_parameters_(param.mu_parameters_),
	PN_regularisation_(param.PN_regularisation_),
	PN_active_constraint_tol_(param.PN_active_constraint_tol_),
	PN_cv_rate_threshold_(param.PN_cv_rate_threshold_),
	PN_alpha_(param.PN_alpha_), PN_gamma_(param.PN_gamma_),
	verbosity_(param.verbosity_), saving_iterations_(param.saving_iterations_) {

	// Init DACE
	DA::init(3, Nx_ + Nu_);
	DA::setEps(1e-90);
}

// Destructor
SolverParameters::~SolverParameters() {}

// Getters
const unsigned int SolverParameters::N() const { return N_; }
const unsigned int SolverParameters::Nx() const { return Nx_; }
const unsigned int SolverParameters::Nu() const { return Nu_; }
const unsigned int SolverParameters::Nineq() const { return Nineq_; }
const unsigned int SolverParameters::Ntineq() const { return Ntineq_; }
const double SolverParameters::ToF() const { return ToF_; }
const bool SolverParameters::with_J2() const { return with_J2_; }
const double SolverParameters::stage_cost_gain() const { return stage_cost_gain_; }
const double SolverParameters::terminal_cost_gain() const { return terminal_cost_gain_; }
const matrixdb SolverParameters::terminal_cost_inv_covariance() const { return terminal_cost_inv_covariance_; }
const matrixdb SolverParameters::navigation_error_covariance() const { return navigation_error_covariance_; }
const double SolverParameters::transcription_beta() const { return transcription_beta_; }
const double SolverParameters::path_quantile() const { return path_quantile_; }
const double SolverParameters::terminal_quantile() const { return terminal_quantile_; }
const double SolverParameters::mass_leak() const { return mass_leak_; }
const double SolverParameters::homotopy_coefficient() const { return homotopy_coefficient_; }
const double SolverParameters::huber_loss_coefficient() const { return huber_loss_coefficient_; }
const vectordb SolverParameters::homotopy_coefficient_sequence() const {
	return homotopy_coefficient_sequence_; }
const vectordb SolverParameters::huber_loss_coefficient_sequence() const {
	return huber_loss_coefficient_sequence_; }
const unsigned int SolverParameters::DDP_type() const { return DDP_type_; }
const double SolverParameters::DDP_tol() const { return DDP_tol_; }
const double SolverParameters::AUL_tol() const { return AUL_tol_; }
const double SolverParameters::PN_tol() const { return PN_tol_; }
const double SolverParameters::LOADS_tol() const { return LOADS_tol_; }
const double SolverParameters::LOADS_max_depth() const { return LOADS_max_depth_; }
const unsigned int SolverParameters::DDP_max_iter() const { return DDP_max_iter_; }
const unsigned int SolverParameters::AUL_max_iter() const { return AUL_max_iter_; }
const unsigned int SolverParameters::PN_max_iter() const { return PN_max_iter_; }
const vector<vectordb> SolverParameters::list_lambda() const { return list_lambda_; }
const vector<vectordb> SolverParameters::list_mu() const { return list_mu_; }
const vectordb SolverParameters::line_search_parameters() const {
	return line_search_parameters_;
}
const bool SolverParameters::backward_sweep_regulation() const {
	return backward_sweep_regulation_;
}
const vectordb SolverParameters::backward_sweep_regulation_parameters() const {
	return backward_sweep_regulation_parameters_;
}
const vectordb SolverParameters::lambda_parameters() const { return lambda_parameters_; }
const vectordb SolverParameters::mu_parameters() const { return mu_parameters_; }
const double SolverParameters::PN_regularisation() const { return PN_regularisation_; }
const double SolverParameters::PN_active_constraint_tol() const {
	return PN_active_constraint_tol_;
}
const double SolverParameters::PN_cv_rate_threshold() const { return PN_cv_rate_threshold_; }
const double SolverParameters::PN_alpha() const { return PN_alpha_; }
const double SolverParameters::PN_gamma() const { return PN_gamma_; }
const unsigned int SolverParameters::verbosity() const { return verbosity_; }
const unsigned int SolverParameters::saving_iterations() const { return saving_iterations_; }

// Setters
void SolverParameters::set_navigation_error_covariance(matrixdb const& navigation_error_covariance) {
	navigation_error_covariance_ = navigation_error_covariance;}
void SolverParameters::set_transcription_beta(double const& transcription_beta) {
	transcription_beta_ = transcription_beta;}
void SolverParameters::set_path_quantile(double const& path_quantile) {
	path_quantile_ = path_quantile;}
void SolverParameters::set_terminal_quantile(double const& terminal_quantile) {
	terminal_quantile_ = terminal_quantile;}
void SolverParameters::set_homotopy_coefficient(double const& homotopy_coefficient) {
	homotopy_coefficient_ = homotopy_coefficient;}
void SolverParameters::set_huber_loss_coefficient(double const& huber_loss_coefficient) {
	huber_loss_coefficient_ = huber_loss_coefficient;
}
void SolverParameters::set_ToF(double const& ToF) {
	ToF_ = ToF;
}
void SolverParameters::set_list_lambda(vector<vectordb> const& list_lambda) {
	list_lambda_ = list_lambda; }
void SolverParameters::set_list_mu(vector<vectordb> const& list_mu) {
	list_mu_ = list_mu; }