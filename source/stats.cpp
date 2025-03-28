/**
	stats.cpp

	Purpose: Implementation of statistical functions.

	@author Thomas Caleb

	@version 1.0 02/10/2024
*/

#include "stats.h"

using namespace DACE;
using namespace std;

// Compares the constant parts' value.
bool compare_DA(DA const& a, DA const& b) {
	return cons(a) < cons(b);
}

// Sorting functions.
vectorDA sort_vector(vectorDA const& list_to_sort) {
	vectorDA list_copy = list_to_sort;
	stable_sort (list_copy.begin(), list_copy.end(), compare_DA);
	return list_copy;
}

// Sorting functions.
// DA version.
vectordb sort_vector(vectordb const& list_to_sort) {
	vectordb list_copy = list_to_sort;
	stable_sort (list_copy.begin(), list_copy.end());
	return list_copy;
}

// Inverse of the Cumulative Distribution Function of a chi squared distribution with d DoF.
double inv_chi_2_cdf(unsigned int const& d, double const& p) {
	// By dichotomy
	double r_min = EPS;
	double r_max = 1e5;
	double r = sqrt(r_min*r_max); // geometric mean

	// Safeguard
	if (d == 0)
		return 0.0;
	if (p <= 0)
		return 0.0;

	// Init
	double value_min = chi_2_cdf(d, r_min);
	double value_max = chi_2_cdf(d, r_max);
	double value = chi_2_cdf(d, r);

	// Dichotomy loop
	while (abs(value-p) > EPS) {
		if (value>p) {
			value_max = value; r_max = r;
		} else if (value<p) {
			value_min = value; r_min = r;
		} else
			return r;

		// Update
		r = sqrt(r_min*r_max); // geometric mean
		value = chi_2_cdf(d, r);
	}
	return r;
}

// Computes the conservatism.
// DOI: WIP
double conservatism(double const& beta, double const& beta_r) {
	return beta/beta_r*sqrt(((1-sqr(beta_r))/(1-sqr(beta))));
}
