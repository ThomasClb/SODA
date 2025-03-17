/**
	stats.h

	Purpose: Implementation of statistical functions.

	@author Thomas Caleb

	@version 1.0 02/10/2024
*/

#ifndef DEF_STATS
#define DEF_STATS

#pragma once

#include <vector>
#include <cmath>

#include <dace/dace_s.h>

#include "linalg.h"
#include "derivation.h"

// Compares the constant parts' value.
bool compare_DA(DACE::DA const& a, DACE::DA const& b);

// Sorting functions.
DACE::vectordb sort_vector(DACE::vectordb const& list_to_sort);
DACE::vectorDA sort_vector(DACE::vectorDA const& list_to_sort); // DA version.

// Regularised incomplete Euler beta-function I_x(a,b).
// See: https://fr.wikipedia.org/wiki/Fonction_b%C3%AAta#Fonction_b%C3%AAta_incompl%C3%A8te
// We only need to compute values of type I_x(a=(d-1)/2, 0.5).
// Method from [DiDonato and Jarnagin 1949].
// See: https://apps.dtic.mil/sti/tr/pdf/AD0642495.pdf
template<typename T> T inc_beta(double const& a, T const& x) {
	double eps = 1e-30;
	if (cons(x) > 1-eps)
		return 1;
	else if (cons(x) < eps)
		return 0;

	T inc_beta;
	if ((int)(2*a)%2 == 0) { // d is impair
		unsigned int k = a;
		T b = 1;
		T sum_b = b;
		for (size_t i=2; i<k+1; i++) {
			b *= ((i-1)*2.0-1)/(2.0*(i-1))*x;
			if (abs_cons(b) < eps)
				break;
			sum_b += b;
		}
		return 1-sqrt(1-x)*sum_b;
	}
	else if ((int)(2*a)%2 == 1) { // d is pair
		unsigned int k = a + 0.5;
		T d = 2/PI;
		T sum_d = d;
		if (k == 1)
			sum_d = 0;
		for (size_t i=2; i<k; i++) {
			d *= ((i-1)*2.0)/(2.0*i-1)*x;
			if (abs_cons(d) < eps)
				break;
			sum_d += d;
		}
		return 2/PI*atan(sqrt(x/(1-x))) - sqrt((1-x)*x)*sum_d;
	}
	return 0;
}

// Cumulative Distribution Function of a chi squared distribution with d DoF.
// It is equal to the Regularised incomplete gamma-function
// Phi_d(r)=gamma_x(d/2,r/2)/Gamma(d/2). Computed by iteration.
// See: https://en.wikipedia.org/wiki/Incomplete_gamma_function#Evaluation_formulae
template<typename T> T chi_2_cdf(
	unsigned int const& d, T const& r) {
	T inc_gamma;

	// Safeguard
	if (d == 0)
		return 0.0;

	// Case disjunction
	T r_div_2_i = r/2.0;
	T exp_r = exp(-r_div_2_i);
	if (d%2 == 1) {
		unsigned int p = (d-1)/2;
		
		// Init inc_gamma(1/2, r/2)
		inc_gamma = sqrt(PI)*erf(sqrt(r_div_2_i));
		T r_div_2_i_pow = sqrt(r_div_2_i)*exp_r;

		// Exact iterative expression
		for (unsigned int i=0; i<p; i++) {
			inc_gamma = (i + 0.5)*inc_gamma - r_div_2_i_pow;
			if (i + 1 != p)
				r_div_2_i_pow *= r_div_2_i;
		}
	}
	else if (d%2==0) { // d is pair
		unsigned int p = d/2; 

		// Init inc_gamma(1, r/2)
		inc_gamma = 1.0 - exp_r;
		r_div_2_i *= exp_r;

		// Exact iterative expression
		for (unsigned int i=1; i<p; i++) {
			inc_gamma = i*inc_gamma - r_div_2_i;
			if (i + 1 != p)
				r_div_2_i *= r_div_2_i;
		}
	}

	// Normalise and return.
	return inc_gamma/std::tgamma((1.0*d)/2.0); 
}

// Inverse of the Cumulative Distribution Function of a chi squared distribution with d DoF.
double inv_chi_2_cdf(unsigned int const& d, double const& p);

// Computes the Cumulative Distribution Function of a Normal distribution over a crown of radii r_1 > r_2.
template<typename T> T delta_phi(unsigned int const& d, T const& r_1, T const& r_2) {
	return chi_2_cdf(d,r_1) - chi_2_cdf(d,r_2);
}

// Computes the conservatism.
// DOI: WIP
double conservatism(double const& beta, double const& beta_r);

// Computes the Mahalanobis distance from constraints to a normal distribution.
// DOI: WIP
template<typename T> DACE::AlgebraicVector<T>  get_list_distance(
	DACE::AlgebraicVector<T> const& mean, DACE::AlgebraicVector<T> const& diag_Sigma) {
	std::size_t d(mean.size());
	DACE::AlgebraicVector<T> list_distance(d);
	for (std::size_t i=0; i<d; i++) {
		if (cons(diag_Sigma[i]) > 0)
			list_distance[i] = -mean[i]/sqrt(diag_Sigma[i]);
		else if (cons(diag_Sigma[i]) == 0)
			list_distance[i] = -mean[i]*1e30;
	}
	return list_distance;
}

// Computes the first order failure risk estimation.
// DOI: WIP
template<typename T> T first_order_risk_estimation(
	DACE::AlgebraicVector<T> const& mean, DACE::AlgebraicVector<T> const& diag_Sigma) {
	// Get dimensions
	size_t d = mean.size();
	double a = (d-1.0)/2.0;

	// Null distribution.
	if (d==0) {
		return 0.0;
	}

	// Get and store list distances.
	DACE::AlgebraicVector<T> list_distance = get_list_distance(mean, diag_Sigma);
	list_distance = sort_vector(list_distance);
	if (cons(list_distance[0]) <= 0)
		return 1;
	return 1-chi_2_cdf(d, DACE::sqr(list_distance[0]));
}

// Computes the d-th order failure risk estimation.
// DOI: WIP
template<typename T> T dth_order_risk_estimation(
	DACE::AlgebraicVector<T> const& mean, DACE::AlgebraicVector<T> const& diag_Sigma) {

	// Get dimensions
	size_t d = mean.size();
	double a = (d-1.0)/2.0;

	// Null distribution.
	if (d==0)
		return 0.0;
	
	// Get and store list distances.
	DACE::AlgebraicVector<T> list_distance = get_list_distance(mean, diag_Sigma);
	list_distance = sort_vector(list_distance);
	if (cons(list_distance[0]) <= 0)
		return 1.0;

	// Compute the cdf for each layer.
	T beta_robust = 1 - chi_2_cdf(d, DACE::sqr(list_distance[0]));
	T layer = 1 - beta_robust;
	
	for (size_t i=1; i<d; i++) {
		if (abs_cons(layer) > 1 - EPS*EPS) { // Checks if there is still stuff to compute.
			break;
		}

		T r_i = list_distance[i];
		T r_im1 = list_distance[i-1];
		T delta_phi_i = delta_phi(d, DACE::sqr(r_i), DACE::sqr(r_im1)); // Total crown cdf.

		layer += delta_phi_i;
		if (abs_cons(delta_phi_i) > EPS*EPS) { // If it is worth computing...
			if (abs_cons(beta_robust) != 0) {
				T sum_sector = 0;
				for (size_t j=1; j<i; j++) {
					if (abs_cons(sum_sector) >= 2.0) { // The next terms will be useless.
						sum_sector = 2.0;
						break;
					}
					sum_sector += inc_beta(a, 1.0 - DACE::sqr(list_distance[j]/r_i));
				}
				if (abs_cons(sum_sector) == 2.0) // The next terms will be useless.
					break;
				else // Improve beta robust
					beta_robust -= delta_phi_i*(1-0.5*sum_sector);
			}
		}
	}
	return beta_robust;
}

#endif
