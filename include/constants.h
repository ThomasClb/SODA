/**
	constants.h

	Purpose: Implementation of the Constants class.

	@author Thomas Caleb

	@version 2.0 02/09/2024
*/

#ifndef DEF_CONSTANTS
#define DEF_CONSTANTS

#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <dace/dace_s.h>

#include "settings.h"


/*

	Constants

*/
class Constants {

// Attributes
protected:
	double mu_; // Mu normalisation unit [km^3.s^-2]
	double lu_; // Lenghts normalisation unit [km]
	double wu_; // Pulsation normalisation unit [s^-1]
	double massu_; // Spacecraft mass normalisation unit [kg]
	double tu_; // Time normalisation unit [s]
	double vu_; // Velocity normalisation unit [m.s^-1]
	double thrustu_; // Spacecraft thrust normalisation unit [m.s^-2]

// Methods
public:
	// Constructors

	// Default constructors
	Constants();

	// Constructor
	Constants(
		double const& mu,
		double const& lu,
		double const& wu,
		double const& massu);

	// Copy constructor
	Constants(
		Constants const& constants);

	// Destructors
	~Constants();

	// Getters
	const double mu() const;
	const double lu() const;
	const double wu() const;
	const double massu() const;
	const double tu() const;
	const double vu() const;
	const double thrustu() const;

	// IO operator
	friend std::ostream& operator<<(std::ostream& os, const Constants& constants);
	friend std::istream& operator>>(std::istream& is, Constants& constants);
};

#endif
