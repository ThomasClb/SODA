/**
	settings.h

	Purpose: Declaration of useful constants and parameters.

	@author Thomas Caleb

	@version 2.0 02/09/2024
*/

#ifndef DEF_SETTINGS
#define DEF_SETTINGS

#pragma once

// Define useful constants
#define PI (4.0*atan(1.0)) // Definiton of pi [rad]
#define DEG_2_RAD (PI/180.0) // Definiton of pi [rad]
#define SIZE_VECTOR 6 // Size of a state vector (3 positions, 3 velocities)
#define SEC2DAYS (1.0/(24.*3600)) // Conversion from seconds to days [days.s^-1]

// Astronomical parameters

// From JPL DE431 ephemerides, publicaly available at:
// https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/gm_de431.tpc
#define MU_SUN (1.327124400419393e11) // Mass parameter of the Sun [km^3.s^-2]
#define MU_EARTH_MOON (4.035032355022598e5) // Mass parameter of the Earth-Moon system [km^3.s^-2]
#define MU_EARTH (3.9860043543609598e5) // Mass parameter of the Earth [km^3.s^-2]
#define MU_MOON (4.9028000661637961e3) // Mass parameter of the Moon [km^3.s^-2]

// Distance between the Sun and the Earth, and astronomical unit [km] 
// From Wikipedia for the AU definition
// https://en.wikipedia.org/wiki/Astronomical_unit
#define SUN_EARTH_DISTANCE (149597870.7)

// Distance between the Earth and the Moon [km]
// Same as in [Caleb et al. 2023]
// DOI: https://doi.org/10.1007/s11071-023-08375-0
#define EARTH_MOON_DISTANCE (384399)

// J2 implementation [-]
// https://en.wikipedia.org/wiki/Geopotential_model
#define J_2 (1.082635854e-3)

// Earth Radius [km]
// From Wikipedia
// https://en.wikipedia.org/wiki/Earth_radius
#define R_EARTH (6378)

// Moon Radius [km]
// From Wikipedia
// https://en.wikipedia.org/wiki/Moon
#define R_MOON (1736)

// Thrust
// https://en.wikipedia.org/wiki/Standard_gravity
#define G_0 (9.81) // Standard gravity [m.s^-1]

// Integration constants
#define MIN_STEP (1e-6) // Minimum step size [s]
#define EPS (1e-13) // Tolerance of the intergration scheme [-]

// Eigenvalue computation constants
#define MAX_ITER_JAC 1000 // Maximum number of iteration of the  matrix sqrt computation algorithm [-]
#define TOL_JAC (1000*EPS) // Tolerance of the matrix sqrt computation algorithm [-]

// Transcription
#define TRANSCRIPTION_METHOD 1 // 0 is spectral radius, 1 is first order.

// GMM from K=3 lambda=1e-3 See [DeMars et al. 2013]
// DOI: https://doi.org/10.2514/1.58987
#define SIGMA_GMM 0.6715664864669252 // Scaling [-]
#define ALPHA_0_GMM 0.5495506294920584 // Central weight [-]
#define ALPHA_1_GMM 0.225224685253970 // Lateral weight [-]
#define MU_GMM 1.057515048576096 // Offset [-]

#endif
