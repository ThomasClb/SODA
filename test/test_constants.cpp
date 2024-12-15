/**
	test_constants.cpp

	Purpose: Test of the implementation of the AULSolver class.

	@author Thomas Caleb

	@version 2.0 12/12/2024
*/

#include "gtest/gtest.h"
#include "../include/soda.h"

using namespace DACE;
using namespace std;

/*

	CONSTANTS

*/

// Constructors
TEST(TestConstants, EmptyConstructor) {
	// Init
	Constants constants;
	double mu = MU_SUN; // [km^3.s^-2]
	double lu = SUN_EARTH_DISTANCE; // [km]
	double wu = sqrt((MU_SUN) / pow(SUN_EARTH_DISTANCE, 3)); // [s^-1]
	double tu = 1/wu; // [s]
	double vu = lu*wu; // [m.s^-1]
	double massu = 1000; // [kg]
	double thrustu = 1000 * vu * massu * wu; // [N]

	// Tests
	EXPECT_EQ(constants.mu(), mu);
	EXPECT_EQ(constants.lu(), lu);
	EXPECT_EQ(constants.wu(), wu);
	EXPECT_EQ(constants.tu(), tu);
	EXPECT_EQ(constants.vu(), vu);
	EXPECT_EQ(constants.massu(), massu);
	EXPECT_EQ(constants.thrustu(), thrustu);
}
TEST(TestConstants, FilledConstructor) {
	// Init
	double mu = MU_MOON/(MU_EARTH + MU_MOON); // [-]
	double lu = EARTH_MOON_DISTANCE; // [km]
	double wu = sqrt((MU_EARTH + MU_MOON) / pow(EARTH_MOON_DISTANCE, 3)); // [s^-1]
	double tu = 1 / wu; // [s]
	double vu = lu * wu; // [m.s^-1]
	double massu = 1000; // [kg]
	double thrustu = 1000 * vu * massu * wu; // [N]
	Constants constants(mu, lu, wu, massu);

	// Tests
	EXPECT_EQ(constants.mu(), mu);
	EXPECT_EQ(constants.lu(), lu);
	EXPECT_EQ(constants.wu(), wu);
	EXPECT_EQ(constants.tu(), tu);
	EXPECT_EQ(constants.vu(), vu);
	EXPECT_EQ(constants.massu(), massu);
	EXPECT_EQ(constants.thrustu(), thrustu);
}
TEST(TestConstants, CopyConstructor) {
	// Init
	double mu = MU_MOON / (MU_EARTH + MU_MOON); // [-]
	double lu = EARTH_MOON_DISTANCE; // [km]
	double wu = sqrt((MU_EARTH + MU_MOON) / pow(EARTH_MOON_DISTANCE, 3)); // [s^-1]
	double tu = 1 / wu; // [s]
	double vu = lu * wu; // [m.s^-1]
	double massu = 1000; // [kg]
	double thrustu = 1000 * vu * massu * wu; // [N]
	Constants constants(mu, lu, wu, massu);
	Constants constants_copy = constants;

	// Tests
	EXPECT_EQ(constants_copy.mu(), mu);
	EXPECT_EQ(constants_copy.lu(), lu);
	EXPECT_EQ(constants_copy.wu(), wu);
	EXPECT_EQ(constants_copy.massu(), massu);
	EXPECT_EQ(constants_copy.tu(), tu);
	EXPECT_EQ(constants_copy.vu(), vu);
	EXPECT_EQ(constants_copy.thrustu(), thrustu);
}
TEST(TestConstants, IOFunctions) {
	// Init
	double mu = MU_SUN; // [km^3/s^2]
	double lu = SUN_EARTH_DISTANCE; // [km]
	double wu = sqrt((MU_SUN) / pow(SUN_EARTH_DISTANCE, 3)); // [s^-1]
	double tu = 1 / wu; // [s]
	double vu = lu * wu; // [m.s^-1]
	double massu = 1000; // [kg]
	double thrustu = 1000 * vu * massu * wu; // [N]
	Constants constants(mu, lu, wu, massu);
	Constants constants_copy; // empty

	// Open file
	string file_name_("../data/constants/test_constants.dat");
	ofstream ofs(file_name_);

	// Store the object to file
	ofs << constants;

	ofs.close();

	// Open file
	ifstream ifs(file_name_);

	// Load data
	ifs >> constants_copy;

	ifs.close();

	// Tests
	EXPECT_EQ(constants_copy.mu(), mu);
	EXPECT_EQ(constants_copy.lu(), lu);
	EXPECT_EQ(constants_copy.wu(), wu);
	EXPECT_EQ(constants_copy.massu(), massu);
	EXPECT_EQ(constants_copy.tu(), tu);
	EXPECT_EQ(constants_copy.vu(), vu);
	EXPECT_EQ(constants_copy.thrustu(), thrustu);
}
