/**
	test_trajectory_split.cpp

	Purpose: Test of the implementation of the TrajectorySplit class.

	@author Thomas Caleb

	@version 1.0 23/01/2025
*/

#include "gtest/gtest.h"
#include "../include/soda.h"

using namespace DACE;
using namespace std;


TEST(TestTrajectorySplit, ConstructorGetSet) {
    // Default constructor
    TrajectorySplit ts_default;
    EXPECT_EQ(ts_default.list_x().size(), 0);
    EXPECT_EQ(ts_default.list_u().size(), 0);
    EXPECT_EQ(ts_default.splitting_history().size(), 0);

    // Constructor with parameters
    std::vector<statedb> list_x = {statedb(), statedb()};
    std::vector<controldb> list_u = {controldb(), controldb()};
    SplittingHistory splitting_history;
    splitting_history.push_back({1, 100});

    TrajectorySplit ts_params(list_x, list_u, splitting_history);
    EXPECT_EQ(ts_params.list_x().size(), list_x.size());
    EXPECT_EQ(ts_params.list_u().size(), list_u.size());
    EXPECT_EQ(ts_params.splitting_history().size(), splitting_history.size());

    // Copy constructor
    TrajectorySplit ts_copy(ts_params);
    EXPECT_EQ(ts_copy.list_x().size(), ts_params.list_x().size());
    EXPECT_EQ(ts_copy.list_u().size(), ts_params.list_u().size());
    EXPECT_EQ(ts_copy.splitting_history().size(), ts_params.splitting_history().size());

    // Copy assignment operator
    TrajectorySplit ts_assigned;
    ts_assigned = ts_params;
    EXPECT_EQ(ts_assigned.list_x().size(), ts_params.list_x().size());
    EXPECT_EQ(ts_assigned.list_u().size(), ts_params.list_u().size());
    EXPECT_EQ(ts_assigned.splitting_history().size(), ts_params.splitting_history().size());

    // Setters
    std::vector<statedb> new_list_x = {statedb()};
    std::vector<controldb> new_list_u = {controldb()};
    SplittingHistory new_splitting_history;
    new_splitting_history.push_back({2, 200});

    ts_default.set_list_x(new_list_x);
    ts_default.set_list_u(new_list_u);
    ts_default.set_splitting_history(new_splitting_history);

    EXPECT_EQ(ts_default.list_x().size(), new_list_x.size());
    EXPECT_EQ(ts_default.list_u().size(), new_list_u.size());
    EXPECT_EQ(ts_default.splitting_history().size(), new_splitting_history.size());
}