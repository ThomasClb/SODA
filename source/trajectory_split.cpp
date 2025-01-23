/**
	trajectory_split.cpp

	Purpose: Implementation of the TrajectorySplit class.

	@author Thomas Caleb

	@version 1.0 23/01/2025
*/

#include "trajectory_split.h"

using namespace DACE;
using namespace std;

// Default constructor
TrajectorySplit::TrajectorySplit()  = default;

// Constructor with parameters
TrajectorySplit::TrajectorySplit(
	vector<statedb> const& list_x,
    vector<controldb> const& list_u,
    SplittingHistory const& splitting_history)
    : list_x_(list_x), list_u_(list_u), splitting_history_(splitting_history) {}

// Copy constructor
TrajectorySplit::TrajectorySplit(
	TrajectorySplit const& other) = default;

// Copy assignment operator
TrajectorySplit& TrajectorySplit::operator=(TrajectorySplit const& other) = default;

// Destructor
TrajectorySplit::~TrajectorySplit() = default;

// Getters
const vector<statedb> TrajectorySplit::list_x() const {
    return list_x_;
}

const vector<controldb> TrajectorySplit::list_u() const {
    return list_u_;
}

const SplittingHistory TrajectorySplit::splitting_history() const {
    return splitting_history_;
}

// Setters
void TrajectorySplit::set_list_x(vector<statedb> const& list_x) {
    list_x_ = list_x;
}

void TrajectorySplit::set_list_u(vector<controldb> const& list_u) {
    list_u_ = list_u;
}

void TrajectorySplit::set_splitting_history(SplittingHistory const& splitting_history) {
    splitting_history_ = splitting_history;
}
