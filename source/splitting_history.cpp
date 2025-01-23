/**
	splitting_history.cpp

	Purpose: Implementation of the SplittingHistory class.

	@author Thomas Caleb

	@version 1.0 23/01/2025
*/

#include "splitting_history.h"

using namespace std;

// Default constructor
SplittingHistory::SplittingHistory() = default;

// Copy constructor
SplittingHistory::SplittingHistory(SplittingHistory const& splitting_history) = default;

// Copy assignment operator
SplittingHistory& SplittingHistory::operator=(SplittingHistory const& other) = default;

// Destructor
SplittingHistory::~SplittingHistory() = default;

const double SplittingHistory::alpha() const {
	double alpha = 1.0;
	for (size_t i=0; i<this->size(); i++) {
		int dir = this->at(i).second;
		if (dir == 0)
			alpha *= ALPHA_0_GMM;
		else if (dir*dir == 1)
			alpha *= ALPHA_1_GMM;
	}
	return alpha;
}

// Checks if two splitting histories can be merged.
const bool SplittingHistory::can_merge(
	SplittingHistory const& other) const {
	size_t n(this->size());
	if (n != other.size())
		return false;
	for (size_t i=0; i<n-1; i++) {
		if (!(this->at(i).first == other[i].first 
			&& this->at(i).second == other[i].second))
			return false;
	}
	if (this->at(n - 1).first != other[n - 1].first)
		return false;
	if (this->at(n - 1).second == other[n - 1].second)
		return false;
	return true;
}