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
SplittingHistory::SplittingHistory(SplittingHistory const& history) = default;

// Copy assignment operator
SplittingHistory& SplittingHistory::operator=(SplittingHistory const& history) = default;

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
	SplittingHistory const& history) const {
	size_t n(this->size());
	if (n != history.size())
		return false;
	for (size_t i=0; i<n-1; i++) {
		if (!(this->at(i).first == history[i].first 
			&& this->at(i).second == history[i].second))
			return false;
	}
	if (this->at(n - 1).first != history[n - 1].first)
		return false;
	if (this->at(n - 1).second == history[n - 1].second)
		return false;
	return true;
}

// Overload the << operator for output
ostream& operator<<(ostream& os, const SplittingHistory& history) {
    for (const auto& entry : history) {
        os << entry.first << ", " << entry.second << "; ";
    }
    cout << endl;
    return os;
}

// Overload the >> operator for input
std::istream& operator>>(std::istream& is, SplittingHistory& history) {
    history.clear(); // Clear the history before reading new data
    unsigned int first;
    int second;
    char comma, semicolon;
    while (is >> first >> comma >> second >> semicolon) {
        if (comma == ',' && semicolon == ';') {
            history.emplace_back(first, second);
        } else {
            is.setstate(std::ios::failbit);
            break;
        }
    }
    return is;
}