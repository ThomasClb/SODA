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

// Returns the probability of ending in that branch.
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

// Checks if a splitting history is a child of this one.
const int SplittingHistory::is_child(
	SplittingHistory const& history_i) const {
	int n(this->size());
	int n_i(history_i.size());
	SplittingHistory hist;
	SplittingHistory hist_i;

	for (int i=0; i<min(n_i, n); i++) {
		// Extract
		hist.assign(this->begin(), this->begin() + i + 1);
		hist_i.assign(history_i.begin(), history_i.begin() + i + 1);
		if (!(hist.can_merge(hist_i))) {
			return n + n_i - 2*i;
		}
	}
	return n + n_i - 2*min(n_i, n);
}

// Print
const string SplittingHistory::to_string() const {
	ostringstream oss;
    for (const auto& entry : *this) {
        oss << entry.first << ", " << entry.second << "; ";
    }
    return oss.str();
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