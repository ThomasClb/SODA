/**
	splitting_history.h

	Purpose: Implementation of the SplittingHistory class.

	@author Thomas Caleb

	@version 1.0 23/01/2025
*/

#ifndef DEF_SPLITTING_HISTORY
#define DEF_SPLITTING_HISTORY

#pragma once

#include <vector>
#include <utility>

#include "settings.h"

class SplittingHistory : public std::vector<std::pair<unsigned int, int>> {
public:
    // Default constructor
    SplittingHistory();

    // Copy constructor
    SplittingHistory(SplittingHistory const& splitting_history);

    // Copy assignment operator
    SplittingHistory& operator=(SplittingHistory const& other);

    // Destructor
    ~SplittingHistory();

    // Getters
    const double alpha() const;

    // Checks if two splitting histories can be merged.
    const bool can_merge(SplittingHistory const& other) const;
};

#endif