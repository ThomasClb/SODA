/**
	trajectory_split.h

	Purpose: Implementation of the TrajectorySplit class.

	@author Thomas Caleb

	@version 1.0 23/01/2025
*/

#ifndef DEF_TRAJECTORY_SPLIT
#define DEF_TRAJECTORY_SPLIT

#pragma once

#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include <dace/dace_s.h>

#include "settings.h"
#include "constants.h"
#include "state_t.h"
#include "control_t.h"
#include "splitting_history.h"


class TrajectorySplit {
private:
    std::vector<statedb> list_x_;
    std::vector<controldb> list_u_;
    SplittingHistory splitting_history_;

public:
    // Default constructor
    TrajectorySplit();

    // Constructor with parameters
    TrajectorySplit(
    	std::vector<statedb> const& list_x,
    	std::vector<controldb> const& list_u,
    	SplittingHistory const& splitting_history);

    // Copy constructor
    TrajectorySplit(TrajectorySplit const& other);

    // Copy assignment operator
    TrajectorySplit& operator=(TrajectorySplit const& other);

    // Destructor
    ~TrajectorySplit();

    // Getters
    const std::vector<statedb> list_x() const;
    const std::vector<controldb> list_u() const;
    const SplittingHistory splitting_history() const;

    // Setters
    void set_list_x(std::vector<statedb> const& list_x);
    void set_list_u(std::vector<controldb> const& list_u);
    void set_splitting_history(SplittingHistory const& splitting_history);

    // Splits a trajectory along a given direction
    std::vector<TrajectorySplit> split(std::size_t const& direction);
};

#endif
