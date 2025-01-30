/**
	robust_trajectory.h

	Purpose: Implementation of the RobustTrajectory class.

	@author Thomas Caleb

	@version 1.0 30/01/2025
*/

#ifndef DEF_ROBUST_TRAJECTORY
#define DEF_ROBUST_TRAJECTORY

#pragma once

#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include <dace/dace_s.h>

#include "settings.h"
#include "trajectory_split.h"


class RobustTrajectory : public std::deque<TrajectorySplit> {
private:
    std::deque<DACE::matrixdb> list_inv_covariance_;

public:
    // Default constructor
    RobustTrajectory() = default;

    // Constructor
    RobustTrajectory(std::deque<TrajectorySplit> const& other);

    // Copy constructor
    RobustTrajectory(RobustTrajectory const& other) = default;

    // Copy assignment operator
    RobustTrajectory& operator=(RobustTrajectory const& other) = default;

    // Destructor
    ~RobustTrajectory() = default;

    // Getter
    const std::deque<DACE::matrixdb> list_inv_covariance() const;

    // Setter
    void set_list_inv_covariance(std::deque<DACE::matrixdb> const& list_inv_covariance);

    // Returns the sorter list of all mahalanobis distances for a given vector x.
    const std::vector<std::pair<double, std::size_t>> get_mahalanobis_distance(DACE::vectordb const& vector) const;
};

#endif
