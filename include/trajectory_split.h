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

#include "ddp_solver.h"
#include "settings.h"
#include "constants.h"
#include "state_t.h"
#include "control_t.h"
#include "loads.h"
#include "splitting_history.h"

class DDPSolver;


class TrajectorySplit {
private:
    std::vector<statedb> list_x_;
    std::vector<controldb> list_u_;
    SplittingHistory splitting_history_;
    std::vector<DACE::vectorDA> list_dynamics_eval_;

public:
    // Default constructor
    TrajectorySplit();

    // Constructor with parameters
    TrajectorySplit(
    	std::vector<statedb> const& list_x,
    	std::vector<controldb> const& list_u,
    	SplittingHistory const& splitting_history);
    TrajectorySplit(
        std::vector<statedb> const& list_x,
        std::vector<controldb> const& list_u,
        SplittingHistory const& splitting_history,
        std::vector<DACE::vectorDA> const& list_dynamics_eval);

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
    const std::vector<DACE::vectorDA> list_dynamics_eval() const; 

    // Setters 
    void set_list_x(std::vector<statedb> const& list_x);
    void set_list_u(std::vector<controldb> const& list_u);
    void set_splitting_history(SplittingHistory const& splitting_history);
    void set_list_dynamics_eval(std::vector<DACE::vectorDA> const& list_dynamics_eval);

    // Find the splitting direction of highest nonlinearity.
    unsigned int find_splitting_direction(
        std::size_t const& index, double const& transcription_beta);

    // Return an updated splitted trajectory.
    TrajectorySplit get_splited_trajectory(
        DACE::vectordb const& modified_x0, DACE::matrixdb const& modified_Sigma,
        DDPSolver const& DDPsolver);

    // Split a TrajectorySplit into 3 along a given direction.
    std::pair<TrajectorySplit,TrajectorySplit> split(
        unsigned int const& dir, DDPSolver const& DDPsolver);

    // Merge 3 TrajectorySplit into 1 along a given direction.
    // The merged split is the central one.
    void merge(
        unsigned int const& dir, DDPSolver const& DDPsolver);


    // Friend functions for stream operators
    friend std::ostream& operator<<(std::ostream& os, const TrajectorySplit& ts);
    friend std::istream& operator>>(std::istream& is, TrajectorySplit& ts);
};

#endif
