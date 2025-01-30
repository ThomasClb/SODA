/**
	robust_trajectory.cpp

	Purpose: Implementation of the RobustTrajectory class.

	@author Thomas Caleb

	@version 1.0 23/01/2025
*/

#include "robust_trajectory.h"

using namespace DACE;
using namespace std;

// Constructor
RobustTrajectory::RobustTrajectory(deque<TrajectorySplit> const& other) : list_inv_covariance_() {
    this->assign(other.begin(), other.end());
    double reg(1e-16);
    for (size_t i=0; i<this->size(); i++) {
        matrixdb Sigma(this->at(i).list_x()[0].Sigma());

        // Compute eigenvalues
        pair<vectordb, matrixdb> eig(jacobi_eigenvalue_(Sigma));
        vectordb eigenvalues(eig.first);
        matrixdb eigenvectors(eig.second);
        for (size_t k=0; k<eigenvalues.size(); k++) {
            double buff = eigenvalues[k];
            if (buff > 0) {
                eigenvalues[k] = 1/buff;
            }
        }
        matrixdb inv_Sigma(eigenvectors*make_diag_matrix_(eigenvalues)*eigenvectors.transpose());
        list_inv_covariance_.push_back(inv_Sigma);
    }
}

// Getter
const deque<DACE::matrixdb> RobustTrajectory::list_inv_covariance() const {
    return list_inv_covariance_;
}

// Setter
void RobustTrajectory::set_list_inv_covariance(deque<matrixdb> const& list_inv_covariance) {
    list_inv_covariance_ = list_inv_covariance;
}

// Returns the sorter list of all mahalanobis distances for a given vector x.
const vector<pair<double, size_t>> RobustTrajectory::get_mahalanobis_distance(vectordb const& x) const {

    vector<pair<double, size_t>> output(this->size());
    for (size_t i=0; i<this->size(); i++) {
        vectordb x_ref(this->at(i).list_x()[0].nominal_state());
        x_ref = x - x_ref;
        output[i].first = sqrt(x_ref.dot(list_inv_covariance_[i]*x_ref));
        output[i].second = i;
    }
    sort(output.begin(), output.end());
    return output;
}