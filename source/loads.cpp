/**
    loads.cpp

    Purpose: Implementation of the Low-Order Automatic Domain Splitting methods.

    @author Thomas Caleb

    @version 1.0 23/10/2024
*/

#include "loads.h"

using namespace std;
using namespace DACE;

// Shortcut function to retrieve eigenvalues knowing robust trajectory.
pair<vectordb, matrixdb> get_eigenvalues(matrixdb const& Sigma_x, matrixdb const& feedback_gain) {
     // Diagonalize Sigma

    // Retrieve state eigenvalues and eigenvectors
    pair<vectordb, matrixdb> diag_x = jacobi_eigenvalue_(Sigma_x);
    vectordb lambda_x = diag_x.first;
    matrixdb V_x = diag_x.second;

    // Retrieve control eigenvalues and eigenvectors
    pair<vectordb, matrixdb> diag_u = jacobi_eigenvalue_(
        feedback_gain*Sigma_x*feedback_gain.transpose());
    vectordb lambda_u = diag_u.first;
    matrixdb V_u = diag_u.second;

    // Merge into a large matrix of size Nu + Nx
    size_t d_x(lambda_x.size()), d_u(lambda_u.size());
    size_t d(d_x + d_u);
    vectordb lambda(d, 0.0);
    matrixdb V(d, d, 0.0);
    for (size_t i=0; i<d_x; i++) {
        lambda[i] = lambda_x[i];
        for (size_t j=0; j<d_x; j++) { V.at(i,j) = V_x.at(i,j); }
    }
    for (size_t i=0; i<d_u; i++) {
        lambda[i + d_x] = lambda_u[i];
        for (size_t j=0; j<d_u; j++) { V.at(i + d_x, j + d_x) = V_u.at(i,j); }
    }

    pair<vectordb, matrixdb> output(lambda, V);
    return output;
}

// Scales a vector based on a robust trajectory.
pair<vectorDA, vectordb> scale(
    vectorDA const& y, matrixdb const& Sigma_x,
    matrixdb const& feedback_gain, double const& transcription_beta) {
    // Diagonalize Sigma
    pair<vectordb, matrixdb> eigen = get_eigenvalues(Sigma_x, feedback_gain);
    vectordb lambda(inv_chi_2_cdf(feedback_gain.ncols(), 1 - transcription_beta)*eigen.first);
    matrixdb V(eigen.second);

    // Make identity vector
    size_t d(feedback_gain.ncols() + feedback_gain.nrows());
    vectorDA dx(d, 0.0);
    for (unsigned int i=0; i<d; i++) {
        if (lambda[i] > 0) {
            dx[i] = sqrt(lambda[i])*DA(i + 1);
        }
    }
    vectorDA dx_scaled = V*dx;
    
    // Scale y
    vectorDA y_scaled = y.eval(dx_scaled);

    return pair<vectorDA, vectordb>(y_scaled, lambda);
}

// Return the GMM decomposition of a given multivariate Gaussian
// Distribution along a specified direction.
// From [DeMars et al. 2013]
// DOI: https://doi.org/10.2514/1.58987
vector<tuple<double, vectordb, matrixdb>> gmm(
    vectordb const& y_mean,
    matrixdb const& Sigma,
    size_t const& direction) {

    // Compute eigenvalues
    pair<vectordb, matrixdb> eig(jacobi_eigenvalue_(Sigma));
    vectordb eigenvalues(eig.first);
    matrixdb eigenvectors(eig.second);

    // Compute Sigma
    eigenvalues[direction] *= SIGMA_GMM*SIGMA_GMM;
    matrixdb Sigma_tilde(eigenvectors*make_diag_matrix_(eigenvalues)*eigenvectors.transpose());

    // Compute nominal_state_tilde
    vector<tuple<double, vectordb, matrixdb>> output(3);
    vectordb d_nominal_state_tilde = (sqrt(eigenvalues[direction])*MU_GMM)*vectordb(eigenvectors.getcol(direction));
    output[0] = {ALPHA_1_GMM, y_mean - d_nominal_state_tilde, Sigma_tilde};
    output[1] = tuple<double,vectordb,matrixdb>{ALPHA_0_GMM, y_mean, Sigma_tilde};
    output[2] = tuple<double,vectordb,matrixdb>{ALPHA_1_GMM, y_mean + d_nominal_state_tilde, Sigma_tilde};

    return output;
}

// Computes the NonLinearity Index (NLI) of a vector given a scaling.
// From [Losacco et al. 2024]
// DOI: https://doi.org/10.2514/1.G007271
double nl_index(vectorDA const& y, vectordb const& lambda) {
    size_t d(DA::getMaxVariables());
    matrixDA J(d, y.size()); // DA expansion of the Jacobian
    for (size_t i=0; i<y.size(); i++) {
        J.setcol(i, y[i].gradient());
        for (size_t j=0; j<d; j++) { // Scaling
            if (lambda[j] > 0)
                J.at(j,i) /= sqrt(lambda[j]);
        }
    }
    matrixdb J_cons(J.cons()); // Linear part (STM)

    matrixdb B(J_cons.nrows(), J_cons.ncols(), 0.0);
    for (size_t i=0; i<J_cons.nrows(); i++) {
        for (size_t j=0; j<J_cons.ncols(); j++) {
            vectordb b(J.at(i, j).linear());
            double buff(0.0);
            for (size_t k=0; k<b.size(); k++) {
                buff += abs(b[k]);
            }
            B.at(i, j) = buff;
        }  
    }
    return frobenius_norm_(B)/frobenius_norm_(J_cons);
}

// Computes the NonLinearity Index (NLI) of a vector given a scaling along each direction.
// From [Losacco et al. 2024]
// DOI: https://doi.org/10.2514/1.G007271
vectordb nl_index_dir(vectorDA const& y, vectordb const& lambda) {
    size_t d(DA::getMaxVariables());
    matrixDA J(d, y.size()); // DA expansion of the Jacobian
    for (size_t i=0; i<y.size(); i++) {
        J.setcol(i, y[i].gradient());
        for (size_t j=0; j<d; j++) { // Scaling
            if (lambda[j] > 0)
                J.at(j,i) /= sqrt(lambda[j]);
        }
    }
    matrixdb J_cons(J.cons()); // Linear part (STM)
    double norm_J(frobenius_norm_(J_cons));

    // Loop on directions
    matrixdb B(J_cons.nrows(), J_cons.ncols(), 0.0);
    vectordb output(d); vectordb b;
    for (size_t dir=0; dir<d; dir++) {
        for (size_t i=0; i<J_cons.nrows(); i++) {
            for (size_t j=0; j<J_cons.ncols(); j++) {
                B.at(i, j) = J.at(i, j).linear()[dir];
            }  
        }
        output[dir] = frobenius_norm_(B)/norm_J;
    }
    return output;
}

// Computes the NonLinearity Index (NLI) of a vector given a robust trajectory step.
// From [Losacco et al. 2024]
// DOI: https://doi.org/10.2514/1.G007271
double nl_index(
    vectorDA const& y, matrixdb const& Sigma_x,
    matrixdb const& feedback_gain, double const& transcription_beta) {    
    // Normalize
    pair<vectorDA, vectordb> scaled_vector(scale(y, Sigma_x, feedback_gain, transcription_beta));
    return nl_index(scaled_vector.first, scaled_vector.second);
}

// Computes the NonLinearity Index (NLI) of a vector given a robust trajectory.
// From [Losacco et al. 2024]
// DOI: https://doi.org/10.2514/1.G007271
vectordb nl_index(
    vector<vectorDA> const& list_dynamics_eval, vector<statedb> const& list_x,
    vector<controldb> const& list_u, double const& transcription_beta) {

    // Init
    vectordb output(list_dynamics_eval.size());

    // First scaling
    pair<vectorDA, vectordb> scaled_vector = scale(
        list_dynamics_eval[0], list_x[0].Sigma(), list_u[0].feedback_gain(), transcription_beta);
    vectorDA x_k = scaled_vector.first;
    output[0] = nl_index(x_k, scaled_vector.second);

    // Loop on all
    for (size_t i=1; i<list_dynamics_eval.size(); i++) {
        // Get correction
        matrixdb K(list_u[i].feedback_gain());
        vectorDA dx_k(x_k - x_k.cons());
        matrixDA dx_k_mat(dx_k.size(), 1); dx_k_mat.setcol(0, dx_k);
        matrixDA du_k_mat(K*dx_k_mat);
        vectorDA delta(dx_k.size() + du_k_mat.nrows());
        for (size_t j=0; j<dx_k.size(); j++) {
            delta[j] = dx_k[j];
        }
        for (size_t j=0; j<du_k_mat.nrows(); j++) {
            delta[j + dx_k.size()] = du_k_mat.at(j, 0);
        }

        // Scale
        x_k = list_dynamics_eval[i].eval(delta); 

        // Diagonalize Sigma
        pair<vectordb, matrixdb> eig(get_eigenvalues(list_x[i].Sigma(), K));
        output[i] = nl_index(x_k, eig.first);
    }
    return output;
}
