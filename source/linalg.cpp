/**
    linalg.cpp

    Purpose: Implementation of the linear algebra methods.

    @author Thomas Caleb

    @version 2.0 12/12/2024
*/

#include "linalg.h"

using namespace std;
using namespace DACE;


// Wraps a value between 0 and mod > 0.
double wrap_mod(double const& value, double const& mod) {
    double val = value;
    while (val > mod || val < 0.0) {
        if (val >= mod)
            val -= mod;
        else if (val < 0)
            val += mod;
    }
    return val;
}

// Wraps a value between 0 and mod > 0.
// DA version
DA wrap_mod(DACE::DA const& value, double const& mod) {
    DA val = value;
    double cons = val.cons();
    val += wrap_mod(cons, mod) - cons;
    return val;
}

// Turns a sym_tridiag_matrixdb to a matrixdb
// upperdiag = true (default) means the sub-diagonal is copied on the upper diagonal
// upperdiag = false means the upper diagonal is 0, therefore, the output is not symmetric
matrixdb sym_tridiag_matrixdb_2_matrixdb_(
    sym_tridiag_matrixdb const& tridiag, bool const& upperdiag) {
    // Unpack
    vector<matrixdb> list_D = tridiag.first;
    vector<matrixdb> list_A = tridiag.second;
    size_t size_diag = list_D.size();

    // Get matrix size
    size_t size_matrix = 0;
    for (size_t i = 0; i < size_diag; i++) { size_matrix += list_D[i].ncols(); }

    // Build matrix
    matrixdb output(size_matrix, size_matrix, 0.0);
    size_t counter = 0;
    for (size_t i = 0; i < size_diag; i++) {
        // Unpack
        matrixdb D_i = list_D[i];
        matrixdb A_i;
        size_t size_D_i = D_i.ncols();
        size_t nrows_A_i = 0;
        if (i != size_diag - 1) {
            A_i = list_A[i];
            nrows_A_i = A_i.nrows();
        }

        // Diagonal term
        for (size_t k = 0; k < size_D_i; k++) {
            for (size_t l = 0; l < size_D_i; l++) {
                output.at(counter + k, counter + l) = D_i.at(k, l);
            }
        }

        // Sub/sur Diag term
        if (i != size_diag - 1) {
            for (size_t k = 0; k < nrows_A_i; k++) {
                for (size_t l = 0; l < size_D_i; l++) {
                    double buffer = A_i.at(k, l);
                    output.at(counter + k + size_D_i, counter + l) = buffer;
                    if (upperdiag)
                        output.at(counter + l, counter + k + size_D_i) = buffer;
                }
            }
        }
        counter += size_D_i;
    }
    return output;
}

// Gets an identity matrix of size n.
matrixdb identity_(size_t const& n) {
    matrixdb A(n, n, 0.0);
    for (size_t i = 0; i < n; i++) { A.at(i, i) = 1.0; }
    return A;
}

// Determines the norm of A - diag(A).
double diagonal_error_(
    DACE::matrixdb const& A) {
    // Unpack and init
    size_t n = A.ncols();
    matrixdb A_test(A);

    // Assign and return
    for (size_t i = 0; i < n; i++) { A_test.at(i, i) = 0.0; }
    return frobenius_norm_(A_test);
}

// Returns a diagonal matrix equal to a given vector
matrixdb make_diag_matrix_(vectordb const& diag) {
    // Unpack and init
    size_t n = diag.size();
    matrixdb A(n, n, 0.0);

    // Assign and return
    for (size_t i = 0; i < n; i++) { A.at(i, i) = diag[i]; }
    return A;
}

// Turns a matrix to a line.
vectordb matrix_to_vector(matrixdb const& M) {
    size_t row(M.nrows()), col(M.ncols());
    vectordb output(row*col);
    for (size_t i=0; i<row; i++) {
        for (size_t j=0; j<col; j++) {
            output[i*col + j] = M.at(i,j);
        }
    }
    return output;
}

// Computes the trace of a matrix.
double trace_(matrixdb const& A) {
	double trace = 0.0;
	size_t n = A.ncols();
	for (size_t i = 0; i < n; i++) { trace += A.at(i, i); }
	return trace;
}

// Checks if a given symmetric matrix is positive-definite.
// Using the eigenvalue criterion: Sp(S) > 0 and S is sym <=> S is sym def pos
bool is_def_pos_(matrixdb const& S) {
    // Compute eigen values
    vectordb eig = jacobi_eigenvalue_(S).first;
    size_t n = eig.size();
    for (size_t i = 0; i < n; i++) {

        // Sp(S) > 0 S in sym <=> S is def pos
        if (eig[i] <=0) 
            return false;
    }
    return true;
}

// Computes the Frobenius norm of a matrix
// That is sqrt(trace(A * A^t))
double frobenius_norm_(matrixdb const& A) {
    // Unpack
    size_t n(A.nrows()), m(A.ncols());

	// Get trace(A * A^t)
    double trace = 0.0;
    for (size_t i = 0; i < n; i++) {
        for (size_t k = 0; k < m; k++) {
            trace += sqr(A.at(i, k));
        }
    }

	// Compute norm
	return sqrt(trace);
}

// Inverts a lower triangular matrix
matrixdb inv_traingluar_matrix_(
    matrixdb const& L) {
    // Unpack
    size_t n = L.nrows();

    // Init
    matrixdb inv_L(n, n);
    for (size_t i = 0; i < n; i++) {
        // Get i-th vector of the cannonical base
        vectordb buff_i(n, 0.0); buff_i[i] = 1.0;

        // Assign forward substitution
        inv_L.setcol(i, forward_substitution_(L, buff_i, i));
    }

    return inv_L;
}

// Computes the Cholesky factorization of
// a symetric positive-definite matrix.
// That is the unique Lower triangular matrix such that:
// S = L * L^t
// See https://en.wikipedia.org/wiki/Cholesky_decomposition
matrixdb cholesky_(matrixdb const& S) {
    size_t n = S.nrows();
    matrixdb L(n, n, 0.0);
    for (size_t i = 0; i < n; ++i) {
        for (size_t k = 0; k < i; ++k) {
            double buff = S.at(i, k);
            for (size_t j = 0; j < k; ++j)
                buff -= L.at(i, j) * L.at(k, j);
            L.at(i, k) = buff / L.at(k, k);
        }
        double buff = S.at(i, i);
        for (size_t j = 0; j < i; ++j)
            buff -= sqr(L.at(i, j));
        L.at(i, i) = sqrt(buff);
    }
    return L;
}

// Computes the Cholesky factorization of
// a symetric positive-definite matrix.
// That is the unique Lower triangular matrix such that:
// S = L * L^t
// See https://en.wikipedia.org/wiki/Cholesky_decomposition
// DA version.
matrixDA cholesky_(matrixDA const& S) {
    size_t n = S.nrows();
    matrixDA L(n, n, 0.0);
    for (size_t i = 0; i < n; ++i) {
        for (size_t k = 0; k < i; ++k) {
            DA buff = S.at(i, k);
            for (size_t j = 0; j < k; ++j)
                buff -= L.at(i, j) * L.at(k, j);
            L.at(i, k) = buff / L.at(k, k);
        }
        DA buff = S.at(i, i);
        for (size_t j = 0; j < i; ++j)
            buff -= sqr(L.at(i, j));
        if (buff.cons() > 0) {
            L.at(i, i) = sqrt(buff);
        }
    }
    return L;
}

// Computes the Cholesky factorization of
// a symetric positive-definite tridiagonal matrix.
// That is the unique Lower triangular matrix such that:
// S = L * L^t
// See https://en.wikipedia.org/wiki/Cholesky_decompositionsym_tridiag_matrix 
// This algorithm adapted to triadiagonal matrices was copied from [Cao et al. 2002]
// DOI: https://doi.org/10.1109/ICPPW.2002.1039748
sym_tridiag_matrixdb cholesky_(sym_tridiag_matrixdb const& tridiag) {
    // Unpack
    vector<matrixdb> list_diag = tridiag.first;
    vector<matrixdb> list_A = tridiag.second;
    size_t N = list_A.size();

    // Init
    vector<matrixdb> list_D, list_E;
    list_D.reserve(N + 1); list_E.reserve(N);
    for (size_t i = 0; i < N; i++) {
        // Unpack
        matrixdb diag_i = list_diag[i];
        matrixdb A_i = list_A[i];

        // Step
        matrixdb D_i = cholesky_(diag_i);
        matrixdb E_i = A_i * (inv_traingluar_matrix_(D_i)).transpose();

        // Update
        list_diag[i + 1] = list_diag[i + 1] - E_i * E_i.transpose();

        // Assign
        list_D.emplace_back(D_i);
        list_E.emplace_back(E_i);
    }

    // Last step
    list_D.emplace_back(cholesky_(list_diag[N]));

    return sym_tridiag_matrixdb(list_D, list_E);
}

// Perfoms lower triangular system solving.
vectordb forward_substitution_(
    matrixdb const& L, vectordb const& b,
    size_t const& starting_index) {
    size_t size_L = L.ncols();
    vectordb x(size_L);
    for (size_t i = starting_index; i < size_L; i++) {
        double buff_z_tilde = b[i];
        for (size_t j = starting_index; j < i; j++) {
            buff_z_tilde -= L.at(i,j) * x[j];
        }
        x[i] = buff_z_tilde / L.at(i, i);
    }
    return x;
}

// Perfoms upper triangular system solving 
// L is a lower triangular matrix, to avoid costly transposition
vectordb backward_substitution_(
    matrixdb const& L, vectordb const& b) {
    size_t size_L = L.ncols();
    vectordb x(size_L);
    for (size_t i = 0; i < size_L; i++) {
        size_t index_i = size_L - 1 - i;
        double buff_z = b[index_i];
        for (size_t j = 0; j < i; j++) {
            size_t index_j = size_L - 1 - j;
            buff_z -= L.at(index_j, index_i) * x[index_j];
        }
        x[index_i] = buff_z / L.at(index_i, index_i);
    }
    return x;
}

// Solves the system: L * L^t * x = b
// Where L is lower triangular, and b is a vector
vectordb solve_cholesky_(
    matrixdb const& L,
    vectordb const& b) {

    // Forward then backward substitution
    return backward_substitution_(L, forward_substitution_(L, b));
}

// Solves the system: L * L^t * X = B
// Where L is lower triangular, and B is a matrix
matrixdb solve_cholesky_( // matrix version
    matrixdb const& L,
    matrixdb const& B) {
    // Init
    size_t n = B.ncols();
    matrixdb X = B;
    vectordb b;

    // Use the vectorial version for each colomn of B
    for (size_t i = 0; i < n; i++) {
        b = solve_cholesky_(L, B.getcol(i));
        X.setcol(i, b);
    }

    return X;
}

// Solves the system: L * L^t * x = b
// Where L is lower triangular and tridiagonal, and b is a vector
vectordb solve_cholesky_(
    sym_tridiag_matrixdb const& tridiag_L,
    vectordb const& b) {
    // Unpack
    vector<matrixdb> list_D = tridiag_L.first;
    vector<matrixdb> list_A = tridiag_L.second;
    size_t N=list_A.size();
    size_t size_L=0;
    for (size_t i = 0; i < N + 1; i++) {size_L += list_D[i].nrows();}

    // The Forward and backward substitutions need to take into account
    // the subdiagonal terms

    // Top to bot step
    vectordb z_tilde(size_L, 0.0);
    size_t counter=0;
    for (size_t i = 0; i < N + 1; i++) {
        matrixdb D_i(list_D[i]), A_i;
        size_t size_D_i = D_i.nrows();
        size_t ncols_A_i;

        if (i != 0) {
            A_i = list_A[i - 1];
            ncols_A_i = A_i.ncols();
        }

        if (i==0) {
            for (size_t k = 0; k < size_D_i; k++) {
                double buff_z_tilde = b[counter + k];

                for (size_t l = 0; l < k; l++) {
                    buff_z_tilde -= D_i.at(k, l) * z_tilde[counter + l];
                }
                z_tilde[counter + k] = buff_z_tilde / D_i.at(k, k);
            }
        }
        else {
            for (size_t k = 0; k < size_D_i; k++) {
                double buff_z_tilde = b[counter + k];
                
                for (size_t l = 0; l < k; l++) { // D_i
                    buff_z_tilde -= D_i.at(k, l) * z_tilde[counter + l];
                }                
                for (size_t l = 0; l < ncols_A_i; l++) { // A_i
                    buff_z_tilde -= A_i.at(k, l) * z_tilde[counter - ncols_A_i + l];
                }

                z_tilde[counter + k] = buff_z_tilde / D_i.at(k, k);
            }
        }
        counter += size_D_i;
    }

    // Bot to Top step
    vectordb z(size_L, 0.0);
    for (size_t i = 0; i < N + 1; i++) {
        matrixdb D_i(list_D[N - i]), A_i;
        size_t size_D_i = D_i.nrows();
        size_t ncols_A_i;
        counter -= size_D_i;

        if (i != 0) {
            A_i = list_A[N - 1 - (i - 1)];
            ncols_A_i = A_i.ncols();
        }

        if (i == 0) {
            for (size_t k = 0; k < size_D_i; k++) {
                size_t index_k = size_D_i - 1 - k;
                vectordb col_D_i_k = D_i.getcol(index_k);
                double buff_z = z_tilde[counter + index_k];

                for (size_t l = 0; l < k; l++) {
                    size_t index_l = size_D_i - 1 - l;
                    buff_z -= col_D_i_k[index_l] * z[counter + index_l];
                }
                z[counter + index_k] = buff_z / col_D_i_k[index_k];
            }
        }
        else {
            for (size_t k = 0; k < size_D_i; k++) {
                size_t index_k = size_D_i - 1 - k;
                vectordb col_D_i_k = D_i.getcol(index_k);
                vectordb col_A_i_k = A_i.getcol(index_k);
                size_t size_col_A_i_k = col_A_i_k.size();
                double buff_z = z_tilde[counter + index_k];

                for (size_t l = 0; l < k; l++) { // D_i
                    size_t index_l = size_D_i - 1 - l;
                    buff_z -= col_D_i_k[index_l] * z[counter + index_l];
                }
                for (size_t l = 0; l < size_col_A_i_k; l++) { // A_i
                    size_t index_l = (size_col_A_i_k - 1) - l;
                    buff_z -= col_A_i_k[index_l] * z[counter + size_D_i + index_l];
                }
                z[counter + index_k] = buff_z / col_D_i_k[index_k];
            }
        }
    }
    
    return z;
}

// Jacobi eigenvalue algorithm
// Only for symetric matrices
// See https://en.wikipedia.org/wiki/Jacobi_eigenvalue_algorithm
pair<vectordb, matrixdb> jacobi_eigenvalue_(matrixdb const& S_) {
    // Unpack
    matrixdb S = S_;
    size_t n = S.ncols();

    // Trival cases n=0, n=1, it is trivial
    if (n == 0) // Null matrix
        return pair<vectordb, matrixdb>();
    else if (n == 1) // Scalar
        return pair<vectordb, matrixdb>(
            vectordb(n, S.at(0, 0)), matrixdb(n, n, 1.0));

    // Check if already diag
    double diag_criterion = diagonal_error_(S);
    if (diag_criterion < TOL_JAC) {
        vectordb eig(n);
        for (size_t i = 0; i < n; i++) {
            eig[i] = S.at(i, i);
        }
        return pair<vectordb, matrixdb>(eig, identity_(n));
    }

    // Init
    matrixdb V = identity_(n);
    vectordb d = get_diag_vector_(S);
    vectordb bw(d); vectordb zw(n, 0.0);

    // Loop
    size_t counter = 0;
    while (counter < MAX_ITER_JAC) {
        counter++;

        //  The convergence threshold is based on the size of the elements in
        //  the strict upper triangle of the matrix.
        double thresh = 0.0;
        for (size_t j = 0; j < n; j++) {
            for (size_t i = 0; i < j; i++) {
                thresh = thresh + S.at(j, i)* S.at(j, i);
            }
        }
        thresh = sqrt(thresh) / (4.0 * n);

        if (thresh == 0.0)
            break;

        for (size_t p = 0; p < n; p++) {
            for (size_t q = p + 1; q < n; q++) {
                double gapq = 10 * abs(S.at(q, p));
                double termp = gapq + abs(d[p]);
                double termq = gapq + abs(d[q]);
                
                //  Annihilate tiny offdiagonal elements.
                if (counter >= 4 && termp == abs(d[p]) && termq == abs(d[q]))
                    S.at(q, p) = 0.0;

                //  Otherwise, apply a rotation.
                else if (thresh <= abs(S.at(q, p))) {
                    double h = d[q] - d[p];
                    double term = abs(h) + gapq;
                    double t, g;
                    if (term == abs(h))
                        t = S.at(q, p) / h;
                    else {
                        double theta = 0.5 * h / S.at(q, p);
                        t = 1.0 / (abs(theta) + sqrt(1.0 + theta * theta));
                        if (theta < 0.0)
                            t = -t;
                    }
                    double c = 1.0 / sqrt(1.0 + t * t);
                    double s = t * c;
                    double tau = s / (1.0 + c);
                    h = t * S.at(q, p);

                    //  Accumulate corrections to diagonal elements.
                    zw[p] -= h;
                    zw[q] += h;
                    d[p] -= h;
                    d[q] += h;
                    S.at(q, p) = 0.0;

                    //  Rotate, using information from the upper triangle of S only.
                    for (size_t j = 0; j < p; j++) {
                        g = S.at(p, j);
                        h = S.at(q, j);
                        S.at(p, j) -= s * (h + g * tau);
                        S.at(q, j) += s * (g - h * tau);
                    }

                    for (size_t j = p + 1; j < q; j++) {
                        g = S.at(j, p);
                        h = S.at(q, j);
                        S.at(j, p) -= s * (h + g * tau);
                        S.at(q, j) += s * (g - h * tau);
                    }

                    for (size_t j = q + 1; j < n; j++) {
                        g = S.at(j, p);
                        h = S.at(j, q);
                        S.at(j, p) -= s * (h + g * tau);
                        S.at(j, q) += s * (g - h * tau);
                    }

                    //  Accumulate information in the eigenvector matrix.
                    for (size_t j = 0; j < n; j++) {
                        g = V.at(p, j);
                        h = V.at(q, j);
                        V.at(p, j) -= s * (h + g * tau);
                        V.at(q, j) += s * (g - h * tau);
                    }
                }
            }
        }

        for (size_t i = 0; i < n; i++) {
            bw[i] = bw[i] + zw[i];
            d[i] = bw[i];
            zw[i] = 0.0;
        }
    }

    if (counter >= MAX_ITER_JAC)
        cout << "Error : no convergence" << endl;
    
    return pair<vectordb, matrixdb>(d, V);
}

// Jacobi eigenvalue algorithm
// Only for symetric matrices
// See https://en.wikipedia.org/wiki/Jacobi_eigenvalue_algorithm
// DA version
pair<vectorDA, matrixDA> jacobi_eigenvalue_(matrixDA const& S_) {
    // Unpack
    matrixDA S = S_;
    size_t n = S.ncols();

    // Trival cases n=0, n=1, it is trivial
    if (n == 0) // Null matrix
        return pair<vectorDA, matrixDA>();
    else if (n == 1) // Scalar
        return pair<vectorDA, matrixDA>(
            vectorDA(n, S.at(0, 0)), matrixDA(n, n, 1.0));

    // Check if already diag
    double diag_criterion = diagonal_error_(S.cons());
    if (diag_criterion < TOL_JAC) {
        vectorDA eig(n);
        for (size_t i = 0; i < n; i++) {
            eig[i] = S.at(i, i);
        }
        return pair<vectorDA, matrixDA>(eig, identity_(n) + 0*DA(1));
    }

    // Init
    matrixDA V = identity_(n) + 0*DA(1);
    vectorDA d = get_diag_vector_(S);
    vectorDA bw(d); vectorDA zw(n, 0.0);

    // Loop
    size_t counter = 0;
    while (counter < MAX_ITER_JAC) {
        counter++;

        //  The convergence threshold is based on the size of the elements in
        //  the strict upper triangle of the matrix.
        double thresh = 0.0;
        for (size_t j = 0; j < n; j++) {
            for (size_t i = 0; i < j; i++) {
                thresh = thresh + sqr(S.at(j, i).cons());
            }
        }
        thresh = sqrt(thresh) / (4.0 * n);

        if (thresh == 0.0)
            break;

        for (size_t p = 0; p < n; p++) {
            for (size_t q = p + 1; q < n; q++) {
                double gapq = 10 * abs_cons(S.at(q, p));
                double termp = gapq + abs_cons(d[p]);
                double termq = gapq + abs_cons(d[q]);

                //  Annihilate tiny offdiagonal elements.
                if (counter >= 4 && termp == abs_cons(d[p]) && termq == abs_cons(d[q]))
                    S.at(q, p) = 0.0;

                //  Otherwise, apply a rotation.
                else if (thresh <= abs_cons(S.at(q, p))) {
                    DA h = d[q] - d[p];
                    double term = abs_cons(h) + gapq;
                    DA t, g;
                    if (term == abs_cons(h)) {
                        t = S.at(q, p) / h;
                    }
                    else {
                        DA theta = 0.5 * h / S.at(q, p);
                        DA abs_theta = theta;
                        if (theta.cons() < 0)
                            abs_theta += 2*abs_cons(theta);
                        t = 1.0 / (abs_theta + sqrt(1.0 + sqr(theta))); // DA
                        if (theta.cons() < 0.0)
                            t = -t;
                    }
                    DA c = 1.0 / sqrt(1.0 + sqr(t));
                    DA s = t * c;
                    DA tau = s / (1.0 + c);
                    h = t * S.at(q, p);

                    //  Accumulate corrections to diagonal elements.
                    zw[p] -= h;
                    zw[q] += h;
                    d[p] -= h;
                    d[q] += h;
                    S.at(q, p) = 0.0;

                    //  Rotate, using information from the upper triangle of S only.
                    for (size_t j = 0; j < p; j++) {
                        g = S.at(p, j);
                        h = S.at(q, j);
                        S.at(p, j) -= s * (h + g * tau);
                        S.at(q, j) += s * (g - h * tau);
                    }

                    for (size_t j = p + 1; j < q; j++) {
                        g = S.at(j, p);
                        h = S.at(q, j);
                        S.at(j, p) -= s * (h + g * tau);
                        S.at(q, j) += s * (g - h * tau);
                    }

                    for (size_t j = q + 1; j < n; j++) {
                        g = S.at(j, p);
                        h = S.at(j, q);
                        S.at(j, p) -= s * (h + g * tau);
                        S.at(j, q) += s * (g - h * tau);
                    }

                    //  Accumulate information in the eigenvector matrix.
                    for (size_t j = 0; j < n; j++) {
                        g = V.at(p, j);
                        h = V.at(q, j);
                        V.at(p, j) -= s * (h + g * tau);
                        V.at(q, j) += s * (g - h * tau);
                    }
                }
            }
        }

        for (size_t i = 0; i < n; i++) {
            bw[i] = bw[i] + zw[i];
            d[i] = bw[i];
            zw[i] = 0.0;
        }
    }
    if (counter >= MAX_ITER_JAC)
        cout << "Error : no convergence" << endl;
    return pair<vectorDA, matrixDA>(d, V);
}
