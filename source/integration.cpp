/**
    integration.cpp

    Purpose: Implementation of the integration
    methods.

    @author Thomas Caleb

    @version 1.0 16/11/2023
*/

#include "integration.h"

using namespace std;
using namespace DACE;

// Declare integration constants
const vector<double> get_A() {
    std::vector<double> A(13);
    A[0] = 0.0; A[1] = 1.0 / 18.0; A[2] = 1.0 / 12.0;
    A[3] = 1.0 / 8.0; A[4] = 5.0 / 16.0; A[5] = 3.0 / 8.0;
    A[6] = 59.0 / 400.0; A[7] = 93.0 / 200.0; A[8] = 5490023248.0 / 9719169821.0;
    A[9] = 13.0 / 20.0; A[10] = 1201146811.0 / 1299019798.0; A[11] = 1.0;
    A[12] = 1.0;
    return A;
}
const vector<double> get_C() {
    std::vector<double> C(13);
    C[0] = 14005451.0 / 335480064.0; C[1] = 0.0; C[2] = 0.0;
    C[3] = 0.0; C[4] = 0.0; C[5] = -59238493.0 / 1068277825.0;
    C[6] = 181606767.0 / 758867731.0; C[7] = 561292985.0 / 797845732.0; C[8] = -1041891430.0 / 1371343529.0;
    C[9] = 760417239.0 / 1151165299.0; C[10] = 118820643.0 / 751138087.0; C[11] = -528747749.0 / 2220607170.0;
    C[12] = 1.0 / 4.0;
    return C;
}
const vector<double> get_D() {
    std::vector<double> D(13);
    D[0] = 13451932.0 / 455176623.0; D[1] = 0.0; D[2] = 0.0;
    D[3] = 0.0; D[4] = 0.0; D[5] = -808719846.0 / 976000145.0;
    D[6] = 1757004468.0 / 5645159321.0; D[7] = 656045339.0 / 265891186.0; D[8] = -3867574721.0 / 1518517206.0;
    D[9] = 465885868.0 / 322736535.0; D[10] = 53011238.0 / 667516719.0; D[11] = 2.0 / 45.0;
    D[12] = 0.0;
    return D;
}
const matrixdb get_B() {
    DACE::matrixdb B(13, 12);
    B.at(0, 0) = 0.0; B.at(0, 1) = 0.0; B.at(0, 2) = 0.0; B.at(0, 3) = 0.0; B.at(0, 4) = 0.0;
    B.at(0, 5) = 0.0; B.at(0, 6) = 0.0; B.at(0, 7) = 0.0; B.at(0, 8) = 0.0; B.at(0, 9) = 0.0;
    B.at(0, 10) = 0.0; B.at(0, 11) = 0.0;

    B.at(1, 0) = 1.0 / 18.0; B.at(1, 1) = 0.0; B.at(1, 2) = 0.0; B.at(1, 3) = 0.0; B.at(1, 4) = 0.0;
    B.at(1, 5) = 0.0; B.at(1, 6) = 0.0; B.at(1, 7) = 0.0; B.at(1, 8) = 0.0; B.at(1, 9) = 0.0;
    B.at(1, 10) = 0.0; B.at(1, 11) = 0.0;

    B.at(2, 0) = 1.0 / 48.0; B.at(2, 1) = 1.0 / 16.0; B.at(2, 2) = 0.0; B.at(2, 3) = 0.0; B.at(2, 4) = 0.0;
    B.at(2, 5) = 0.0; B.at(2, 6) = 0.0; B.at(2, 7) = 0.0; B.at(2, 8) = 0.0; B.at(2, 9) = 0.0;
    B.at(2, 10) = 0.0; B.at(2, 11) = 0.0;

    B.at(3, 0) = 1.0 / 32.0; B.at(3, 1) = 0.0; B.at(3, 2) = 3.0 / 32.0; B.at(3, 3) = 0.0; B.at(3, 4) = 0.0;
    B.at(3, 5) = 0.0; B.at(3, 6) = 0.0; B.at(3, 7) = 0.0; B.at(3, 8) = 0.0; B.at(3, 9) = 0.0;
    B.at(3, 10) = 0.0; B.at(3, 11) = 0.0;

    B.at(4, 0) = 5.0 / 16.0; B.at(4, 1) = 0.0; B.at(4, 2) = -75.0 / 64.0; B.at(4, 3) = 75.0 / 64.0; B.at(4, 4) = 0.0;
    B.at(4, 5) = 0.0; B.at(4, 6) = 0.0; B.at(4, 7) = 0.0; B.at(4, 8) = 0.0; B.at(4, 9) = 0.0;
    B.at(4, 10) = 0.0; B.at(4, 11) = 0.0;

    B.at(5, 0) = 3.0 / 80.0; B.at(5, 1) = 0.0; B.at(5, 2) = 0.0; B.at(5, 3) = 3.0 / 16.0; B.at(5, 4) = 3.0 / 20.0;
    B.at(5, 5) = 0.0; B.at(5, 6) = 0.0; B.at(5, 7) = 0.0; B.at(5, 8) = 0.0; B.at(5, 9) = 0.0;
    B.at(5, 10) = 0.0; B.at(5, 11) = 0.0;

    B.at(6, 0) = 29443841.0 / 614563906.0; B.at(6, 1) = 0.0; B.at(6, 2) = 0.0; B.at(6, 3) = 77736538.0 / 692538347.0; B.at(6, 4) = -28693883.0 / 1125000000.0;
    B.at(6, 5) = 23124283.0 / 1800000000.0; B.at(6, 6) = 0.0; B.at(6, 7) = 0.0; B.at(6, 8) = 0.0; B.at(6, 9) = 0.0;
    B.at(6, 10) = 0.0; B.at(6, 11) = 0.0;

    B.at(7, 0) = 16016141.0 / 946692911.0; B.at(7, 1) = 0.0; B.at(7, 2) = 0.0; B.at(7, 3) = 61564180.0 / 158732637.0; B.at(7, 4) = 22789713.0 / 633445777.0;
    B.at(7, 5) = 545815736.0 / 2771057229.0; B.at(7, 6) = -180193667.0 / 1043307555.0; B.at(7, 7) = 0.0; B.at(7, 8) = 0.0; B.at(7, 9) = 0.0;
    B.at(7, 10) = 0.0; B.at(7, 11) = 0.0;

    B.at(8, 0) = 39632708.0 / 573591083.0; B.at(8, 1) = 0.0; B.at(8, 2) = 0.0; B.at(8, 3) = -433636366.0 / 683701615.0; B.at(8, 4) = -421739975.0 / 2616292301.0;
    B.at(8, 5) = 100302831.0 / 723423059.0; B.at(8, 6) = 790204164.0 / 839813087.0; B.at(8, 7) = 800635310.0 / 3783071287.0; B.at(8, 8) = 0.0; B.at(8, 9) = 0.0;
    B.at(8, 10) = 0.0; B.at(8, 11) = 0.0;

    B.at(9, 0) = 246121993.0 / 1340847787.0; B.at(9, 1) = 0.0; B.at(9, 2) = 0.0; B.at(9, 3) = -37695042795.0 / 15268766246.0; B.at(9, 4) = -309121744.0 / 1061227803.0;
    B.at(9, 5) = -12992083.0 / 490766935.0; B.at(9, 6) = 6005943493.0 / 2108947869.0; B.at(9, 7) = 393006217.0 / 1396673457.0; B.at(9, 8) = 123872331.0 / 1001029789.0; B.at(9, 9) = 0.0;
    B.at(9, 10) = 0.0; B.at(9, 11) = 0.0;

    B.at(10, 0) = -1028468189.0 / 846180014.0; B.at(10, 1) = 0.0; B.at(10, 2) = 0.0; B.at(10, 3) = 8478235783.0 / 508512852.0; B.at(10, 4) = 1311729495.0 / 1432422823.0;
    B.at(10, 5) = -10304129995.0 / 1701304382.0; B.at(10, 6) = -48777925059.0 / 3047939560.0; B.at(10, 7) = 15336726248.0 / 1032824649.0; B.at(10, 8) = -45442868181.0 / 3398467696.0; B.at(10, 9) = 3065993473.0 / 597172653.0;
    B.at(10, 10) = 0.0; B.at(10, 11) = 0.0;

    B.at(11, 0) = 185892177.0 / 718116043.0; B.at(11, 1) = 0.0; B.at(11, 2) = 0.0; B.at(11, 3) = -3185094517.0 / 667107341.0; B.at(11, 4) = -477755414.0 / 1098053517.0;
    B.at(11, 5) = -703635378.0 / 230739211.0; B.at(11, 6) = 5731566787.0 / 1027545527.0; B.at(11, 7) = 5232866602.0 / 850066563.0; B.at(11, 8) = -4093664535.0 / 808688257.0; B.at(11, 9) = 3962137247.0 / 1805957418.0;
    B.at(11, 10) = 65686358.0 / 487910083.0; B.at(11, 11) = 0.0;

    B.at(12, 0) = 403863854.0 / 491063109.0; B.at(12, 1) = 0.0; B.at(12, 2) = 0.0; B.at(12, 3) = -5068492393.0 / 434740067.0; B.at(12, 4) = -411421997.0 / 543043805.0;
    B.at(12, 5) = 652783627.0 / 914296604.0; B.at(12, 6) = 11173962825.0 / 925320556.0; B.at(12, 7) = -13158990841.0 / 6184727034.0; B.at(12, 8) = 3936647629.0 / 1978049680.0; B.at(12, 9) = -160528059.0 / 685178525.0;
    B.at(12, 10) = 248638103.0 / 1413531060.0; B.at(12, 11) = 0.0;
    return B;
}

// Utilities

// Computes the minimum
double min(double const& a, double const& b) {

    double res = a;
    if (a > b) { res = b; }
    return res;
}

// Computes the maximum
double max(double const& a, double const& b) {

    double res = a;
    if (b > a) { res = b; }
    return res;
}

// Computes the infinite norm
double normtmp(int const& N, vectordb const& X) {

    double res = 0.0;
    for (size_t I = 0; I < N; I++) {
        res = max(res, abs(X[I]));
    }

    return res;
}

// Integrators

// Propagates a state using a order 7(8) Runge - Kutta - Fehlberg embedded method(RK78)
// See [Fehlberg 1968] Classical fifth-, sixth-, seventh-, and eighthorder
// Runge - Kutta formulas with stepsize control.
// The integrated function takes a DA state, a DA control, a time and a SpacecraftParameters
// as inputs and provides a DA derivative of the same size.
// Adapted from: https://github.com/dacelib/dace/blob/master/Tutorials/Tutorial1/Example11DARK78.cpp
vectorDA RK78(
    vectorDA(*f)(
        vectorDA const&, vectorDA const&,
        double const&,
        SpacecraftParameters const&,
        Constants const&,
        SolverParameters const&),
    vectorDA Y0, vectorDA const& U,
    double const& X0, double const& DX,
    SpacecraftParameters const& spacecraft_parameters,
    Constants const& constants,
    SolverParameters const& solver_parameters) {

    double X1 = X0 + DX;
    int N = Y0.size();
    double ERREST;
    int MAX_ITER = 10000; int NB_ITER = 0;
    double H0 = MIN_STEP/constants.tu(); double HS = 0.1; double H1 = 100.0;
    double BS = 20 * EPS;

    matrixDA Z(N, 16);
    vectorDA Y1(N);
    std::vector<double> Y1cons(N);

    double VIHMAX = 0.0, X, H;
    int I, J, K;
    int sign = 1; if (X0 > X1) { sign = -1; }
    double RFNORM, HH0, HH1;
    double HSQR = 1.0 / 9.0;
    
    for (I = 0; I < N; I++) { Z.at(I, 0) = Y0[I]; Z.at(I, 1) = 0.0; }

    H = sign*abs(HS); HH0 = sign*abs(H0); HH1 = sign*abs(H1);
    X = X0; RFNORM = 0.0; ERREST = 0.0;
    while (X != X1 && NB_ITER < MAX_ITER) {
        NB_ITER++;

        // compute new stepsize
        if (RFNORM != 0) { H = H * min(4.0, exp(HSQR * log(EPS / RFNORM))); }
        if (abs(H) > abs(HH1)) { H = HH1; }
        else if (abs(H) < abs(HH0) * 0.99) {
            H = HH0;
        }

        if ((X + H - X1) * H > 0) { H = X1 - X; }

        for (J = 0; J < 13; J++) {

            for (I = 0; I < N; I++) {

                Y0[I] = 0.0; // EVALUATE RHS AT 13 POINTS

                for (K = 0; K < J; K++) { Y0[I] = Y0[I] + Z.at(I, K + 3) * B.at(J, K); }

                Y0[I] = H * Y0[I] + Z.at(I, 0);
            }

            Y1 = f(Y0, U, X + H * A[J], spacecraft_parameters, constants, solver_parameters);

            for (I = 0; I < N; I++) { Z.at(I, J + 3) = Y1[I]; }
        }

        for (I = 0; I < N; I++) {

            Z.at(I, 1) = 0.0; Z.at(I, 2) = 0.0; // EXECUTE 7TH,8TH ORDER STEPS

            for (J = 0; J < 13; J++) {
                Z.at(I, 1) = Z.at(I, 1) + Z.at(I, J + 3) * D[J];
                Z.at(I, 2) = Z.at(I, 2) + Z.at(I, J + 3) * C[J];
            }

            Y1[I] = (Z.at(I, 2) - Z.at(I, 1)) * H;
            Z.at(I, 2) = Z.at(I, 2) * H + Z.at(I, 0);
        }


        for (I = 0; I < N; I++) { Y1cons[I] = cons(Y1[I]); }

        RFNORM = normtmp(N, Y1cons); // ESTIMATE ERROR AND DECIDE ABOUT BACKSTEP

        if ((RFNORM > BS) && (abs(H / H0) > 1.2)) {
            H = H / 3.0;
            RFNORM = 0;
        }
        else {
            for (I = 0; I < N; I++) { Z.at(I, 0) = Z.at(I, 2); }
            X = X + H;
            VIHMAX = max(VIHMAX, sign*H);
            ERREST = ERREST + RFNORM;
        }
    }

    for (I = 0; I < N; I++) { Y1[I] = Z.at(I, 0); }

    return Y1;
}

// Propagates a state using a order 7(8) Runge - Kutta - Fehlberg embedded method(RK78)
// See [Fehlberg 1968] Classical fifth-, sixth-, seventh-, and eighthorder
// Runge - Kutta formulas with stepsize control.
// The integrated function takes a double state, a double control, a time and a SpacecraftParameters
// as inputs and provides a double derivative of the same size.
// Adapted from: https://github.com/dacelib/dace/blob/master/Tutorials/Tutorial1/Example11DARK78.cpp
vectordb RK78(
    vectordb(*f)(
        vectordb const&, vectordb const&,
        double const&,
        SpacecraftParameters const&,
        Constants const&,
        SolverParameters const&),
    vectordb Y0, vectordb const& U,
    double const& X0, double const& DX,
    SpacecraftParameters const& spacecraft_parameters,
    Constants const& constants,
    SolverParameters const& solver_parameters) {

    double X1 = X0 + DX;
    int N = Y0.size();
    double ERREST;
    int MAX_ITER = 10000; int NB_ITER = 0;
    double H0 = MIN_STEP / constants.tu(); double HS = 0.1; double H1 = 100.0;
    double BS = 20 * EPS;

    matrixdb Z(N, 16);
    vectordb Y1(N);

    double VIHMAX = 0.0, X, H;
    int I, J, K;
    int sign = 1; if (X0 > X1) { sign = -1; }
    double RFNORM, HH0, HH1;

    double HSQR = 1.0 / 9.0;
    for (I = 0; I < N; I++)
    {
        Z.at(I, 0) = Y0[I];
        Z.at(I, 1) = 0.0;
    }

    H = sign * abs(HS); HH0 = sign * abs(H0); HH1 = sign * abs(H1);
    X = X0; RFNORM = 0.0; ERREST = 0.0;
    while (X != X1 && NB_ITER < MAX_ITER) {

        NB_ITER++;

        // compute new stepsize
        if (RFNORM != 0) { H = H * min(4.0, exp(HSQR * log(EPS / RFNORM))); }
        if (abs(H) > abs(HH1)) { H = HH1; }
        else if (abs(H) < abs(HH0) * 0.99) {
            H = HH0;
        }

        if ((X + H - X1) * H > 0) { H = X1 - X; }

        for (J = 0; J < 13; J++) {

            for (I = 0; I < N; I++) {

                Y0[I] = 0.0; // EVALUATE RHS AT 13 POINTS

                for (K = 0; K < J; K++) { Y0[I] = Y0[I] + Z.at(I, K + 3) * B.at(J, K); }

                Y0[I] = H * Y0[I] + Z.at(I, 0);
            }

            Y1 = f(Y0, U, X + H * A[J], spacecraft_parameters, constants, solver_parameters);

            for (I = 0; I < N; I++) { Z.at(I, J + 3) = Y1[I]; }
        }

        for (I = 0; I < N; I++) {

            Z.at(I, 1) = 0.0; Z.at(I, 2) = 0.0; // EXECUTE 7TH,8TH ORDER STEPS

            for (J = 0; J < 13; J++) {
                Z.at(I, 1) = Z.at(I, 1) + Z.at(I, J + 3) * D[J];
                Z.at(I, 2) = Z.at(I, 2) + Z.at(I, J + 3) * C[J];
            }

            Y1[I] = (Z.at(I, 2) - Z.at(I, 1)) * H;
            Z.at(I, 2) = Z.at(I, 2) * H + Z.at(I, 0);
        }


        RFNORM = normtmp(N, Y1); // ESTIMATE ERROR AND DECIDE ABOUT BACKSTEP

        if ((RFNORM > BS) && (abs(H / H0) > 1.2)) {
            H = H / 3.0;
            RFNORM = 0;
        }
        else {
            for (I = 0; I < N; I++) { Z.at(I, 0) = Z.at(I, 2); }
            X = X + H;
            VIHMAX = max(VIHMAX, sign*H);
            ERREST = ERREST + RFNORM;
        }
    }

    for (I = 0; I < N; I++) { Y1[I] = Z.at(I, 0); }

    return Y1;
}
