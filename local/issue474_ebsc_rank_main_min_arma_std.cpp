
// The eBsc package on CRAN (at version 4.17) reports an error under the Intel compiler.
// This error now also appears with RcppArmadillo 15.0.0. In the example for plot.eBsc,
// the function drbasis() is called with n=250 and a second parameter p that varies from
// 1 to 6. In the case of 5, NA values return leading an rank failure. The drbasis()
// function calls an underlying C function drbasis from file drbasisC.cpp which has case
// statements for the differenct values of p.  We extrace the one for 5 here.

#if defined(USE_ARMA)
#include <armadillo>
#else
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>
#endif

std::vector<double> compute(int nn) {
    const double E = std::exp(1);
    std::vector<double> t(nn);
    std::vector<double> k(nn);
    typedef std::complex<double> dcomp;
    double Pi = M_PI;
    std::complex<double> comp0 (0,1);
    std::complex<double> comp1 (-1,0);
    std::complex<double> comp2 (0,-0.5);
    std::complex<double> comp3 (1,-1);
    std::complex<double> comp4 (-1,1);
    std::complex<double> comp5 (0,0.5);
    std::complex<double> comp6 (0,0.25);
    int const0 = 1;
    double const1 = -std::sqrt(3);
    double const2 = 2*std::sqrt(3);
    double const3 = -std::sqrt(5);
    double const4 = 6*std::sqrt(5);
    double const5 = -6*std::sqrt(5);
    double const6 = -std::sqrt(7);
    double const7 = 12*std::sqrt(7);
    double const8 = -30*std::sqrt(7);
    double const9 = 20*std::sqrt(7);
    double const10 = 3;
    double const11 = -60;
    double const12 = 270;
    double const13 = -420;
    double const14 = 210;
    /*
    double const15 = -std::sqrt(11);
    double const16 = 30*std::sqrt(11);
    double const17 = -210*std::sqrt(11);
    double const18 = 560*std::sqrt(11);
    double const19 = 630*std::sqrt(11);
    double const20 = 252*std::sqrt(11);
    */

    for (int j=0;j<nn;j++) {
        t[j]= j/double((nn-1));
    }
    for (int j=0;j<nn;j++) {
        k[j]= j+1;
    }

    std::vector<double> M5(nn);
    constexpr int i = 0;
    for (int j=0;j<nn;j++) {
        if (j==0) {
            M5[j] = const0;
        } else if (j==1) {
            M5[j] = const1 + const2*t[i];
        } else if (j==2) {
            M5[j] = const3 + const4*t[i] + const5* std::pow(t[i],2);
        } else if (j==3) {
            M5[j] = const6 + const7*t[i] + const8* std::pow(t[i],2) + const9* std::pow(t[i],3);
        } else if (j==4) {
            M5[j] = const10 + const11*t[i] + const12* std::pow(t[i],2) + const13* std::pow(t[i],3) + const14* std::pow(t[i],4);
        } else {
            dcomp ans = (std::sqrt(2)*(1 + std::sqrt(5)/2) -
                         (comp5*(std::sqrt(10 - 2*std::sqrt(5)) +
                                 std::sqrt(2*(5 + std::sqrt(5)))))/std::sqrt(2))*(std::pow(comp1,1 + k[j])/std::pow(E,std::pow(comp1,0.1)* (dcomp(-3 + k[j]))*Pi*(1 - t[i])) + std::pow(E,-(std::pow(comp1,0.1)*(dcomp(-3 + k[j]))* Pi*t[i]))) +
                (-(1/std::sqrt(2)) -
                 (comp5*(std::sqrt(10 - 2*std::sqrt(5)) + std::sqrt(2*(5 + std::sqrt(5)))))/
                 std::sqrt(2))*(std::pow(comp1,1 + k[j])/std::pow(E,std::pow(comp1,0.3)*(dcomp (-3 + k[j]))*Pi*(1 - t[i])) + std::pow(E,-(std::pow(comp1,0.3)*(dcomp (-3 + k[j]))*Pi*t[i]))) +
                (std::sqrt(2)*(1 + std::sqrt(5)/2) + (comp5*(std::sqrt(10 - 2*std::sqrt(5)) + std::sqrt(2*(5 + std::sqrt(5)))))/std::sqrt(2))
                *(std::pow(comp1,1 + k[j])/ std::pow(E,(std::sqrt(0.625 + std::sqrt(5)/8) - comp6*(-1 + std::sqrt(5)))*(dcomp (-3 + k[j]))*Pi*
                                           (1 - t[i])) +std::pow(E,-((std::sqrt(0.625 + std::sqrt(5)/8.) - comp6*(-1 + std::sqrt(5)))*(dcomp (-3 + k[j]))*
                                                                Pi*t[i]))) +(-(1/std::sqrt(2)) + (comp5*(std::sqrt(10 - 2*std::sqrt(5)) + std::sqrt(2*(5 + std::sqrt(5)))))/std::sqrt(2))
                *(std::pow(comp1,1 + k[j])/std::pow(E,(dcomp(std::sqrt(0.625 - std::sqrt(5)/8)) - comp6*(1 + std::sqrt(5)))*
                                          (dcomp (-3 + k[j]))*Pi*(1 - t[i])) +std::pow(E,-((std::sqrt(0.625 - std::sqrt(5)/8) - comp6*(1 + std::sqrt(5)))* (dcomp (-3 + k[j]))*Pi*t[i])))- std::sqrt(2)*cos((-3 + k[j])*Pi*t[i]);
            M5[j] =ans.real();
        }
        M5[j] = M5[j]/std::sqrt(nn);
    }
    return M5;
}

int main(int argc, char *argv[]) {
    constexpr int nn = 250;
    std::vector<double> M = compute(nn);
    for (int i=240; i<250; i++) std::cout << M[i] << " ";
    std::cout << std::endl;
    exit(0);
}
