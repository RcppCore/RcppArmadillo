
// The eBsc package on CRAN (at version 4.17) reports an error under the Intel compiler.
// This error now also appears with RcppArmadillo 15.0.0. In the example for plot.eBsc,
// the function drbasis() is called with n=250 and a second parameter p that varies from
// 1 to 6. In the case of 5, NA values return leading an rank failure. The drbasis()
// function calls an underlying C function drbasis from file drbasisC.cpp which has case
// statements for the differenct values of p.  We extrace the one for 5 here.

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo/RcppArmadillo>

// [[Rcpp::export]]
arma::mat compute(int nn) {
    using namespace arma;
    using namespace std;
    //Rcpp::List out;
    const double E = exp(1);
    arma::vec t(nn);
    arma::vec k(nn);
    typedef complex<double> dcomp;
    double Pi = M_PI;
    complex<double> comp0 (0,1);
    complex<double> comp1 (-1,0);
    complex<double> comp2 (0,-0.5);
    complex<double> comp3 (1,-1);
    complex<double> comp4 (-1,1);
    complex<double> comp5 (0,0.5);
    complex<double> comp6 (0,0.25);
    int const0 = 1;
    double const1 = -sqrt(3);
    double const2 = 2*sqrt(3);
    double const3 = -sqrt(5);
    double const4 = 6*sqrt(5);
    double const5 = -6*sqrt(5);
    double const6 = -sqrt(7);
    double const7 = 12*sqrt(7);
    double const8 = -30*sqrt(7);
    double const9 = 20*sqrt(7);
    double const10 = 3;
    double const11 = -60;
    double const12 = 270;
    double const13 = -420;
    double const14 = 210;
    double const15 = -sqrt(11);
    double const16 = 30*sqrt(11);
    double const17 = -210*sqrt(11);
    double const18 = 560*sqrt(11);
    double const19 = 630*sqrt(11);
    double const20 = 252*sqrt(11);

    for (int j=0;j<nn;j++) {
        t(j)= j/double((nn-1));
    }
    for (int j=0;j<nn;j++) {
        k(j)= j+1;
    }


    arma::mat M5(nn,nn);
    for (int i=0;i<nn; i++) {
        for (int j=0;j<nn;j++) {
            if (j==0) {
                M5(i,j)=const0;
            } else if (j==1) {
                M5(i,j)= const1 + const2*t(i);
            } else if (j==2) {
                M5(i,j)= const3 + const4*t(i) + const5* pow(t(i),2);
            } else if (j==3) {
                M5(i,j)= const6 + const7*t(i) + const8* pow(t(i),2) + const9* pow(t(i),3);
            } else if (j==4) {
                M5(i,j)= const10 + const11*t(i) + const12* pow(t(i),2) + const13* pow(t(i),3) + const14* pow(t(i),4);
            } else {
                dcomp ans = (sqrt(2)*(1 + sqrt(5)/2) -
                             (comp5*(sqrt(10 - 2*sqrt(5)) +
                                     sqrt(2*(5 +
                                             sqrt(5)))))/sqrt(2))*(pow(comp1,1 + k(j))/pow(E,pow(comp1,0.1)* (dcomp(-3 + k(j)))*Pi*(1 - t(i))) +
                                                                   pow(E,-(pow(comp1,0.1)*(dcomp(-3 + k(j)))* Pi*t(i)))) +
                    (-(1/sqrt(2)) -
                     (comp5*(sqrt(10 - 2*sqrt(5)) + sqrt(2*(5 + sqrt(5)))))/
                     sqrt(2))*(pow(comp1,1 + k(j))/pow(E,pow(comp1,0.3)*(dcomp (-3 + k(j)))*Pi*(1 - t(i))) + pow(E,-(pow(comp1,0.3)*(dcomp (-3 + k(j)))*Pi*t(i)))) +
                    (sqrt(2)*(1 + sqrt(5)/2) + (comp5*(sqrt(10 - 2*sqrt(5)) + sqrt(2*(5 + sqrt(5)))))/sqrt(2))
                    *(pow(comp1,1 + k(j))/ pow(E,(sqrt(0.625 + sqrt(5)/8) - comp6*(-1 + sqrt(5)))*(dcomp (-3 + k(j)))*Pi*
                                               (1 - t(i))) +pow(E,-((sqrt(0.625 + sqrt(5)/8.) - comp6*(-1 + sqrt(5)))*(dcomp (-3 + k(j)))*
                                                                    Pi*t(i)))) +(-(1/sqrt(2)) + (comp5*(sqrt(10 - 2*sqrt(5)) + sqrt(2*(5 + sqrt(5)))))/sqrt(2))
                    *(pow(comp1,1 + k(j))/pow(E,(dcomp(sqrt(0.625 - sqrt(5)/8)) - comp6*(1 + sqrt(5)))*
                                              (dcomp (-3 + k(j)))*Pi*(1 - t(i))) +pow(E,-((sqrt(0.625 - sqrt(5)/8) - comp6*(1 + sqrt(5)))* (dcomp (-3 + k(j)))*Pi*t(i))))- sqrt(2)*cos((-3 + k(j))*Pi*t(i));
                M5(i,j) =ans.real();
            }
            M5(i,j) = M5(i,j)/sqrt(nn);
        }
    }
    /*
    arma::mat vectors=M5;
    arma::mat vectorsQR=Turn(M5); arma::vec values=eigenvalues(nn,5);
    out = Rcpp::List::create(Rcpp::Named("eigenvectors") = vectors,
                             Rcpp::Named("eigenvectorsQR")   = vectorsQR,
                             Rcpp::Named("eigenvalues") = values,
                             Rcpp::Named("x") = t) ;
    return out;
    */
    return M5;
}

/*** R
res <- compute(250)
if (!any(is.na(res))) {
  cat("Kappa: ")
  kappa(res)
} else {
  cat("Nb NAs: ")
  sum(is.na(res))
  print(res[1:10,241:250])
}
*/
