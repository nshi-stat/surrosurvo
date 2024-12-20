#include <Rcpp.h>

using namespace Rcpp;

extern DataFrame surrosurvoCpp(const DataFrame, const int);

extern DataFrame confintCpp(const DataFrame, const int, const int);

extern void sso(const NumericVector, const NumericVector,
                const NumericVector, const NumericVector,
                double *, double *, double *, double *, double *);

extern void jsso(const NumericVector, const NumericVector,
                 const NumericVector, const NumericVector,
                 int, double *, double *, double *, double *, double *);

extern NumericVector array2nvec(const double *, const int);
