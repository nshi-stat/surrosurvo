#include <Rcpp.h>

using namespace Rcpp;

extern DataFrame surrosurvoCpp(const DataFrame, const int, const int);

extern DataFrame confintCpp(const DataFrame, const int, const int, const int);

extern void sso(const NumericVector, const NumericVector,
                const NumericVector, const NumericVector,
                const NumericVector, const NumericVector,
                const int, const int,
                double *, double *, double *, double *);

extern void jsso(const NumericVector, const NumericVector,
                 const NumericVector, const NumericVector,
                 const NumericVector, const NumericVector,
                 const int, const int, const int,
                 double *, double *, double *, double *);

extern NumericVector array2nvec(const double *, const int);
