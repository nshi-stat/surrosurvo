#include <Rcpp.h>

using namespace Rcpp;

extern DataFrame surrosurvoCpp(const DataFrame, const int, const int);

extern DataFrame confintCpp(const DataFrame, const int, const int, const int);

extern NumericVector jsso(const NumericVector, const NumericVector,
                 const NumericVector, const NumericVector,
                 const NumericVector, const NumericVector,
                 const int, const int, const int);
