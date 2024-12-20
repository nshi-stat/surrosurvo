#include <Rcpp.h>
#include <Rmath.h>
#include "surrosurvo.h"

#ifdef _OPENMP
#include <omp.h>
// [[Rcpp::plugins("openmp")]]
#endif

using namespace Rcpp;

// [[Rcpp::export]]
DataFrame surrosurvoCpp(const DataFrame df, const int n) {

  double *tauo = new double[1];
  double *taumo1 = new double[1];
  double *taumo2 = new double[1];
  double *tauso1 = new double[1];
  double *tauso2 = new double[1];

  NumericVector y = df[0];
  NumericVector c = df[1];
  NumericVector x = df[2];
  NumericVector p = df[3];
  sso(y, c, x, p, &tauo[0], &taumo1[0], &taumo2[0], &tauso1[0], &tauso2[0]);

  DataFrame out = DataFrame::create(
      Named("tauo") = array2nvec(tauo, 1),
      Named("taumo1") = array2nvec(taumo1, 1),
      Named("taumo2") = array2nvec(taumo2, 1),
      Named("tauso1") = array2nvec(tauso1, 1),
      Named("tauso2") = array2nvec(tauso2, 1)
  );

  delete [] tauo;
  delete [] taumo1;
  delete [] taumo2;
  delete [] tauso1;
  delete [] tauso2;

  return out;

}

// [[Rcpp::export]]
DataFrame confintCpp(const DataFrame df, const int n, const int nthread) {

  int i;
  double *tauo = new double[n];
  double *taumo1 = new double[n];
  double *taumo2 = new double[n];
  double *tauso1 = new double[n];
  double *tauso2 = new double[n];
  NumericVector y = df[0];
  NumericVector c = df[1];
  NumericVector x = df[2];
  NumericVector p = df[3];

  #ifdef _OPENMP
  // Rprintf("nthread %d\n", nthread);
  // Rprintf("maxthread %d\n", omp_get_max_threads());
  // omp_set_num_threads(std::min(nthread, omp_get_max_threads()));
  omp_set_num_threads(nthread);
  #endif
  #pragma omp parallel default(shared), private(i)
  {
    #pragma omp for nowait
    for (i = 0; i < n; i++) {
      jsso(y, c, x, p, i, &tauo[i], &taumo1[i], &taumo2[i], &tauso1[i], &tauso2[i]);
    }
  }

  DataFrame out = DataFrame::create(
    Named("jktauo") = array2nvec(tauo, n),
    Named("jktaumo1") = array2nvec(taumo1, n),
    Named("jktaumo2") = array2nvec(taumo2, n),
    Named("jktauso1") = array2nvec(tauso1, n),
    Named("jktauso2") = array2nvec(tauso2, n)
  );

  delete [] tauo;
  delete [] taumo1;
  delete [] taumo2;
  delete [] tauso1;
  delete [] tauso2;

  return out;

}

void sso(const NumericVector y, const NumericVector c,
         const NumericVector x, const NumericVector p,
         double *tauo, double *taumo1,
         double *taumo2, double *tauso1, double *tauso2) {

  int l, a1, a2, b1, b2, denom0;
  int n = y.size();
  double denom0c, p2;
  int tx = 0;
  int ty = 0;
  double tau1 = 0;
  double mtau1 = 0;
  double mtau2 = 0;
  double denom1 = 0;

  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      l = (std::min(y[i], y[j]) < std::min(c[i], c[j]));
      a1 = 2*(y[i] - y[j] > 0) - 1;
      b1 = 2*(x[i] - x[j] > 0) - 1;
      a2 = a1*(y[i] != y[j]);
      b2 = b1*(x[i] != x[j]);
      tx += (y[i] == y[j]);
      ty += (x[i] == x[j]);
      if (i < j) {
        tau1 += l*a1*b1;
        p2 = pow(std::max(p[i], p[j]), 2);
        if (p2 > 0) {
          mtau1 += l*a1*b1/p2;
          mtau2 += l*a2*b2/p2;
          denom1 += abs(l*a1*b1/p2);
        }
      }
    }
  }
  denom0 = R::choose(n, 2);
  denom0c = sqrt(n*(n - 1)/2 - tx)*sqrt(n*(n - 1)/2 - ty);

  *tauo = tau1/denom0;
  *taumo1 = mtau1/denom0;
  *taumo2 = mtau1/denom1;
  *tauso1 = mtau1/denom0c;
  *tauso2 = mtau2/denom0c;

}

void jsso(const NumericVector y, const NumericVector c,
         const NumericVector x, const NumericVector p,
         const int r, double *tauo, double *taumo1,
         double *taumo2, double *tauso1, double *tauso2) {

  int l, a1, a2, b1, b2, denom0;
  int n = y.size();
  double denom0c, p2;
  int tx = 0;
  int ty = 0;
  double tau1 = 0;
  double mtau1 = 0;
  double mtau2 = 0;
  double denom1 = 0;

  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      if (i != r && j != r) {
        l = (std::min(y[i], y[j]) < std::min(c[i], c[j]));
        a1 = 2*(y[i] - y[j] > 0) - 1;
        b1 = 2*(x[i] - x[j] > 0) - 1;
        a2 = a1*(y[i] != y[j]);
        b2 = b1*(x[i] != x[j]);
        tx += (y[i] == y[j]);
        ty += (x[i] == x[j]);
        if (i < j) {
          tau1 += l*a1*b1;
          p2 = pow(std::max(p[i], p[j]), 2);
          if (p2 > 0) {
            mtau1 += l*a1*b1/p2;
            mtau2 += l*a2*b2/p2;
            denom1 += abs(l*a1*b1/p2);
          }
        }
      }
    }
  }
  n = n - 1;
  denom0 = R::choose(n, 2);
  denom0c = sqrt(n*(n - 1)/2 - tx)*sqrt(n*(n - 1)/2 - ty);

  *tauo = tau1/denom0;
  *taumo1 = mtau1/denom0;
  *taumo2 = mtau1/denom1;
  *tauso1 = mtau1/denom0c;
  *tauso2 = mtau2/denom0c;

}

NumericVector array2nvec(const double *ary, const int n) {

  NumericVector nvec(n);
  for (int i = 0; i < n; i++) {
    nvec[i] = ary[i];
  }
  return nvec;

}
