#include <Rcpp.h>
#include <Rmath.h>
#include "surrosurvo.h"

#ifdef _OPENMP
#include <omp.h>
// [[Rcpp::plugins("openmp")]]
#endif

using namespace Rcpp;

// [[Rcpp::export]]
DataFrame surrosurvoCpp(const DataFrame df, const int n, const int censtypen) {

  double *tauo = new double[1];
  double *taumo1 = new double[1];
  double *taumo2 = new double[1];
  double *tauso = new double[1];

  NumericVector y = df[0];
  NumericVector c = df[1];
  NumericVector x = df[2];
  NumericVector d = df[3];
  NumericVector p = df[4];
  NumericVector q = df[5];
  sso(y, c, x, d, p, q, n, censtypen,
      &tauo[0], &taumo1[0], &taumo2[0], &tauso[0]);

  DataFrame out = DataFrame::create(
    Named("tauo") = array2nvec(tauo, 1),
    Named("taumo1") = array2nvec(taumo1, 1),
    Named("taumo2") = array2nvec(taumo2, 1),
    Named("tauso") = array2nvec(tauso, 1)
  );

  delete [] tauo;
  delete [] taumo1;
  delete [] taumo2;
  delete [] tauso;

  return out;

}

// [[Rcpp::export]]
DataFrame confintCpp(const DataFrame df, const int n, const int censtypen,
                     const int nthread) {

  int i;
  double *tauo = new double[n];
  double *taumo1 = new double[n];
  double *taumo2 = new double[n];
  double *tauso = new double[n];
  NumericVector y = df[0];
  NumericVector c = df[1];
  NumericVector x = df[2];
  NumericVector d = df[3];
  NumericVector p = df[4];
  NumericVector q = df[5];

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
    jsso(y, c, x, d, p, q, n, censtypen, i,
         &tauo[i], &taumo1[i], &taumo2[i], &tauso[i]);
  }
}

DataFrame out = DataFrame::create(
  Named("jktauo") = array2nvec(tauo, n),
  Named("jktaumo1") = array2nvec(taumo1, n),
  Named("jktaumo2") = array2nvec(taumo2, n),
  Named("jktauso") = array2nvec(tauso, n)
);

delete [] tauo;
delete [] taumo1;
delete [] taumo2;
delete [] tauso;

return out;

}

void sso(const NumericVector y, const NumericVector c,
         const NumericVector x, const NumericVector d,
         const NumericVector p, const NumericVector q,
         const int n, const int censtypen,
         double *tauo, double *taumo1, double *taumo2, double *tauso) {

  int l, a1, a2, b1, b2;
  double p2;
  int tx = 0;
  int ty = 0;
  double tau1 = 0;
  double mtau1 = 0;
  double mtau2 = 0;
  double denom1 = 0;
  double n1 = n*(n - 1);

  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      l = (std::min(y[i], y[j]) < std::min(c[i], c[j]))*
        (std::min(x[i], x[j]) < std::min(d[i], d[j]));
      a1 = 2*(y[i] - y[j] > 0) - 1;
      b1 = 2*(x[i] - x[j] > 0) - 1;
      a2 = a1*(y[i] != y[j]);
      b2 = b1*(x[i] != x[j]);
      ty += (y[i] == y[j]);
      tx += (x[i] == x[j]);
      tau1 += l*a1*b1;
      if (censtypen == 1) {
        p2 = pow(std::min(std::max(p[i], p[j]), std::max(q[i], q[j])), 2);
      } else {
        p2 = pow(std::max(p[i], p[j])*std::max(q[i], q[j]), 2);
      }
      if (p2 > 0) {
        mtau1 += l*a1*b1/p2;
        mtau2 += l*a2*b2/p2;
        denom1 += abs(l*a1*b1/p2);
      }
    }
  }

  *tauo = 2*tau1/n1;
  *taumo1 = 2*mtau1/n1;
  *taumo2 = mtau1/denom1;
  *tauso = mtau2/sqrt(n1*0.5 - tx)/sqrt(n1*0.5 - ty);

}

void jsso(const NumericVector y, const NumericVector c,
          const NumericVector x, const NumericVector d,
          const NumericVector p, const NumericVector q,
          const int n, const int censtypen, const int r,
          double *tauo, double *taumo1, double *taumo2, double *tauso) {

  int l, a1, a2, b1, b2;
  double p2;
  int tx = 0;
  int ty = 0;
  double tau1 = 0;
  double mtau1 = 0;
  double mtau2 = 0;
  double denom1 = 0;
  double n1 = (n - 1)*(n - 2);

  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      if (i != r && j != r) {
        l = (std::min(y[i], y[j]) < std::min(c[i], c[j]))*
          (std::min(x[i], x[j]) < std::min(d[i], d[j]));
        a1 = 2*(y[i] - y[j] > 0) - 1;
        b1 = 2*(x[i] - x[j] > 0) - 1;
        a2 = a1*(y[i] != y[j]);
        b2 = b1*(x[i] != x[j]);
        tx += (y[i] == y[j]);
        ty += (x[i] == x[j]);
        tau1 += l*a1*b1;
        if (censtypen == 1) {
          p2 = pow(std::min(std::max(p[i], p[j]), std::max(q[i], q[j])), 2);
        } else {
          p2 = pow(std::max(p[i], p[j])*std::max(q[i], q[j]), 2);
        }
        if (p2 > 0) {
          mtau1 += l*a1*b1/p2;
          mtau2 += l*a2*b2/p2;
          denom1 += abs(l*a1*b1/p2);
        }
      }
    }
  }

  *tauo = 2*tau1/n1;
  *taumo1 = 2*mtau1/n1;
  *taumo2 = mtau1/denom1;
  *tauso = mtau2/sqrt(n1*0.5 - tx)/sqrt(n1*0.5 - ty);

}

NumericVector array2nvec(const double *ary, const int n) {

  NumericVector nvec(n);
  for (int i = 0; i < n; i++) {
    nvec[i] = ary[i];
  }
  return nvec;

}
