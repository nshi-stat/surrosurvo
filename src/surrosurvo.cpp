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

  const NumericVector y = df[0];
  const NumericVector c = df[1];
  const NumericVector x = df[2];
  const NumericVector d = df[3];
  const NumericVector p = df[4];
  const NumericVector q = df[5];

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

  DataFrame out = DataFrame::create(
    Named("tauo") = 2*tau1/n1,
    Named("taumo1") = 2*mtau1/n1,
    Named("taumo2") = mtau1/denom1,
    Named("tauso") = mtau2/sqrt(n1*0.5 - tx)/sqrt(n1*0.5 - ty)
  );

  return out;

}

// [[Rcpp::export]]
DataFrame confintCpp(const DataFrame df, const int n, const int censtypen,
                     const int nthread) {

  int i;
  NumericMatrix tau(n, 4);
  const NumericVector y = df[0];
  const NumericVector c = df[1];
  const NumericVector x = df[2];
  const NumericVector d = df[3];
  const NumericVector p = df[4];
  const NumericVector q = df[5];

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
    tau(i, _) = jsso(y, c, x, d, p, q, n, censtypen, i);
  }
}

DataFrame out = DataFrame::create(
  Named("jktauo") = tau(_, 0),
  Named("jktaumo1") = tau(_, 1),
  Named("jktaumo2") = tau(_, 2),
  Named("jktauso") = tau(_, 3)
);

return out;

}

NumericVector jsso(const NumericVector y, const NumericVector c,
                   const NumericVector x, const NumericVector d,
                   const NumericVector p, const NumericVector q,
                   const int n, const int censtypen, const int r) {

  int l, a1, a2, b1, b2;
  double p2;
  int tx = 0;
  int ty = 0;
  double tau1 = 0;
  double mtau1 = 0;
  double mtau2 = 0;
  double denom1 = 0;
  double n1 = (n - 1)*(n - 2);
  NumericVector out(4);

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

  out(0) = 2*tau1/n1;
  out(1) = 2*mtau1/n1;
  out(2) = mtau1/denom1;
  out(3) = mtau2/sqrt(n1*0.5 - tx)/sqrt(n1*0.5 - ty);

  return out;

}
