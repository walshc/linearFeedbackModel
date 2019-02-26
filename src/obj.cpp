#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector quasiDifference(NumericVector theta, IntegerMatrix idx, int nT,
    NumericMatrix data)
{
  int K = data.ncol();
  int N = data.nrow();

  int n = N / (nT - 1);

  int i;

  /* Find mu for every observation: */
  NumericVector xBeta(N);
  NumericVector mu(N);
  if (K > 2) {
    for (i = 0; i < N; i++) {
      xBeta[i] = 0.0;
        for (int k = 2; k < K; k++) {
          xBeta[i] += theta[k - 1] * data(i, k);
        }
      mu[i] = exp(xBeta[i]);
    }
  } else {
    for (i = 0; i < N; i++) {
      mu[i] = 1.0;
    }
  }

  /* Get the lag of mu and the 2nd lag of the dependent variable: */
  IntegerVector id  = idx(_, 0);
  IntegerVector t   = idx(_, 1);
  NumericVector y   = data(_, 0);
  NumericVector l_y = data(_, 1);
  NumericVector l2_y(N);
  NumericVector l_mu(N);
  for (i = 0; i < N; i++) {
    if (t[i] == 2) {
      l2_y[i] = 0.0;
      l_mu[i] = 0.0;
    } else {
      l2_y[i] = l_y[i - 1];
      l_mu[i] = mu[i - 1];
    }
  }

  /* Apply the Wooldridge quasi-differencing transformation to get q: */
  NumericVector q(n * (nT - 2));
  int j = 0;
  for (i = 0; i < N; i++) {
    if (t[i] > 2) {
      q[j] = (y[i] - theta[0] * l_y[i]) * l_mu[i] / mu[i] - \
          (l_y[i] - theta[0] * l2_y[i]);
      j += 1;
    }
  }
  return q;
}

// [[Rcpp::export]]
double GMM(NumericVector theta, IntegerMatrix idx, int nT,
    NumericMatrix data, NumericMatrix Z, NumericMatrix W)
{
  /* Perform the Wooldridge quasi differencing: */
  NumericVector q = quasiDifference(theta, idx, nT, data);

  int L = Z.ncol();  // Number of instruments
  int N = q.size();  // Number of observations (N * (T - 2))*/

  int i, j, l;

  /* Calculate each of the moments: */
  NumericVector g(L);
  for (l = 0; l < L; l++) {
    for (i = 0; i < N; i++) {
        g[l] += q[i] * Z(i, l);
    }
  }

  /* Get the objective using the moments and weighting matrix: */
  double obj;
  for (i = 0; i < L; i++) {
    for (j = 0; j < L; j++) {
      obj += g[i] * W(i, j) * g[j];
    }
  }
  return obj;
}
