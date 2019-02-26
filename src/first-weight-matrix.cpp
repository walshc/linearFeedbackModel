#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix firstWeightMatrix(NumericMatrix Z)
{

  int L = Z.ncol();  // Number of instruments
  int N = Z.nrow();  // Number of observations (N * (T - 2))*/

  int i, j, k;

  NumericMatrix W(L, L);

  for (j = 0; j < L; j++) {
    for (k = 0; k < L; k++) {
      for (i = 0; i < N; i++) {
          W(j, k) += Z(i, j) * Z(i, k) / N;
      }
    }
  }

  return W;
}
