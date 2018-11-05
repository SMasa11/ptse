#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
NumericVector getQuantileCpp(NumericVector &FY, NumericVector &kI, NumericVector &kU) {
  NumericVector QY(kU.length());
  FY.sort();

  int curj = 1;
  for (int i = 0; i < kU.length(); i++)
  {
    for (int j = curj; j < kI.length(); j++)
    {
      if (FY[j] >= kU[i])
      {
        QY[i] = kI[j];
        curj = j;
        break;
      }
    }
  }
  return QY;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R

*/
