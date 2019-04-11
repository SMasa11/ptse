#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]

//Prototypes
NumericVector seqCppFloatQ1(double &first, double &last, int &length);

DataFrame appxQ1(NumericVector x, NumericVector y, NumericVector xout);

//Function calll
// [[Rcpp::export]]
NumericVector predQ1Cpp(NumericMatrix &FY1Wpred,NumericMatrix &FY0Wpred,NumericVector &itrY, NumericVector &itrU, int &nGridsFine, int n, int discY, int UB)
{
  //double first {itrY.begin()};
  //double last {itrY.end()};
  NumericVector itrI = seqCppFloatQ1(itrY[0],itrY[itrY.length()-1],nGridsFine);
  NumericMatrix QY1XCF(n,itrI.length());

  NumericVector FY1WpredInt(FY1Wpred.ncol());
  NumericVector FY0WpredInt(FY0Wpred.ncol());
  for (int i = 0; i < n; i++)
  {
    DataFrame dfFY1WpredInt = appxQ1(itrY,FY1Wpred(i,_),itrI);
    FY1WpredInt = dfFY1WpredInt["y"];
    //reordering
    FY1WpredInt.sort();
    DataFrame dfFY0WpredInt = appxQ1(itrY,FY0Wpred(i,_),itrI);
    FY0WpredInt = dfFY0WpredInt["y"];
    FY0WpredInt.sort();

    int curl = 0;
    // for each j in grids of Y0 to evaluate
    for (int j = 0; j < itrI.length(); j++)
    {
      bool flagEnd = 0;
      bool flagPassProcess = 0;
      // if discY is on, then additional check is required
      if (discY == 1)
      {
        // we need to check if FY1(0) <= FY0(j)
        if (FY1WpredInt[0] > FY0WpredInt[j])
        {
          // in this case, bound must be applied
          flagPassProcess = 1;
          if (UB == 1)
          {
            // upper bound is FY1(0), i.e., QY1 = 0
            QY1XCF(i,j) = 0;
          } else {
            // lower bound is 0-. This must be retrieved later as -1 is not just latent var.
            QY1XCF(i,j) = -1;
          }
        }
      }

      if (flagPassProcess == 0)
      {
        // find the grids of Y1, l such that FY1(l) = FY0(j)
        for (int l = curl; l < itrI.length(); l++)
        {
          if (FY1WpredInt[l] >= FY0WpredInt[j])
          {
            QY1XCF(i,j) = itrI[l];
            curl = l;
            break;
          }
          if (l == itrI.length()-1)
          {
            flagEnd = 1;
          }
        }
        if (flagEnd == 1)
        {
          QY1XCF(i,j) = itrI[itrI.length()-1];
        }
      }
    }
  }

  return QY1XCF;
}

//Inner function definitions
DataFrame appxQ1(NumericVector x, NumericVector y, NumericVector xout)
{
  Function f("approx");

  return f(Named("x")=x,Named("y")=y,Named("xout")=xout,Named("method")="linear",Named("ties")="ordered");
}

NumericVector seqCppFloatQ1(double &first, double &last, int &length) {
  NumericVector y(length);
  double dif = (last - first) / (length - 1);

  double iFl = 0.0;
  for (int i = 0; i < length; i++)
  {
    y[i] = first + iFl*dif;
    iFl++;
  }
  return y;
}

/*** R
itrY <- seq(0.1,11.2,length=10)
itrU <- seq(0.1,0.9,length=10)
n <- 2000
FY1Wpred <- matrix(n*10,n,10)
FY0Wpred <- matrix(n*10,n,10)
for (i in 1:n)
{
  FY1Wpred[i,] <- sort(runif(10,0,1))
  FY0Wpred[i,] <- sort(runif(10,0,1))
}

nGridsFine <- 1000

QY1XCF <- predQ1Cpp(FY1Wpred,FY0Wpred,itrY,itrU,nGridsFine,n)
*/
