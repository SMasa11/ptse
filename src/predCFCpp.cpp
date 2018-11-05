#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]

//Prototypes
NumericVector seqCppFloat(double &first, double &last, int &length);

DataFrame appx(NumericVector x, NumericVector y, NumericVector xout);

//Function calll
// [[Rcpp::export]]
NumericVector predCFCpp(NumericMatrix &FY1XWpred,NumericMatrix &QY1XCF,NumericVector &itrY, NumericVector &itrU, NumericVector &weights, int &nGridsFine, int n, int weighted, int discY)
{
  //double first {itrY.begin()};
  //double last {itrY.end()};
  NumericVector itrI = seqCppFloat(itrY[0],itrY[itrY.length()-1],nGridsFine);
  NumericMatrix FY0XCF(n,itrI.length());

  NumericVector FY1XWpredInt(FY1XWpred.ncol());
  for (int i = 0; i < n; i++)
  {
    DataFrame dfFY1XWpredInt = appx(itrY,FY1XWpred(i,_),itrI);
    FY1XWpredInt = dfFY1XWpredInt["y"];
    FY1XWpredInt.sort();

    int curl = 0;
    // for every percentile j in [0,1]
    for (int j = 0; j < itrI.length(); j++)
    {
      bool flagPassProcess = 0;

      // retrieve the case of -1 for LB
      if (discY == 1)
      {
        if (QY1XCF(i,j) == -1)
        {
          flagPassProcess = 1;
          FY0XCF(i,j) = 0;
        }
        //if (QY1XCF(i,j) == 0)
        //{
        //  flagPassProcess = 1;
        //  FY0XCF(i,j) = FY1XWpredInt[0];
        //}
      }

      if (flagPassProcess == 0)
      {
        bool flagEnd = 0;
        // search for the value of Y0, l
        for (int l = curl; l < itrI.length(); l++)
        {
          // such that QY1(j) = l
          if (itrI[l] >= QY1XCF(i,j))
          {
            // then that is the value of FY0 to match FY1
            FY0XCF(i,j) = FY1XWpredInt[l];
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
          FY0XCF(i,j) = 1;
        }
      }
    }
    //rearrangement
    NumericVector FY0XCFI(FY0XCF.ncol());
    FY0XCFI = FY0XCF(i,_);
    FY0XCFI.sort();
    FY0XCF(i,_) = FY0XCFI;
  }

  NumericVector FY0XCFM(itrI.length());
  for (int j = 0; j < itrI.length(); j++)
  {
    double meanVal = 0;
    for (int i = 0; i < n; i++)
    {
      if (weighted == 1)
      {
        meanVal += FY0XCF(i,j)*weights[i]/n;
      } else {
        meanVal += FY0XCF(i,j)/n;
      }
    }
    FY0XCFM[j] = meanVal;
  }
  FY0XCFM.sort();
  return FY0XCFM;
}

//Inner function definitions
DataFrame appx(NumericVector x, NumericVector y, NumericVector xout)
{
  Function f("approx");

  return f(Named("x")=x,Named("y")=y,Named("xout")=xout,Named("method")="linear",Named("ties")="ordered");
}

NumericVector seqCppFloat(double &first, double &last, int &length) {
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
nGridsFine <- 1000
FY1Wpred <- matrix(n*10,n,10)
QY1XCF <- matrix(n*nGridsFine,n,nGridsFine)
for (i in 1:n)
{
  FY1Wpred[i,] <- sort(runif(10,0,1))
  QY1XCF[i,] <- sort(rnorm(nGridsFine,0,1))
}
weights <- rexp(n)

FY0XCF <- predCFCpp(FY1Wpred,QY1XCF,itrY,itrU,weights,nGridsFine,n,0)
*/
