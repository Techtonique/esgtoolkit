#include <Rcpp.h>
using namespace Rcpp;


double crossprod_cpp(NumericVector x, NumericVector y)
{
  unsigned long int n = x.size();
  if (y.size() != n) {
    ::Rf_error("you must have x.size() == y.size()");
  }
  double res = 0;
  
  for(int i = 0; i < n; i++) {
    res += x(i)*y(i);
  }
  return(res);
}

double prod(NumericVector x)
{
  const unsigned long int n = x.size();
  double res = 0;
  for (unsigned long int i = 0; i < n; i++)
  {
    res *= x(i);
  }
  return(res);
}

NumericVector squareVector(NumericVector x)
{
  const unsigned long int n = x.size();
  NumericVector res(n);
  for (unsigned long int i = 0; i < n; i++)
  {
    res(i) = std::pow(x(i), 2);
  }
  return(res);
}

NumericVector vectorSlice(NumericVector x, const unsigned long int start, 
                          const unsigned long int end)
{
  if (start > end) {
    ::Rf_error("start index must be <= end index");
  }
  return(x[Rcpp::Range(start, end)]);
}

NumericVector whichSup(NumericVector zi, NumericVector zj)
{
  const unsigned long int n = zi.size();
  NumericVector res(n);
  
  if (zj.size() != n) {
    ::Rf_error("inputs must have the same size");
  }
  
  for (unsigned long int i = 0; i < n; i++)
  {
    if (zi(i) < zj(i))
    {
      res(i) = 1;   
    } else {
      res(i) = 0;   
    }
  }
  return(res);
}

// adapted from vrtest::Mammen
// [[Rcpp::export]]
NumericVector Mammen_cpp(const unsigned long int n)
{
  double p = 0;
  double sqrt5 = std::sqrt(5);
  NumericVector zmat(n);
  NumericVector u(n);
  
  p = (sqrt5 + 1)/(2 * sqrt5);
  zmat = rep(-(sqrt5 - 1)/2, n);
  u = runif(n, 0, 1);
  for (unsigned long int i = 0; i < n; i++)
  {
    if (u(i) > p)
    {
      zmat[i] = 0.5*(sqrt5 + 1);
    }
  }
  return(zmat);
}

// https://stackoverflow.com/questions/47246200/how-to-slice-rcpp-numericvector-for-elements-2-to-101
// [[Rcpp::export]]
List DLtest_cpp(NumericVector y, const unsigned long int p)
{
  List res = List::create(Named("Cpstat") =  0, _["Kpstat"] = 0);
  const unsigned long int n = y.size();
  // Rprintf("n %i \n", n);
  // Rprintf("the value of y[%i] : %f \n", 0, y[0]);
  // Rprintf("the value of y[%i] : %f \n", 1, y[1]);
  // Rprintf("the value of y[%i] : %f \n", 2, y[2]);
  NumericVector ym(n);
  double s2 = 0;
  double sum1, sum2, Cp, Kp, zi, zj, tem1;
  NumericVector sum3(n - p);
  unsigned short int indicate = 0;
  
  ym = y - mean(y);
  sum2 = 0;
  s2 = crossprod_cpp(ym, ym)/(n - p);
  // Rprintf("the value of ym[%i] : %f \n", 0, ym[0]);
  // Rprintf("the value of ym[%i] : %f \n", 1, ym[1]);
  // Rprintf("the value of ym[%i] : %f \n", 2, ym[2]);
  //Rprintf("s2 %f \n", s2);
  
  for(unsigned long int j = p; j < n; j++) {
    sum1 = 0;
    for(unsigned long int i = p; i < n; i++) {
      indicate = 0;
      zi = vectorSlice(ym, i - 1, i - p)(0); //ym((i-1):(i-p));
      zj = vectorSlice(ym, j - 1, j - p)(0); //ym((j-1):(j-p));
      tem1 = (zi <= zj)*1;
      if(tem1 == 1){
        indicate = 1;
      }
      sum1 += ym(i)*indicate;
    }
    sum2 += std::pow(sum1, 2);
    sum3(j - p) = std::abs(sum1/std::sqrt(n - p));
  }
  
  // Rprintf("sum2 %f \n", sum2);
  // Rprintf("s2 %f \n", s2);
  
  Cp = sum2/(s2*std::pow((n-p), 2));
  Kp = max(sum3)/sqrt(s2);
  res["Cpstat"] = Cp;
  res["Kpstat"] = Kp;
  
  return(res);
}

/*** R
#set.seed(123); Mammen_cpp(10)
*/