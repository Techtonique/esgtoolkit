#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix rnormESGcpp(const int N, const int M) 
{
  // N : Number of independent simulations
  // M : Number of dates for projection
  NumericMatrix mat(M, N);
  
  RNGScope scope;  
   for(int j = 0; j < N; j++) {
    mat(_, j) = rnorm(M) ;
  }
  
  return mat;
}


// [[Rcpp::export]]
NumericMatrix rOUESGcpp(const int N, const int horizon, const double Delta, const double x0, 
                        NumericVector theta, NumericMatrix eps) 
{
  // N : Number of independent simulations

    double theta1 = theta(0);
    double theta2 = theta(1);
    double theta3 = theta(2);
    int nbdates = horizon/Delta + 1;
    NumericMatrix out(nbdates, N);
    int nbdates_1 = nbdates - 1;
    double sqrt_Delta = sqrt(Delta); 

    for(int j = 0; j < N; j++)
    {
      out(0, j) = x0;        
    }
        
    if (theta1 == 0) 
    { 
      if (theta2 == 0)
      {
        for(int i = 0; i < nbdates_1; i++)
        {
          out((i+1), _) = out(i, _) + theta3*sqrt_Delta*eps(i, _);
        }
      }
      else
      {
        for(int i = 0; i < nbdates_1; i++)
        {
          out((i+1), _) = out(i, _)*(1 - theta2*Delta) + theta3*sqrt_Delta*eps(i, _);
        }
      }      
    } else {      
              for(int i = 0; i < nbdates_1; i++)
              {
                  out((i+1), _) = out(i, _) + (theta1 - theta2*out(i, _))*Delta + theta3*sqrt_Delta*eps(i, _); 
              }
            }

  return out;
}

// [[Rcpp::export]]
NumericMatrix rOUESGcppexact(const int N, const int horizon, const double Delta, const double x0, 
                        NumericVector theta, NumericMatrix eps) 
{
  // N : Number of independent simulations

    double theta1 = theta(0);
    double theta2 = theta(1);
    double theta3 = theta(2);
    
    
    int nbdates = horizon/Delta + 1;
    NumericMatrix out(nbdates, N);
    int nbdates_1 = nbdates - 1;
    double sqrt_Delta = sqrt(Delta); 

    for(int j = 0; j < N; j++)
    {
      out(0, j) = x0;        
    }
    
    /*
    exp_a <- exp(-a_opt*delta_t)
    sqrt_exp_2a <- sqrt((1-exp(-2*a_opt*delta_t))/(2*a_opt))
    x_up <- exp_a*x_down + sigma_opt*sqrt_exp_2a*Zmat1[,j] 
    */
    
    if (theta1 == 0) 
    { 
      if (theta2 == 0)
      {
        for(int i = 0; i < nbdates_1; i++)
        {
          out((i+1), _) = out(i, _) + theta3*sqrt_Delta*eps(i, _);
        }
      }
      else
      {
        double sqrt_exp_2a = sqrt((1-exp(-2*theta2*Delta))/(2*theta2));
        double exp_a = exp(-theta2*Delta);
        
        for(int i = 0; i < nbdates_1; i++)
        {
          out((i+1), _) = exp_a*out(i, _) + theta3*sqrt_exp_2a*eps(i, _);
        }
      }      
    } else {    
              double sqrt_exp_2a = sqrt((1-exp(-2*theta2*Delta))/(2*theta2));
              double exp_a = exp(-theta2*Delta);
              
              for(int i = 0; i < nbdates_1; i++)
              {
                  out((i+1), _) = exp_a*out(i, _) + (theta1/theta2)*(1 - exp_a) + theta3*sqrt_exp_2a*eps(i, _);
              }
            }

  return out;
}

// [[Rcpp::export]]
NumericMatrix rCIRESGcpp(const int N, const int horizon, const double Delta, const double x0, 
                         NumericVector theta, NumericMatrix eps) 
{
  // N : Number of independent simulations
    double theta1 = theta(0);
    double theta2 = theta(1);
    double theta3 = theta(2);
    double theta3squared = pow(theta3, 2);
    int nbdates = horizon/Delta + 1;
    NumericMatrix mat(nbdates, N);
    NumericMatrix out(nbdates, N);
    int nbdates_1 = nbdates - 1;
    double sqrt_Delta = sqrt(Delta); 
    
    for(int j = 0; j < N; j++)
    {
      out(0, j) = sqrt(x0);
      mat(0, j) = x0;
    }

      for(int i = 0; i < nbdates_1; i++)
      {
        out((i+1), _) = out(i, _) + 0.5*((theta1 - theta2*pow(out(i, _), 2)) - 0.25*theta3squared)*Delta/out(i, _) + 
        0.5*theta3*sqrt_Delta*eps(i, _);
        mat((i+1), _) = pow(out((i+1), _), 2);
      }
              
  return mat;
}

// [[Rcpp::export]]
NumericMatrix rCIRESGcppexact(const int N, const int horizon, const double Delta, const double x0, 
                         NumericVector theta, NumericMatrix eps) 
{
  // N : Number of independent simulations
    double theta1 = theta(0);
    double theta2 = theta(1);
    double theta3 = theta(2);
    double theta3squared = pow(theta3, 2);
    double d = 4*theta1/theta3squared;
    int nbdates = horizon/Delta + 1;
    NumericMatrix mat(nbdates, N);
    int nbdates_1 = nbdates - 1;
    NumericVector X(N);
    NumericVector lambda(N);
    int NPoisson;
    
    for(int j = 0; j < N; j++)
    {
      mat(0, j) = x0;
    }
    
    double exp_a = exp(-theta2*Delta);
    double c = theta3squared*(1 - exp_a)/(4*theta2);
    
    if (d > 1)
    {
      for(int i = 0; i < nbdates_1; i++)
      { 
        lambda = mat(i, _)*exp_a/c;
        X = rchisq(N, d - 1);
        mat((i+1), _) = c*(X + (eps(i, _) + sqrt(lambda))*(eps(i, _) + sqrt(lambda)));
      }
    }
     else
    {
      for(int i = 0; i < nbdates_1; i++)
      {
        lambda = mat(i, _)*exp_a/c;
         for(int j = 0; j < N; j++)
        { 
          NPoisson = rpois(1, 0.5*lambda(j))[1];
          mat((i+1), j) = c*rchisq(1, d + 2*NPoisson)[1];
        } 
      }
    }
              
  return mat;
}

// [[Rcpp::export]]
NumericMatrix rGBMESGcpp(const int N, const int horizon, const double Delta, const double x0, 
                        NumericMatrix theta1, NumericMatrix theta2, NumericMatrix eps) 
{
  // N : Number of independent simulations
    int nbdates = horizon/Delta + 1;
    NumericMatrix mat(nbdates, N);
    NumericMatrix out(nbdates, N);
    int nbdates_1 = nbdates - 1;
    double sqrt_Delta = sqrt(Delta); 
    
    for(int j = 0; j < N; j++)
    {
      out(0, j) = log(x0);        
      mat(0, j) = x0;
    }
        
      for(int i = 0; i < nbdates_1; i++)
      {
        out((i+1), _) = out(i, _) + (theta1(i, _) - 0.5*pow(theta2(i, _), 2))*Delta + theta2(i, _)*sqrt_Delta*eps(i, _);
        mat((i+1), _) = exp(out((i+1), _));
      }
              
  return mat;
}

// [[Rcpp::export]]
NumericMatrix rGBMjumpsnormESGcpp(const int N, const int horizon, const double Delta,
                             const double x0, NumericMatrix theta1, NumericMatrix theta2,
                             const double lambda, const double mu, const double sigma, 
                             NumericMatrix eps) 
{
    // N : Number of independent simulations
    int nbdates = horizon/Delta + 1;
    NumericMatrix mat(nbdates, N);
    NumericMatrix out(nbdates, N);
    NumericVector NJ(N); 
    NumericVector MJ(N); 
    int nbdates_1 = nbdates - 1;
    double sqrt_Delta = sqrt(Delta); 
    double lambdadelta = lambda*Delta;
    double mulogJ = log(1+mu) - 0.5*pow(sigma, 2);
    
   for(int j = 0; j < N; j++)
    {
      out(0, j) = log(x0);        
      mat(0, j) = x0;
    }

      RNGScope scope;  
      for(int i = 0; i < nbdates_1; i++)
      {
        NJ = rpois(N, lambdadelta);            
        MJ = mulogJ*NJ + sigma*sqrt(NJ)*rnorm(N, 0, 1);              
        out((i+1), _) = out(i, _) + (theta1(i, _) - 0.5*pow(theta2(i, _), 2))*Delta + 
        theta2(i, _)*sqrt_Delta*eps(i, _) + MJ; 
        mat((i+1), _) = exp(out((i+1), _));
      }             
  return mat;
}

// [[Rcpp::export]]
NumericMatrix rGBMjumpskouESGcpp(const int N, const int horizon, const double Delta,
                             const double x0, NumericMatrix theta1, NumericMatrix theta2,
                             const double lambda, const double eta_up, const double eta_down,
                             const double p, NumericMatrix eps) 
{
    // N : Number of independent simulations
    int nbdates = horizon/Delta + 1;
    NumericMatrix mat(nbdates, N);
    NumericMatrix out(nbdates, N);
    NumericVector NJ(N); 
    double NJtemp;
    double K;
    double MJ; 
    int nbdates_1 = nbdates - 1;
    double sqrt_Delta = sqrt(Delta); 
    double lambdadelta = lambda*Delta;
    
    for(int j = 0; j < N; j++)
    {
      out(0, j) = log(x0);        
      mat(0, j) = x0;
    }
      
      RNGScope scope;  
      for(int i = 0; i < nbdates_1; i++)
      {
        NJ = rpois(N, lambdadelta);
        for(int j = 0; j < N; j++)
        {
          NJtemp = NJ(j);
            if (NJtemp == 0)
            {
              MJ = 0;
            }
            else
            {
              K =  rbinom(1, NJtemp, p)[1];
              MJ = rgamma(1, K, 1/eta_up)[1] - rgamma(1, NJtemp - K, 1/eta_down)[1];  
            }
          out((i+1), j) = out(i, j) + (theta1(i, j) - 0.5*pow(theta2(i, j), 2))*Delta + 
          theta2(i, j)*sqrt_Delta*eps(i, j) + MJ;
        }       
        mat((i+1), _) = exp(out((i+1), _));
      }              
  return mat;
}

// [[Rcpp::export]]
NumericVector TAGcorecpp(NumericVector sim, NumericVector sj_down, 
              NumericVector sj_up, const int n, const int p)
{
int i, j, divsomme;
double sim_j, somme;
NumericVector out(p);
  
  for(i = 0; i < p; i++)
  {
    somme = 0;
    divsomme = 0;  
    for(j = 0; j < n; j++)
    {
      sim_j = sim(j);    
      if((sim_j >= sj_down(i)) & (sim_j < sj_up(i))) 
      {
        somme = somme + sim_j;
        divsomme = divsomme + 1;
      } 
    }
    out(i) = somme/divsomme;
  }

return out;
}
