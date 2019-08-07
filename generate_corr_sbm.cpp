#include <Rcpp.h>
using namespace Rcpp;

#include <time.h>
#include <stdlib.h>

// Initialization, should only be called once.

// [[Rcpp::export]]
NumericMatrix simulate_corr_SBM(double a,double b,NumericMatrix Aer,double rho,double p) {
  srand(time(NULL)); 
  int n = Aer.nrow();
  NumericMatrix A(n,n);
  int br = 0;
  
  for (int i = 0; i < n; i++) {
    for (int j = i; j <= n ; j++) {
      double r = ((double) rand() / (RAND_MAX));
      
      if (i < n/2 && j < n/2) { 
        double param1 = a +(rho/p)* pow((p*a*(1-p)*(1-a)),.5);
        double param2 = (a - param1*p)/(1-p);
        if (Aer(i,j) == 1) {
          if (r >= 1 - param1) {
            br = 1;
          }
        } else {
          if (r >= 1 - param2) {
            br = 1;
          }
        }
      } else { 
        double param1 = b +(rho/p)*pow((p*a*(1-p)*(1-a)),.5);
        double param2 = (b - param1*p)/(1-p);
        if (Aer(i,j) == 1) {
          if (r >= 1 - param1) {
            br = 1;
          }
        } else {
          if (r >= 1 - param2) {
            br = 1;
          }
        }
      }
      A(i,j) = br;
      A(j,i) = A(i,j);
      br = 0;
    }
  }

  return A;
  
}



    
