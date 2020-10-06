#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

NumericVector fillvec(int n, NumericMatrix H)
{
  int m = n*n;
  NumericVector ah(m);
  int i,j;
  
  for(i=0; i<n; i++){
    for(j=0; j<n; j++){
      if(i==j) 
        ah[i*n+j] = 0;
      else
      	ah[i*n+j] = H[j*n+i];
     }  
  }
        
	return ah;
}