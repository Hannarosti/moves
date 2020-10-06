#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

List fillmat(int n, NumericMatrix Ht)
{
  
  NumericVector ii(n*n*n);  //   i <- numeric(n^2*n)
  NumericVector jj(n*n*n);  //   j <- numeric(n^2*n)
  NumericVector zz(n*n*n);  //   zz <- numeric(n^2*n)
  
  int k;
  int h;
  int z;
  int hh=0;
  double valor;
  
  for(k=0; k<n; k++){
    for(h=0; h<n; h++){
      for(z=0; z<n; z++){
        if(h==k || z==k){
          if(h==z){
            valor=1;
          }
          else valor=0;
        }
        else valor= Ht[h*n+z];
        jj[hh]  = k*n+h +1;
        ii[hh]  = k*n+z +1;
        zz[hh]  = valor; //T taux(ii,jj,valor);
        hh++;
      }
    }
  }
  return DataFrame::create(Named("ii") = ii, Named("jj")=jj, Named("zz")=zz);
}