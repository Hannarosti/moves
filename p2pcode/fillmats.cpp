#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

List fillmats(int n, NumericMatrix Ht, NumericMatrix UH)
{

NumericVector ii(n*n*n);
NumericVector jj(n*n*n);
NumericVector zz(n*n*n);
NumericVector zu(n*n*n);

int k;
int h;
int z;
int hh=0;
double valor;
double valoru;

  for(k=0; k<n; k++){
    for(h=0; h<n; h++){
      for(z=0; z<n; z++){
        if(h==k || z==k){
          if(h==z){
            valor=1;
            valoru=1;
          }
          else valor=0, valoru =0;
          }
          else valor= Ht[h*n+z];
          valoru = UH[h*n+z];
          jj[hh]  = k*n+h +1;
          ii[hh]  = k*n+z +1;
          zz[hh]  = valor;
          zu[hh] = valoru;
          hh++;
      }
    }
  }
	return DataFrame::create(Named("ii") = ii, Named("jj")=jj, Named("zz")=zz, Named("zu")=zu);
}
