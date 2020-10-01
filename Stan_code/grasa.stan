
functions{
  //Function that returns the log pdf of the wrapped-Cauchy
  real wrappedCauchy_lpdf(real y, real mu, real rho) {
    return(- log(2*pi()) + log((1-rho^2)/(1+rho^2-2*rho*cos(y-mu))));
  }
}

data {
    int<lower=0> T; // length of the time series
    vector[T] steps; // step lengths
    vector[T] turns; // turning turns
    vector[T] m; // direction to patch
    vector[T] phi; // movement direction
    int<lower=1> N; // number of states
    vector[T] idx; // in patch indicator
}

parameters {
  //vector<lower=0, upper=1>[N] rho;
  //vector<lower=0>[N] scale;
  //positive_ordered[N] scale;
  //vector<lower=0>[N] shape;
  real<lower=0> scale;
  real<lower=0> shape;
  real<lower=0, upper=1> rho;
  vector<lower=0>[N] theta;
  real<lower=1> gmax;
  real g0;
}  

transformed parameters {
  vector[N] log_theta[T];
  vector[T] g;
  
  g[1] = g0;

  for(i in 2:T){
    g[i] = fmin(g[i-1] + theta[1]*idx[i-1] - theta[2]*(1-idx[i-1]), gmax);
  }  
  
  for(t in 1:T){
    log_theta[t,1] = log(Phi_approx(g[t]));
    log_theta[t,2] = log(1-Phi_approx(g[t]));
  }
}

model {
  g0 ~ normal(0.5,0.1);
  rho ~ beta(2, 2);
  shape ~ gamma(12, 6);
  scale ~ normal(0,5);
  gmax ~ normal(5,0.1);
  theta ~ normal(0,1);
    // likelihood computation
    for (t in 1:T) {
      if(steps[t]>=0 && turns[t]>-pi() ){
        vector[N] lps = log_theta[t];
        //for(n in 1:N){
          lps[1] += wrappedCauchy_lpdf(turns[t] | 0, rho);
          lps[2] += wrappedCauchy_lpdf(phi[t] | m[t], rho);
          target += log_sum_exp(lps);
          target += weibull_lpdf(steps[t] | shape, scale);
        //}
    }
  }
}

generated quantities {
  matrix[T,N] stateProbs;
  real log_p[T,N];
  for(t in 1:T){
    stateProbs[t,1] = negative_infinity();
    stateProbs[t,2] = negative_infinity();
    log_p[t,1] = 0;
    log_p[t,2] = 0;
    if(steps[t]>=0 && turns[t]>-pi() ){
      log_p[t,1] = log_theta[t,1] +
       wrappedCauchy_lpdf(turns[t] | 0, rho);
       log_p[t,2] = log_theta[t,1] + 
       wrappedCauchy_lpdf(phi[t] | m[t], rho);
    }
    for(i in 1:N){
      if(steps[t]>=0 && turns[t]>-pi()){
        stateProbs[t,i] = exp(log_p[t,i])/(exp(log_sum_exp(log_p[t,])));
      }
    }
  }
}
