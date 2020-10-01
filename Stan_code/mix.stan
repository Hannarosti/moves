// j.m.morales June 2018

functions{
  //Function that returns the log pdf of the wrapped-Cauchy
  real wrappedCauchy_lpdf(real y, real mu, real rho) {
    return(- log(2*pi()) + log((1-rho^2)/(1+rho^2-2*rho*cos(y-mu))));
  }
}

data {
    int<lower=0> T; // length of the time series
    //int ID[T]; // track identifier
    vector[T] steps;
    vector[T] turns;
    int<lower=1> N;  // number of states
    real lb; // lower bound for shape
}

parameters {
  vector<lower=-1, upper=4.5>[N] mu; // this allows mu to be centered around pi
  vector<lower=0, upper=1>[N] rho;
  positive_ordered[N] mu_step;
  vector<lower=lb>[N] shape;
  simplex[N] theta;
}  

transformed parameters {
  vector<lower=0>[N] scale;
  
  // get scale from mean (Weibull)
  for(n in 1:N){
    scale[n] = exp(log(mu_step[n]) + lgamma(1 + 1/shape[n]));
  }
  
}

model {
  vector[N] log_theta = log(theta); 
  // priors
  rho ~ beta(1, 1);
  mu[1] ~ normal(pi(),0.5);
  mu[2] ~ normal(0,0.5);
  shape ~ gamma(12, 6);
  mu_step ~ normal(0, 5);
    
    // likelihood computation
    for (t in 1:T) {
      if(steps[t]>=0 && turns[t]>-pi() ){
        vector[N] lps = log_theta;
        for(n in 1:N){
          lps[n] += weibull_lpdf(steps[t] | shape[n], scale[n]) +
            wrappedCauchy_lpdf(turns[t] | mu[n], rho[n]);
        }
      target += log_sum_exp(lps);
    }
  }
}

generated quantities {
    matrix[T,N] stateProbs;
    real log_p[T,N];
    vector[N] log_theta = log(theta);
    for(t in 1:T){
      for(n in 1:N){
        stateProbs[t,n] = negative_infinity();
        log_p[t,n] = 0;
        if(steps[t]>=0 && turns[t]>-pi() ){
          log_p[t,n] = log_theta[n] +
          weibull_lpdf(steps[t] | shape[n], scale[n]) +
          wrappedCauchy_lpdf(turns[t] | mu[n], rho[n]);
        }
      }
      for(i in 1:N){
        if(steps[t]>=0 && turns[t]>-pi()){
          stateProbs[t,i] = exp(log_p[t,i])/(exp(log_sum_exp(log_p[t,])));
        }
        }
      }
    }

