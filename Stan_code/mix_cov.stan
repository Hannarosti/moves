// j.m.morales June 2018
// priors set according to one simulated example!

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
    int nCovs; // number of covariates
    matrix[T, nCovs] covs; // covariates
    real lb; // lower bound for shape
}

parameters {
  vector<lower=-1, upper=4.5>[N] mu; // this allows mu to be centered around pi
  vector<lower=0, upper=1>[N] rho;
  positive_ordered[N] mu_step;
  vector<lower=lb>[N] shape;
  //simplex[N] p;
  real alpha;
  vector[nCovs] beta;
}  

transformed parameters {
  vector<lower=0>[N] scale;
  vector[N] log_p[T];
  vector[T] p;
  
  // get scale from mean (Weibull)
  for(n in 1:N){
    scale[n] = exp(log(mu_step[n]) + lgamma(1 + 1/shape[n]));
  }
  
  p = inv_logit(covs * beta + alpha);
  
  for(t in 1:T){
    log_p[t,1] = log(1-p[t]);
    log_p[t,2] = log(p[t]);
  }
  
}

model {
  //vector[N] log_p; 
  // priors
  rho ~ beta(1, 1);
  mu[1] ~ normal(pi(),0.5);
  mu[2] ~ normal(0,0.5);
  shape ~ gamma(12, 6);
  mu_step[1] ~ normal(1, 1);
  mu_step[2] ~ normal(5, 1);
  alpha ~ normal(1,2);
  beta ~ normal(0,1);
    // likelihood computation
    for (t in 1:T) {
      if(steps[t]>=0 && turns[t]>-pi() ){
        vector[N] lps = log_p[t];
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
    real lp[T,N];
//    vector[N] log_p = log(p);
    for(t in 1:T){
      for(n in 1:N){
        stateProbs[t,n] = negative_infinity();
        lp[t,n] = 0;
        if(steps[t]>=0 && turns[t]>-pi() ){
          lp[t,n] = log_p[t,n] +
          weibull_lpdf(steps[t] | shape[n], scale[n]) +
          wrappedCauchy_lpdf(turns[t] | mu[n], rho[n]);
        }
      }
      for(i in 1:N){
        if(steps[t]>=0 && turns[t]>-pi()){
          stateProbs[t,i] = exp(lp[t,i])/(exp(log_sum_exp(lp[t,])));
        }
      }
    }
  }

