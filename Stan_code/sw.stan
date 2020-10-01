// j.m.morales June 2018
// part of this code came from Ben Lambert, Theo Micehlot and Vianey Leos-Barajas

functions{
  //Function that returns the log pdf of the wrapped-Cauchy
  real wrappedCauchy_lpdf(real y, real mu, real rho) {
    return(- log(2*pi()) + log((1-rho^2)/(1+rho^2-2*rho*cos(y-mu))));
  }
}

data {
    int<lower=0> T; // length of the time series
    int ID[T]; // track identifier
    vector[T] steps; // step lengths
    vector[T] turns; // turning turns
    int<lower=1> N; // number of states
    int nCovs; // number of covariates
    matrix[T,nCovs+1] covs; // covariates
    real lb; // lower bound for shape
}

parameters {
  vector<lower=-1, upper=5>[N] mu; // this allows mu to be centered around pi
  vector<lower=0, upper=1>[N] rho;
  positive_ordered[N] mu_step;
  vector<lower=lb>[N] shape;
  // regression coefficients for transition probabilities
  matrix[N*(N-1),nCovs+1] beta;
}  

transformed parameters {
  vector<lower=0>[N] scale;
  matrix[N,N] log_theta[T];
  matrix[N,N] log_theta_tr[T];
  
  // get scale from mean (Weibull)
  for(n in 1:N){
    scale[n] = exp(log(mu_step[n]) + lgamma(1 + 1/shape[n]));
  }
  
  // derive array of (log-)transition probabilities
  {
      matrix[N,N] theta[T];
        for(t in 1:T) {
            int aCount = 1;
            for(n_from in 1:N) {
                for(n in 1:N) {
                    if(n_from==n) {
                        theta[t,n_from,n] = 1;
                    } else {
                        theta[t,n_from,n] = exp(beta[aCount] * to_vector(covs[t]));
                        aCount = aCount + 1;
                    }
                }
            }
            for(n_from in 1:N)
                log_theta[t][n_from] = log(theta[t][n_from]/sum(theta[t][n_from]));
        }
      
        // transpose
        for(t in 1:T)
            for(n_from in 1:N)
                for(n  in 1:N)
                    log_theta_tr[t, n, n_from] = log_theta[t,n_from,n];
    }
}

model {
  vector[N] lp; // for forward variables
  vector[N] lp_p1; // for forward variables
  
  // priors
  rho ~ beta(2, 2);
  mu[1] ~ normal(pi(), 0.25);
  mu[2] ~ normal(0, 0.25);
  shape ~ gamma(12, 6);
  mu_step ~ normal(0, 5);
    
    for(i in 1:(N*(N-1))) {
        beta[i,1] ~ normal(-2,2);
        if(nCovs>0) {
            for(j in 2:(nCovs+1))
                beta[i,j] ~ normal(0,1);
        }
    }
    
    // likelihood computation
    for (t in 1:T) {
        if(t==1 || ID[t]!=ID[t-1])
            lp = rep_vector(-log(N), N); // init dist = rep(1,N)/N
            
        for (n in 1:N) {
            lp_p1[n] = log_sum_exp(to_vector(log_theta_tr[t,n]) + lp);
            if(steps[t]>=0)
                lp_p1[n] = lp_p1[n] + weibull_lpdf(steps[t] | shape[n], scale[n]);
            if(turns[t]>=(-pi()))
                 lp_p1[n] = lp_p1[n] + wrappedCauchy_lpdf(turns[t] | mu[n], rho[n]); 

        }
        lp = lp_p1;
        
        if(t==T || ID[t+1]!=ID[t])
            target += log_sum_exp(lp);
    }
}

generated quantities {
    int<lower=1,upper=N> viterbi[T];
    real stateProbs[T,N];
    vector[N] lp;
    vector[N] lp_p1;
    vector[N] g;
    
    for(i in 1:(N*(N-1))) {
        g[i] = exp(beta[i,1]);
    }

    // Viterbi algorithm (most likely state sequence)
    {
        real max_logp;
        int back_ptr[T, N];
        real best_logp[T, N];

        for (t in 1:T) {
            if(t==1 || ID[t]!=ID[t-1]) {
                for(n in 1:N)
                    best_logp[t, n] = weibull_lpdf(steps[t] | shape[n], scale[n]);
            } else {
                for (n in 1:N) {
                    best_logp[t, n] = negative_infinity();
                    for (j in 1:N) {
                        real logp;
                        logp = best_logp[t-1, j] + log_theta[t,j,n];
                        if(steps[t]>0)
                            logp = logp + weibull_lpdf(steps[t] | shape[n], scale[n]);
                        if(turns[t]>(-pi()))
                            logp = logp + wrappedCauchy_lpdf(turns[t] | mu[n], rho[n]);

                        if (logp > best_logp[t, n]) {
                            back_ptr[t, n] = j;
                            best_logp[t, n] = logp;
                        }
                    }
                }
            }
        }

        for(t0 in 1:T) {
            int t = T - t0 + 1;
            if(t==T || ID[t+1]!=ID[t]) {
                max_logp = max(best_logp[t]);

                for (n in 1:N)
                    if (best_logp[t, n] == max_logp)
                        viterbi[t] = n;
            } else {
                viterbi[t] = back_ptr[t+1, viterbi[t+1]];
            }
        }
    }

    // forward-backward algorithm (state probabilities)
    {
        real logalpha[T,N];
        real logbeta[T,N];
        real llk;

        // log alpha probabilities
        for(t in 1:T) {
            if(t==1 || ID[t]!=ID[t-1]) {
                for(n in 1:N)
                    lp[n] = -log(N);
            }

            for (n in 1:N) {
                lp_p1[n] = log_sum_exp(to_vector(log_theta_tr[t,n]) + lp);
                if(steps[t]>=0)
                    lp_p1[n] = lp_p1[n] + weibull_lpdf(steps[t] | shape[n], scale[n]);
                if(turns[t]>=(-pi())) {
                    lp_p1[n] = lp_p1[n] + wrappedCauchy_lpdf(turns[t] | mu[n], rho[n]);
                }
                logalpha[t,n] = lp_p1[n];
            }
            lp = lp_p1;
        }

        // log beta probabilities
        for(t0 in 1:T) {
            int t = T - t0 + 1;

            if(t==T || ID[t+1]!=ID[t]) {
                for(n in 1:N)
                    lp_p1[n] = 0;
            } else {
                for(n in 1:N) {
                    lp_p1[n] = log_sum_exp(to_vector(log_theta_tr[t+1,n]) + lp);
                    if(steps[t+1]>=0)
                        lp_p1[n] = lp_p1[n] + weibull_lpdf(steps[t+1] | shape[n], scale[n]);
                    if(turns[t+1]>=(-pi()))
                        lp_p1[n] = lp_p1[n] + wrappedCauchy_lpdf(turns[t+1] | mu[n], rho[n]);
                }
            }
            lp = lp_p1;
            for(n in 1:N)
                logbeta[t,n] = lp[n];
        }

        // state probabilities
        for(t0 in 1:T) {
            int t = T - t0 + 1;
            if(t==T || ID[t+1]!=ID[t])
                llk = log_sum_exp(logalpha[t]);
            for(n in 1:N)
                stateProbs[t,n] = exp(logalpha[t,n] + logbeta[t,n] - llk);
        }
    }
}
