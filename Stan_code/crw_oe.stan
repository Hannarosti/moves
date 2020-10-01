functions{
  //Function that returns the log pdf of the wrapped-Cauchy
  real wrappedCauchy_lpdf(real y, real rho, real mu) {
    return(log(1/(2*pi())) + log((1-rho^2)/(1+rho^2-2*rho*cos(y-mu))));
    }
}

data {
  int<lower=0> T;   // length of the time series
  vector[T] xo;
  vector[T] yo;
}

parameters {
  real<lower=0, upper=1> rho;
  real<lower=-pi(), upper=pi()> m;
  real<lower=0> shape;
  real<lower=0> scale;
  real<lower=0> sigma;
  vector<lower=0>[T] steps;
  vector[T] h;
  real x1;
  real y1;
  
}

transformed parameters{
  vector[T] x;
  vector[T] y;
  x[1] = x1; //~ normal(0,1);
  y[1] = y1; //~ normal(0,1);
  for(i in 2:T){
    x[i] = x[i-1] + cos(h[i]) * steps[i-1];
    y[i] = y[i-1] + sin(h[i]) * steps[i-1];
    }
}

model {
  h[1] ~ uniform(0, 2*pi());
  x1 ~ normal(xo[1], sigma);
  y1 ~ normal(yo[1], sigma);
  xo[1] ~ normal(x[1], sigma);
  yo[1] ~ normal(y[1], sigma);
  m ~ normal(0,0.5);
  rho ~ beta(4,1);
  scale ~ normal(0,1);
  shape ~ gamma(2,1);
  
  for(i in 2:T){
    steps[i-1] ~ weibull(shape, scale);
    h[i] ~ wrappedCauchy(rho, h[i-1] + m);
    }
    
    for (t in 2:T) {
      xo[t] ~ normal(x[t], sigma);
      yo[t] ~ normal(y[t], sigma);
      }
}
