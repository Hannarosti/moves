p2pmcmc <- function(like, npars, niter = 4000, adapt = 1000, adapt1 = 2000, 
                    nch = 3, maxtry = 1, opt = 0){
  
  primu <- rep(0, npars)  # mean for all priors
  prisd <- rep(1, npars)  # SD for all priors
  
  # a function to calculate sum of log priors
  pri <- function(pars) sum(dnorm(pars, mean=primu, sd=prisd, log=TRUE))
  
  p.s  <- array(0, c(niter, nch, npars))   # matrix to hold posterior values
  liks <- array(0, c(niter, nch))  # vector to hold likelihood values
  pos  <- array(0, c(niter, nch))  # vector to hold posterior values
  
  # initial values for proposal distributions
  mu <- rep(0, npars)
  ac <- array(0, c(nch, npars, 2))      # for acceptance rates
  kk <- array(0.5, c(nch, npars))       # sd for the proposal distributions
  
  # initial value for multivariate sigma
  prop_sigma <- array(0, c(nch, npars, npars))
  for(i in 1:nch) prop_sigma[i, , ] <- diag(npars)
  
  #-------- Initiate Chains -----------------------------------------------------#
  li <- numeric(nch) # to hold current likelihood
  
  # define initial values for each parameter
  for(j in 1:nch){
    p.s[1, j, ] = rnorm(npars, primu, prisd)
  }
  
  for(j in 1:nch){
    li[j] <- like(newp = p.s[1, j, ]) 
    while(!is.finite(li[j])){ # if li is not definied try other combinations
      p.s[1, j, ] <- rnorm(npars, primu, prisd)
      li[j] <- like(newp = p.s[1, j, ])
    }
    pos[1, j] <- li[j] + pri(pars = p.s[1, j, ]) 
  }
  
  if(opt != 1){
    # create progress bar
    pb <- txtProgressBar(min = 0, max = niter, style = 3)
  }
  #--------------------------------- MCMC run -----------------------------------#
  for(i in 2:niter){
    for(j in 1:nch){
      p.s[i, j, ] <- p.s[i-1, j, ]
      liks[i, j]  <- liks[i-1, j ]
      pos[i, j]   <- pos[i-1, j ]
      
      if(i < adapt){
        for(k in 1:npars){
          newp <- p.s[i, j, ]
          ac[j, k, 1] <- ac[j, k, 1] + 1
          dvs <- rnorm(1, mean = 0, sd = kk[j, k])
          newp[k] <- newp[k] + dvs
          nli <- like(newp)
          log_alpha <- nli + pri(pars = newp) - li[j] - pri(pars = p.s[i, j, ])
          if(is.finite(nli) & log(runif(1)) < log_alpha){
            ac[j, k, 2] <- ac[j, k, 2] + 1
            p.s[i, j, ] <- newp
            li[j] <- nli
            liks[i, j] <- nli
            pos[i, j]  <- nli + pri(pars = newp)
          }
        }
      } else{
        ac[j, , 1] <- ac[j, , 1] + 1   
        newp <-  p.s[i-1, j, ] + mvrnorm(1, mu, Sigma = prop_sigma[j, , ])
        nli <- like(newp)
        log_alpha <- nli + pri(pars = newp) - li[j] - pri(pars = p.s[i, j, ])
        if(is.finite(nli) & log(runif(1)) < log_alpha){
          ac[j, , 2] <- ac[j, , 2] + 1
          p.s[i, j, ] <- newp
          li[j] <- nli
          liks[i, j] <- nli
          pos[i, j] <- nli + pri(pars = newp)
        }
      }
      
      if (i< adapt & i%%50 == 0){ # tune sd of proposal distr
        arate <- ac[j, , 2] / ac[j, , 1]
        ac[j, , 1] <- rep(0, npars)
        ac[j, , 2] <- rep(0, npars)
        kk[j, ] <- kk[j, ] * 2^(arate - 0.23)
      }
      if(i >= adapt & i%%20 == 0 & i < adapt1){ # tune multivariate sigma
        len <- floor(i * 0.5):i
        x <- p.s[len, j, ]
        N <- length(len)
        p_sigma <- (N - 1) * var(x) / N
        p_sigma <- makePositiveDefinite(p_sigma)
        if (!(0 %in% p_sigma)) prop_sigma[j, , ] <- p_sigma
      }
    }
    
    if(opt == 1){
      op=par(mfrow = c(3, 3), mar = c(3, 2, 2, 1))
      for(jj in 1:npars){
        plot((p.s[1:i, 1, jj]),type = "l", ylab = "", 
             ylim = c(min(p.s[1:i, , jj]), max(p.s[1:i, , jj])))    
        for(ii in 2:nch) lines((p.s[1:i, ii, jj]), type = "l", ylab = "", col = ii)
      }
      par(op)
    } else{
      if(i%%100 == 0){
        Sys.sleep(0.1)
        # update progress bar
        setTxtProgressBar(pb, i)
      }
    }
    if(opt != 1) close(pb)
  }
  return(list(p.s[adapt1:niter,,], pos[adapt1:niter,]))
}