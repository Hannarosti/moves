p2psims <- function(map, nind = 30, kappa = 0, m = 1e-05, mscale = 0.01, per = 0.01){ 
  
  x  <- map$x
  y  <- map$y
  rs <- map$r
  xv <- t(rbind(map$x, map$y))
  distance <- as.matrix(dist(xv, upper = TRUE, diag = TRUE)) 
  
  mshape = 2
 # mscale = 0.01
  gshape = 1
  gscale = 10
  meanspeed = 0.1
  # per = mscale * gamma(1 + 1/mshape) # 0.01  #1e-2
  
  zs <- NULL
  cts <- NULL
  id <- NULL
  tt <- NULL
  
  for(i in 1:nind){
    id <- c(id,i)
    tt <- c(tt, 0) # travel time
    cT <- 0
    TT <- 0
    cts <- c(cts,0) #numeric(nm) # cummulative time
    nm <- ceiling(rexp(1,m)) + 1
    R <- rweibull(nm, shape=mshape, scale=mscale)
    Tn <- rvonmises(nm, mu = circular(0), kappa = kappa)
    X <- numeric(nm)
    Y <- numeric(nm)
    B <- numeric(nm)
    #  Z <- numeric(nm)
    
    ini <- sample(1:npatches, 1, replace=TRUE) # starting at a random patch
    init <- runif(1, 0,2*pi)             # initial movement direction
    B[1] <- init
    X[1] <- x[ini] + cos(init) * rs[ini] 
    Y[1] <- y[ini] + sin(init) * rs[ini]
    zs <- c(zs,ini)
    for(j in 2:nm){
      B[j] <- (B[j-1] + Tn[j]) %% (2*pi) # keep directions between 0 and 2*pi
      X[j] <- X[j-1] + R[j]*cos(B[j])
      Y[j] <- Y[j-1] + R[j]*sin(B[j])
      #Z[j] <- Z[j-1]
      cT <- cT + meanspeed*R[j]
      TT <- TT + meanspeed*R[j]
      dst <- sqrt( (x-X[j])^2 + (y-Y[j])^2 )
      dst <- dst - (rs + per)
      if(any(dst <= 0)){
        idx <- which(dst==min(dst))
        if(idx != zs[length(zs)]){
          zs <- c(zs,idx)
          id <- c(id,i)
          tmp <- runif(1,0,2*pi)
          X[j] <- x[idx] + cos(tmp) * rs[idx] 
          Y[j] <- y[idx] + sin(tmp) * rs[idx] 
          B[j] <- tmp
          cts <- c(cts,cT)
          tt <- c(tt, TT)
          TT <- 0
          cT <- cT + rgamma(1,shape=gshape,scale=gscale)
        }
      }  
    }
    zs <- c(zs,npatches+1)
    cts <- c(cts, 0)
    tt <- c(tt, 0)
    id <- c(id,i)
  }
  
  return(data.frame(id, zs, cts, tt))
} 