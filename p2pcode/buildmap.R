overlap <- function(x,y,r, xn,yn,rn, tol=0){
  dst <- sqrt( (x-xn)^2 + (y-yn)^2 )
  sum(dst < (r+rn+tol))
} 

buildmap <- function(npatches=50, tol=2, minr=0.01, maxr=0.25, sc=10){
  x <- numeric(npatches)
  y <- numeric(npatches)
  r <- numeric(npatches)
  r[1] <- runif(1,minr,maxr) #sample(c(0.01, 0.2),1,prob=c(0.75,0.25)) #
  
  count = 1
  
  for(i in 2:npatches){
    idx <- sample(1:count, 1) 
    ds <- rweibull(1, shape=2, scale=sc)
    an <- runif(1,0, 2*pi)
    xn <- x[idx] + cos(an)*ds
    yn <- y[idx] + sin(an)*ds
    rn <- runif(1,minr,maxr)
    
    while(overlap(x[1:i],y[1:i],r[1:i], xn,yn,rn, tol)>0){
      idx <- sample(1:count, 1) 
      ds <- rweibull(1, shape=2, scale=sc)
      an <- runif(1,0, 2*pi)
      xn <- x[idx] + cos(an)*ds
      yn <- y[idx] + sin(an)*ds
      rn <- runif(1,minr,maxr)
    }
    count = count+1
    x[i] = xn
    y[i] = yn
    r[i] = rn
  } 
  return(data.frame(x,y,r))
}

