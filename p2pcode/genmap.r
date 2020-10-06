overlap <- function(x, y, r, xn, yn, rn, tol = 0){
 dst <- sqrt( (x-xn)^2 + (y-yn)^2 )
 sum(dst < (r+rn+tol))
} 

genmap <- function(npatches=20, tol=0.05, minr=0.01, maxr=0.2, minx=0, maxx=1){
 x <- numeric(npatches)
 y <- numeric(npatches)
 r <- numeric(npatches)

 x[1] <- runif(1,minx,maxx)
 y[1] <- runif(1,minx,maxx)
 r[1] <- runif(1,minr,maxr) #sample(c(0.01, 0.2),1,prob=c(0.75,0.25)) #

 for(i in 2:npatches){
  xn <- runif(1,minx,maxx)
  yn <- runif(1,minx,maxx)
  rn <- runif(1,minr,maxr) #s sample(c(0.01, 0.2),1,prob=c(0.75,0.25)) #runif(1,minr,maxr)

 while(overlap(x[1:i],y[1:i],r[1:i], xn,yn,rn, tol)>0){
  xn <- runif(1,minx,maxx)
  yn <- runif(1,minx,maxx)
  rn <- runif(1,minr,maxr) #s sample(c(0.01, 0.2),1,prob=c(0.75,0.25)) # runif(1,minr,maxr)
 }
 x[i] = xn
 y[i] = yn
 r[i] = rn
} 
return(data.frame(x,y,r))
}
#symbols(x,y,r, inches=F, asp=1)
