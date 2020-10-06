# jmm 01-17-2017
# ----------------- probability function --------------------------------------#
dp2p <- function(n, distance, radius, a1, a2, q1, q2, b1, b2, q3, square = FALSE){
  P <- matrix(NA,n,n+1)
  pp <- array(P)
  pmat <- matrix(NA,n,n)
  g <- matrix(NA,n,n)
  
  H <- besselK((a1 * radius)^q1 + (a2*distance)^q2, 0) / besselK((a1 * radius)^q1, 0)
  U <- (1 + b1 * radius) * (b2 * distance)^q3

  diag(U) <- 0
  
  if(!any( !is.finite(H))){
    ah  <- fillvec(n,H)
    au  <- fillvec(n,U)
    tmp <- fillmats(n,H,U*H)
    
    As <- sparseMatrix(tmp$ii, tmp$jj, x = tmp$zz)
    Us <- sparseMatrix(tmp$ii, tmp$jj, x = tmp$zu)
    pp <- try(solve(As, ah)) # try to solve for the p_ij
    
    # now solve for the time components
    ahu <- ah * au - Us %*% pp
    pg  <- try(solve(As, ahu))
    g   <- pg/pp
    
    G <- matrix(as.numeric(g), n, n, byrow=TRUE)
    diag(G) = 0
    pmat <- matrix(as.numeric(pp), n, n, byrow=TRUE)
    
    # discard if results don't make sense
    if(any(rowSums(pmat) > 1 )) pmat <- matrix(NA, n, n)
    if(any(pp < 0 )) pmat <- matrix(NA, n, n)
    
    P <- cbind(pmat, pmax(0, 1 - rowSums(pmat))) # the last row allows for 
    #the possibility of dying or emigrating from the patch network
    if(square == TRUE) P <- pmat / rowSums(pmat)
  }
  return(list(P, G))
}