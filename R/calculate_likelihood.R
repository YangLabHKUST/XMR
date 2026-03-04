# The function for calculating likelihood
cal_likeli <- function(hat.b, A, Var, S, C1, C2, pi0){

  S2 <- mapply(function(M1){M1 +  A %*% Var %*% t(A) }, M1=S, SIMPLIFY = FALSE)

  P1 =  mapply(mvtnorm::dmvnorm, x= hat.b,  sigma=S2)
  P2 =  mapply(mvtnorm::dmvnorm, x= hat.b,  sigma=S)

  #return likelihood
  return(sum(log(pi0 * P1/2/pnorm(C1) + (1-pi0) * P2/2/pnorm(C2))))

}