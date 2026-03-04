# The function for calculating evidence lower bound
cal_elbo <- function(A, pi0, m, Pi, S, inv.S, Var, inv.Var, muj, Sigmaj, hat.b, C_1, C_2, Var_q){
  
      # Elog(p(\hat gamma_j, \hat Gamma_j | \gamma_j, \alpha_j, Z_j))

      - 3 * log(2*pi) * m - 1/2*sum(log(unlist(lapply(S, det)))) -
      1/2 *sum(mapply(function(M1, vec1){t(vec1) %*% M1 %*% vec1}, vec1=hat.b, M1=inv.S)) +
      sum(Pi *mapply(function(M1, vec1, vec2){t(vec1) %*% M1 %*% A %*% vec2}, vec1=hat.b, M1=inv.S, vec2=muj)) -
      1/2 * sum(Pi*mapply(function(M1, vec1, M2){sum(diag(t(A) %*% M1 %*% A %*% (vec1 %*% t(vec1) + M2)))},
                          vec1=muj, M1=inv.S, M2=Sigmaj)) +
      # E(log p((Z_j))
      sum(Pi * log(pi0 + (pi0==0)) + (1-Pi)*log(1-pi0 + (pi0==1))) +
    
      # E(log p(gamma_j, alpha_j))

      - 1/2*m*log(det(Var)) -
      1/2 * sum(Pi * mapply(function(vec1, M2){sum(diag(inv.Var %*% (vec1 %*% t(vec1) + M2)))},
                            vec1=muj, M2=Sigmaj)) -

      1/2 * sum(1-Pi)*sum(diag(inv.Var %*% Var_q)) +
    
      # -E(log(q(Z_j)))
      sum(-Pi * log(Pi + (Pi==0)) - (1-Pi) * log(1-Pi + (Pi==1))) +
      # -E(log(q(gamma_j,alpha_j|Z_j)))

      3/2 * (1 + log(2*pi)) * m + 
      sum(Pi * (1/2 *log(unlist(lapply(Sigmaj, det))))) +
      sum((1-Pi) * (1/2 *log(det(Var_q)))) -
      sum(Pi *log(2*pnorm(C_1)) + (1-Pi)*log(2*pnorm(C_2)))
  
}


XMR_EM_func <- function(data = NULL,
                        beta=NULL,
                        SigmaX = NULL,
                        tau.sq = NULL,
                        pi0 = NULL,
                        fix.beta = F,
                        fix.tau = F,
                        fix.SigmaX = F,
                        C = diag(3),
                        Omega0 = matrix(0,3,3),
                        tol=1e-06,
                        Threshold=1,
                        ELBO=F,
                        min_thres=1e-4,
                        max_iter=5000){
  
  if(is.null(data)){
    message("No data for MR testing")
    return(NULL)
  }
  
  
  m=nrow(data)
  Cr <- abs(qnorm(Threshold/2))
  
  if(is.null(data$L2.pop1)) data$L2.pop1 = matrix(1, nrow=m, ncol=1)
  if(is.null(data$L2.pop2)) data$L2.pop2 = matrix(1, nrow=m, ncol=1)
  if(is.null(data$L12)) data$L12 = matrix(1, nrow=m, ncol=1)
  
  data$L2.pop1 = ifelse(data$L2.pop1<1, 1, data$L2.pop1)
  data$L2.pop2 = ifelse(data$L2.pop2<1, 1, data$L2.pop2)
  
  LDsc = data[, c("L2.pop1", "L12", "L12", "L12", "L2.pop2", "L2.pop2","L12", "L2.pop2", "L2.pop2")]
  
  # genome-wide shared
  Omega = matrix(unlist(LDsc), nrow=m, ncol=9) * (matrix(1, nrow=m, ncol=1) %*% matrix(as.vector(Omega0), nrow = 1,ncol=9))
  Omega = lapply(split(Omega, row(Omega)), function(x){matrix(x, nrow=3, ncol=3)})
  
  se = as.matrix(data.frame(data[, c("se.exp.pop1", "se.exp.pop2", "se.out.pop2")]))
  se.diags = lapply(split(se, row(se)), function(x){diag(x)})
  S =  mapply(function(M1, M2){M1 + M2 %*% C %*% M2}, M1 = Omega, M2 = se.diags)
  S = lapply(split(S, col(S)), function(x){matrix(x, nrow=3, ncol=3)})
  inv.S <- lapply(S, solve)
  S11 = data$se.exp.pop1^2*C[1,1] + data$L2.pop1*Omega0[1,1]
  
  # initialize
  if(is.null(beta)){
    beta = 0
  }
  
  if(is.null(SigmaX)) SigmaX = var(data[, c("b.exp.pop1", "b.exp.pop2")])
  if(is.null(tau.sq)) tau.sq = var(data[, c("b.out.pop2")])
  
  if(is.null(pi0)) pi0 = 0.5
  
  V1 = matrix(c(0,0,0,0,0,0,0,1,0), 3, 3, byrow = T)
  
  hat.b = t(as.matrix(data[,c("b.exp.pop1","b.exp.pop2","b.out.pop2")]))
  
  hat.b = lapply(split(hat.b, col(hat.b)), function(x){c(x)})
  
  likelis = NULL
  elbos = NULL
  elbo_comps = NULL
  
  s11 =  data$se.exp.pop1^2 
  C_2 = - Cr*sqrt(s11)/sqrt(S11)
  
  ###############################A
  sigma.sq = drop(SigmaX[1,1])
    
  Var =  as.matrix(Matrix::bdiag(SigmaX, tau.sq))
  inv.Var = as.matrix(Matrix::bdiag(solve(SigmaX), 1/tau.sq))
    
  A = matrix(c(1, 0, 0, 0, 1, 0, 0, beta, 1), 3, 3, byrow = T)
    
  C_1 = - Cr*sqrt(s11)/sqrt(S11 + sigma.sq)
    
  for(i in 1:max_iter){
    
    ###################################### E-step #########################################
    
    #cat("Round ", i, ":")
      
    inv.Sigmaj <-  mapply(function(M1, M2){t(A) %*% M1 %*% A + inv.Var}, M1=inv.S)
    
    inv.Sigmaj <- lapply(split(inv.Sigmaj, col(inv.Sigmaj)), function(x){matrix(x, nrow=3, ncol=3)})
     
    Sigmaj <- lapply(inv.Sigmaj, solve)
    
    muj <- mapply(function(M1, M2, vec1){M2 %*% t(A) %*% M1 %*% vec1},
                  vec1=hat.b,
                  M1=inv.S,
                  M2=Sigmaj)
    
    muj = lapply(split(muj, col(muj)),
                  function(x){matrix(x, nrow=3, ncol=1)})
      
    Var_q = Var
    
    # E_p(Z_j)
    if(pi0 == 1) Pi = rep(1, m)
    
    if(pi0 !=1)  {
        
      bj = 1/2 *mapply(function(M1, vec1){ t(vec1) %*% M1 %*% vec1}, vec1 = muj, M1=inv.Sigmaj) +
        log(pi0/(1-pi0)) + 1/2* mapply(function(M1, M2){ log(det(M1)) - log(det(Var))}, M1 = Sigmaj) +
        log(pnorm(C_2)/pnorm(C_1))
      
      Pi = 1/(1+exp(-bj))
      
      Pi = ifelse(Pi<min_thres, min_thres, Pi)
      Pi = ifelse(Pi>0.9999, 0.9999, Pi)
      
    }
    
    likeli <- cal_likeli(hat.b, A, Var, S, C_1, C_2, pi0)
    likelis <-  c(likelis, likeli)
    # #if(!ELBO){
    if(i>1 && abs((likelis[i]-likelis[(i-1)])/likelis[(i-1)]) < tol)  break
    # #}else{
    elbo <- cal_elbo(A, pi0, m, Pi, S, inv.S, Var, inv.Var, muj, Sigmaj, hat.b, C_1, C_2, Var_q)
    elbos <-  c(elbos, elbo)
      
    elbo_comps <-  c(elbo_comps, elbo)   
      
    #if(i>1 && abs((elbos[i]-elbos[(i-1)])/elbos[(i-1)]) < tol)  break
    
    #cat("beta", beta, "\t", "pi", pi0, "\t" ,"sigma.sq", sigma.sq, "\t", "tau.sq", tau.sq, "\n")

    #cat("Marginal Likelihood: ", likeli, "\t", "ELBO:", elbo,"\t", "Difference: ",likeli-elbo ,"\n")
    #cat("beta", beta, "\t", "pi", pi0, "\t" ,"sigma.sq", sigma.sq, "\t", "tau.sq", tau.sq, "\n")

    ############################################## M-step #################################################################
    if(fix.beta==F){
      
      # update beta
      temp1 = sum(Pi * mapply(function(vec1, M1,vec2)
        return(t(vec1)%*% t(V1) %*% M1 %*% vec2),
        vec1=muj, M1=inv.S, vec2= hat.b)) -
        sum(Pi * mapply(function(M1,M2,vec1)
          return(sum(diag(t(V1) %*% M1 %*% (M2 + vec1 %*% t(vec1))))),
          M1 =inv.S , M2 = Sigmaj, vec1 = muj))
      
      temp2 =
        sum(Pi * mapply(function(M1, M2, vec1)
          return(sum(diag(t(V1) %*% M1 %*% V1 %*% (M2 + vec1 %*% t(vec1))))),
          M1 =inv.S , M2 = Sigmaj, vec1 = muj))
      
      beta = temp1/temp2
      
    }
    
    ## update pi0 when pi0!=1
    if(pi0!=1){
      
      pi0 = sum(Pi)/m
      # set lower bound and upper bound of $\pi_0$ to avoid log(0)
      if(pi0<min_thres)  pi0 = min_thres
      if(pi0>0.9999)  pi0 = 0.9999
      
    }
    
    ## update sigma^2
    # update SigmaX
    muj.X <- lapply(muj, function(x) {x[1:2]})
    Sigmaj.X <- lapply(Sigmaj, function(x) {x[(1:2),(1:2)]})
    
    if(!fix.SigmaX){
      
      temp0 = mapply(function(vec1, M1) vec1 %*% t(vec1) + M1,
                     vec1 = muj.X, M1 = Sigmaj.X)
      
      tempL = matrix(temp0 %*% (Pi), byrow=F, nrow = 2, ncol =2)
      
      ### i.e. t==0
      if(Threshold==1) {
          
        SigmaX = (tempL + sum(1 - Pi) * SigmaX)/m
          
      }
      
      if(Threshold!=1) {
        
        E = matrix(0, 2, 2)
        E[1,1] = 1
        
        ######################
          
        tempR = m * solve(SigmaX) +  sum(Pi * dnorm(C_1)/pnorm(C_1) * Cr *
                                             sqrt(s11) * (S11 + sigma.sq)^(-3/2)) * E
          
        L = t(chol(tempR))
        inv.L = solve(L)
        
        SigmaX = t(inv.L) %*% expm::sqrtm(t(L) %*% (tempL + sum(1 - Pi) * SigmaX) %*% L) %*% inv.L
      }
    }
    
    muj.Y <- lapply(muj, function(x) {x[3]})
    Sigmaj.Y <- lapply(Sigmaj, function(x) {x[3,3]})
    
    ############################
    if(!fix.tau)  tau.sq = sum(Pi * (unlist(muj.Y)^2 + unlist(Sigmaj.Y)) + (1 - Pi) * tau.sq)/m
      
    sigma.sq = drop(SigmaX[1,1])
    Var =  as.matrix(Matrix::bdiag(SigmaX, tau.sq))
    inv.Var = as.matrix(Matrix::bdiag(solve(SigmaX), 1/tau.sq))
    A = matrix(c(1, 0, 0, 0, 1, 0, 0, beta, 1), 3, 3, byrow = T)
    C_1 = - Cr*sqrt(s11)/sqrt(S11 + sigma.sq)
      
    elbo_comps <-  c(elbo_comps, cal_elbo(A, pi0, m, Pi, S, inv.S, Var, inv.Var, muj, Sigmaj, hat.b, C_1, C_2, Var_q))  
    
  }
  
  muj.x <- lapply(muj, function(x) {x[2]})
  Sigmaj.x <- lapply(Sigmaj, function(x) {x[2,2]})
  
  return(list(beta = beta,
              SigmaX=SigmaX,
              tau.sq =tau.sq,
              pi0 = pi0,
              fix.tau = fix.tau,
              post = list(mu = muj.x, Pi = Pi, IVsignal.sum =  sum(Pi * unlist(muj.x)^2 + Pi*unlist(Sigmaj.x))),
              likeli = likeli,
              likelis = likelis,
              elbos = elbos,
              elbo_comps = elbo_comps,
              log_elbo=elbos[length(elbos)],
              Threshold = Threshold))
}




