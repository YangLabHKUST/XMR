utils::globalVariables(c("b.exp"))

ci_normal <- function(type, mean, se, alpha){
  x <- 1 - alpha/2
  
  if(type == "l") return(mean - qnorm(x)*se)
  else if (type == "u") return(mean + qnorm(x)*se)
}

ci_t <- function(type, mean, se, df, alpha){
  x <- 1 - alpha/2
  
  if(type == "l") return(mean - qt(x, df = df)*se)
  else if (type == "u") return(mean + qt(x, df = df)*se)
}

mr_conmix <- function(snps, Bx ,By, Bxse, Byse, psi = 0, CIMin  = NA, CIMax  = NA, CIStep = 0.01,
                      alpha = 0.05)
{
  
  nsnps = length(Bx)
  
  ratio = By/Bx; ratio.se = abs(Byse/Bx);
  if (is.na(CIMin)) { CIMin = min((By-2*Byse)/Bx) }
  if (is.na(CIMax)) { CIMax = max((By+2*Byse)/Bx) }
  if (psi < 0 | psi == 0) {   psi = 1.5*sd(ratio)  }
  theta = seq(from = CIMin, to = CIMax, by = CIStep)
  iters = length(theta)
  lik=NULL
  for (j1 in 1:iters) {
    lik.inc = -(theta[j1]-ratio)^2/2/ratio.se^2 - log(sqrt(2*pi*ratio.se^2))
    lik.exc = -ratio^2/2/(psi^2+ratio.se^2) - log(sqrt(2*pi*(psi^2+ratio.se^2)))
    valid = (lik.inc>lik.exc)*1
    lik[j1] = sum(c(lik.inc[valid==1], lik.exc[valid==0]))
    if (which.max(lik)==length(lik)) { valid.best = valid }
  }
  phi = ifelse(sum(valid.best)<1.5, 1,
               max(sqrt(sum(((ratio[valid.best==1]-weighted.mean(ratio[valid.best==1],
                                                                 ratio.se[valid.best==1]^-2))^2*
                               ratio.se[valid.best==1]^-2))/(sum(valid.best)-1)), 1))
  loglik = lik
  
  lik.inc0 = -ratio^2/2/ratio.se^2 - log(sqrt(2*pi*ratio.se^2))
  lik.exc0 = -ratio^2/2/(psi^2+ratio.se^2) - log(sqrt(2*pi*(psi^2+ratio.se^2)))
  valid = (lik.inc0>lik.exc0)*1
  loglik0 = sum(c(lik.inc0[valid==1], lik.exc0[valid==0]))
  
  
  whichin = which(2*loglik>(2*max(loglik)-qchisq(1-alpha, df=1)*phi^2))
  # provides an index of estimate values in the 95% confidence interval
  betaConMix = CIMin+CIStep*(which.max(loglik)-1)
  # modal estimate
  CIRange    = CIMin+CIStep*(whichin-1);
  CILower <- c(min(CIRange), CIRange[which(diff(CIRange)>1.01*CIStep)+1])
  CIUpper <- c(CIRange[which(diff(CIRange)>1.01*CIStep)], max(CIRange))
  
  Pvalue = pchisq(2*(max(loglik)-loglik0)*phi^2, df=1, lower.tail=FALSE)
  
  return(list(Psi = as.numeric(psi),
              Estimate = as.numeric(betaConMix),
              CIRange  = as.numeric(CIRange),
              CILower  = as.numeric(CILower),
              CIUpper  = as.numeric(CIUpper),
              CL = min(CIRange),
              CU = max(CIRange),
              CIMin    = as.numeric(CIMin),
              CIMax    = as.numeric(CIMax),
              CIStep   = as.numeric(CIStep),
              Valid    = as.numeric(which(valid.best==1)),
              ValidSNPs= as.character(snps[which(valid.best==1)]),
              Pvalue   = as.numeric(Pvalue),
              SNPs = nsnps,
              Alpha = alpha))
}

mr_lasso <- function(snps, Bx ,By, Bxse, Byse, distribution = "normal", alpha = 0.05, 
                     lambda = c(seq(from=0.01, to=2, by=0.01), seq(from=2.2, to=10, by=0.2))){
  nsnps = length(Bx)

  if(distribution %in% c("normal", "t-dist")){
    
    S = diag(Byse^-2)
    b = S^(1/2) %*% Bx
    Pb = b %*% solve(t(b) %*% b, t(b))
    xlas = (diag(nsnps) - Pb) %*% S^(1/2)
    ylas = (diag(nsnps) - Pb) %*% S^(1/2) %*% By
    
    if (length(lambda) != 0) {
      las_fit = glmnet::glmnet(xlas, ylas, intercept = FALSE, lambda = lambda)
      las_mod = list(fit = las_fit$beta[, 1], lambda = lambda)
    } else {
      las_fit = glmnet::glmnet(xlas, ylas, intercept = FALSE)
      lamseq = sort(las_fit$lambda)
      lamlen = length(lamseq)
      rse = sapply(1:lamlen, function(j){
        av = which(las_fit$beta[, (lamlen - j + 1)] == 0)
        mod = lm(S[av, av]^(1/2) %*% By[av] ~ S[av, av]^(1/2) %*% Bx[av] - 1)
        c(sqrt(t(mod$residuals) %*% (mod$residuals) / (mod$df.residual)), length(av))
      })
      rse_inc = rse[1, 2:lamlen] - rse[1, 1:(lamlen-1)]
      het = which(rse[1, 2:lamlen] > 1 & rse_inc > ((qchisq(0.95, 1) / rse[2, 2:lamlen])))
      if (length(het) == 0){
        lam_pos = lamlen
      } else {
        lam_pos = min(het)
      }
      num_valid = rse[2, ]
      min_lam_pos = min(which(num_valid > 1))
      if (lam_pos < min_lam_pos){lam_pos = min_lam_pos}
      las_mod = list(fit = las_fit$beta[, (lamlen - lam_pos + 1)], lambda = lamseq[lam_pos])
    }
    a = las_mod$fit
    e = By - a
    est = solve(t(Bx) %*% S %*% Bx, t(Bx) %*% S %*% e)
    v = which(a == 0)
    print(length(which(a != 0)))
    if (length(v) > 1){
      post_mod = summary(lm(By[v] ~ Bx[v] - 1, weights = Byse[v]^-2))
      post_est = post_mod$coef[, 1]
      post_se = post_mod$coef[, 2] / min(post_mod$sigma, 1)
    } else {
      post_est = NA
      post_se = NA
      cat("Specified value of lambda results in fewer than two valid instruments. Post-lasso method cannot be performed.")
    }
    
    if(distribution == "normal"){
      ciLower <- ci_normal("l", post_est, post_se, alpha)
      ciUpper <- ci_normal("u", post_est, post_se, alpha)
    } else if (distribution == "t-dist"){
      ciLower <- ci_t("l", post_est, post_se, length(v)-1, alpha)
      ciUpper <- ci_t("u", post_est, post_se, length(v)-1, alpha)
    }
    
    if (distribution == "normal") { pvalue = 2*pnorm(-abs(post_est/post_se)) }
    if (distribution == "t-dist") { pvalue = 2*pt(-abs(post_est/post_se), df=length(v)-1) }
    
    return(list(Estimate = as.numeric(post_est),
                StdError = as.numeric(post_se),
                CILower =  as.numeric(ciLower),
                CIUpper = as.numeric(ciUpper),
                Alpha = alpha,
                Pvalue = as.numeric(pvalue),
                SNPs = nsnps,
                RegEstimate = as.numeric(est),
                RegIntercept = as.numeric(a),
                Valid = length(v),
                ValidSNPs= as.character(snps[v]),
                Lambda = las_mod$lambda))
  } else {
    cat("Distribution must be one of : normal, t-dist. \n")
    cat("See documentation for details. \n")
  }
}

runMRmethods <- function(dat, set.seed = T, ldmat=NULL, Threshold=5e-08, exposure="exposure", outcome="outcome",
                         methods.list = c("IVW","dIVW", "RAPS","Egger", "MRMix",
                                          "Weighted-median", "Weighted-mode", "BWMR","MR-PRESSO","MR-Lasso"),
                         ConMixonly=FALSE){
  if(set.seed == T){
    set.seed(1234)
  }
  
  cat("Pair: ", exposure,"~", outcome,"\n")
  dat = subset(dat, b.exp!=0)
  
  indx = which(dat$pval.exp <= Threshold)
  
  nsnp = length(indx)
  
  res= NULL
  
  
  if(nsnp>3){
    
    # IVW-fixed
    
    if("IVW_fe" %in% methods.list & ConMixonly==FALSE){
      cat("IVW \n")
      res.IVW <- try(TwoSampleMR::mr_ivw_fe(dat$b.exp[indx], dat$b.out[indx], dat$se.exp[indx], dat$se.out[indx]))
      
      
      if(!inherits(res.IVW,"try-error")){
        cat("IVW-fe: beta.hat=", res.IVW$b, " se=", res.IVW$se, "pval=", res.IVW$pval, "\n")
        
        IVW_res <- data.frame(exposure = exposure, outcome = outcome, method =  "IVW_fe", Threshold = Threshold,
                              nsnp = nsnp, beta =  res.IVW$b, se = res.IVW$se, pval = res.IVW$pval)
        res = rbind(res, IVW_res)
      }
    }
    
    
    # IVW-random
    if("IVW" %in% methods.list & ConMixonly==FALSE){
      cat("IVW \n")
      res.IVW <- try(TwoSampleMR::mr_ivw(dat$b.exp[indx], dat$b.out[indx], dat$se.exp[indx], dat$se.out[indx]))
      
      
      if(!inherits(res.IVW,"try-error")){
        cat("IVW: beta.hat=", res.IVW$b, " se=", res.IVW$se, "pval=", res.IVW$pval, "\n")
        
        IVW_res <- data.frame(exposure = exposure, outcome = outcome, method =  "IVW", Threshold = Threshold,
                              nsnp = nsnp, beta =  res.IVW$b, se = res.IVW$se, pval = res.IVW$pval)
        res = rbind(res, IVW_res)
      }
    }
    
    
    # Egger
    if("Egger" %in% methods.list & ConMixonly==FALSE){
      
      res.egger <- try(TwoSampleMR::mr_egger_regression(dat$b.exp[indx], dat$b.out[indx],
                                                        dat$se.exp[indx], dat$se.out[indx]))
      if(!inherits(res.egger,"error")){
        
        cat("Egger: beta.hat=", res.egger$b, " se=", res.egger$se, "pval=", res.egger$pval, "\n")
        
        egger_res <- data.frame(exposure = exposure, outcome = outcome, method =  "Egger", Threshold = Threshold,
                                nsnp = nsnp, beta =  res.egger$b, se = res.egger$se, pval = res.egger$pval)
        res = rbind(res, egger_res)
      }
      
    }
    
    # MRmix 
    if("MRMix" %in% methods.list & ConMixonly==FALSE){
      
      res.MRMix <- try(MRMix::MRMix(dat[indx,]$b.exp, dat[indx,]$b.out, dat[indx,]$se.exp, dat[indx,]$se.out))
      
      if(!inherits(res.MRMix,"try-error")){
        
        
        cat("MRMix: beta.hat=", res.MRMix$theta, " se=", res.MRMix$SE_theta, "pval=", res.MRMix$pvalue_theta, "\n")
        
        
        MRMix_res <- data.frame(exposure = exposure, outcome = outcome, method = "MRMix", Threshold = Threshold,
                                nsnp = nsnp, beta = res.MRMix$theta, se = res.MRMix$SE_theta, pval = res.MRMix$pvalue_theta)
        
        res = rbind(res, MRMix_res)
        
      }
    }
    
    # Weighted-median
    if("Weighted-median" %in% methods.list & ConMixonly==FALSE){
      
      res.median <- try(TwoSampleMR::mr_weighted_median(dat$b.exp[indx], dat$b.out[indx], dat$se.exp[indx], dat$se.out[indx]))
      
      cat("Weighted-median: beta.hat=", res.median$b, " se=", res.median$se, "pval=", res.median$pval, "\n")
      
      if(!inherits(res.median,"try-error")) {
        median_res <- data.frame(exposure = exposure, outcome = outcome, method =  "Weighted-median", Threshold = Threshold,
                                 nsnp = nsnp, beta = res.median$b, se = res.median$se, pval = res.median$pval)
        res = rbind(res, median_res)
      }
      
    }
    
    # Weighted-mode
    if("Weighted-mode" %in% methods.list & ConMixonly==FALSE){
      
      res.mode <- try(TwoSampleMR::mr_weighted_mode(dat$b.exp[indx], dat$b.out[indx], dat$se.exp[indx], dat$se.out[indx]))
      
      cat("Weighted-mode: beta.hat=", res.mode$b, " se=", res.mode$se, "pval=", res.mode$pval, "\n")
      
      if(!inherits(res.mode,"try-error")){
        mode_res <- data.frame(exposure = exposure, outcome = outcome, method =  "Weighted-mode", Threshold = Threshold,
                               nsnp = nsnp, beta = res.mode$b, se = res.mode$se, pval = res.mode$pval)
        res = rbind(res, mode_res)
      }
      
      
    }
    

    # MR-PRESSO
    if("MR-PRESSO" %in% methods.list & ConMixonly==FALSE){
      
      NbDistribution = 1000
      if(nrow(dat[indx,]) >= NbDistribution ) NbDistribution = 1.1 * nrow(dat[indx,])
      res.presso <- try(MRPRESSO::mr_presso(BetaOutcome = "b.out", BetaExposure = "b.exp",
                                            SdOutcome = "se.out", SdExposure = "se.exp",
                                            OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                                            data = dat[indx,], NbDistribution = NbDistribution,  
                                            SignifThreshold = 0.05, seed=1234))
      save(res.presso, file="CTSD.PRESSO")
      if(!inherits(res.presso, "try-error")){
        
        if(is.na(res.presso[[1]][2,3])){
          presso_res <- data.frame(exposure = exposure, outcome = outcome, method =  "MR-PRESSO", Threshold = Threshold,
                                   nsnp = nsnp, beta = res.presso[[1]][1,3], se = res.presso[[1]][1,4], pval = res.presso[[1]][1,6]) 
          cat("MR-PRESSO: beta.hat=", res.presso[[1]][1,3], " se=", res.presso[[1]][1,4], "pval=", res.presso[[1]][1,6], "\n")
          
        }else{
          presso_res <- data.frame(exposure = exposure, outcome = outcome, method =  "MR-PRESSO", Threshold = Threshold,
                                   nsnp = nsnp, beta = res.presso[[1]][2,3], se = res.presso[[1]][2,4], pval = res.presso[[1]][2,6]) 
          cat("MR-PRESSO: beta.hat=", res.presso[[1]][2,3], " se=", res.presso[[1]][2,4], "pval=", res.presso[[1]][2,6], "\n")
          
        }
        
        res = rbind(res, presso_res)
      }
    }
    
    # MR-Robust
    if("MR-Robust"  %in% methods.list & ConMixonly==FALSE){
      
      cat("MR-Robust")
          
      rob =  try(robustbase::lmrob(dat$b.out[indx]~dat$b.exp[indx]-1, 
                                   weights=dat$se.out[indx]^-2, k.max=500))
      
      if(!inherits(rob,"try-error")) {
        fitrob = summary(rob)
        betaIVW.robust = fitrob$coef[1]
        # sebetaIVW.robust.fixed = summary(lmrob(betaYG~betaXG-1, weights=sebetaYG^-2, k.max=500))$coef[1,2]/
        #   summary(lmrob(betaYG~betaXG-1, weights=sebetaYG^-2, k.max=500))$sigma
        sebetaIVW.robust.random = fitrob$coef[1,2]/min(fitrob$sigma,1)
        IVW.robust.random_res <- data.frame(exposure = exposure, outcome = outcome, 
                                            method =  "MR-Robust", Threshold = Threshold,
                                            nsnp = nsnp, beta =  betaIVW.robust, se = sebetaIVW.robust.random,
                                            pval = pchisq(betaIVW.robust^2/sebetaIVW.robust.random^2, 1, lower.tail = F))
        res = rbind(res, IVW.robust.random_res)
        
      }
    }
    
    # MR-Lasso
    if("MR-Lasso"  %in% methods.list & ConMixonly==FALSE){
      
      Lasso_fit = try(mr_lasso(snps=dat$SNP[indx],dat$b.exp[indx], dat$b.out[indx], dat$se.exp[indx], dat$se.out[indx], lambda=NULL))
      
      if(!inherits(Lasso_fit,"try-error")) {
        cat("MR-Lasso: beta.hat=", Lasso_fit$Estimate, " pval=", Lasso_fit$Pvalue,"\n")
        
        Lasso_res = data.frame(exposure = exposure, outcome = outcome, method = "MR-Lasso",Threshold = Threshold,
                               nsnp = nsnp, beta =  Lasso_fit$Estimate, se = Lasso_fit$StdError, pval = Lasso_fit$Pvalue)
        
        res = rbind(res,  Lasso_res)
      }
    }
    
    
    # cML-MA-DP
    if("cML-MA-DP" %in% methods.list& ConMixonly==FALSE){
      final_dat = dat[indx,]
      n= min(median(1/final_dat$se.exp^2), median(1/final_dat$se.out^2))
      res_cML = try(suppressWarnings(MRcML::mr_cML_DP(final_dat$b.exp,
                                                      final_dat$b.out,
                                                      final_dat$se.exp,
                                                      final_dat$se.out,
                                                      n = n,
                                                      random_start = 10,
                                                      random_start_pert = 10,
                                                      random_seed = 1,
                                                      num_pert = 200)))
      cat("cML-MA: beta.hat=", res_cML$MA_BIC_DP_theta, " se=", res_cML$MA_BIC_DP_se, "pval=", res_cML$MA_BIC_DP_p, "\n")
      
      
      cML_res = data.frame(exposure = exposure, outcome = outcome, method =  "cML-MA", Threshold = Threshold,
                           nsnp = nsnp, beta =  res_cML$MA_BIC_DP_theta, se = res_cML$MA_BIC_DP_se, pval = res_cML$MA_BIC_DP_p)
      
      res = rbind(res, cML_res)
    }   
    
    # cML-MA
    if("cML-MA" %in% methods.list & ConMixonly==FALSE){
      final_dat = dat[indx,]
      n= min(median(1/final_dat$se.exp^2), median(1/final_dat$se.out^2))
      res_cML = try(suppressWarnings(MRcML::mr_cML(final_dat$b.exp,
                                                   final_dat$b.out,
                                                   final_dat$se.exp,
                                                   final_dat$se.out,
                                                   n = n,
                                                   random_start = 10,
                                                   random_seed = 1)))
      cat("cML-MA: beta.hat=", res_cML$MA_BIC_theta, " se=", res_cML$MA_BIC_se, "pval=", res_cML$MA_BIC_p, "\n")
      
      
      cML_res = data.frame(exposure = exposure, outcome = outcome, method =  "cML-MA", Threshold = Threshold,
                           nsnp = nsnp, beta =  res_cML$MA_BIC_theta, se = res_cML$MA_BIC_se, pval = res_cML$MA_BIC_p)
      
      res = rbind(res, cML_res)
    }   
    
    
    # MR-ConMix
    if(ConMixonly){
      
      ConMix_fit = try(mr_conmix(snps=dat$SNP[indx],dat$b.exp[indx], dat$b.out[indx], dat$se.exp[indx], dat$se.out[indx], CIMin  = -1.5, CIMax  = 1.5))
      
      if(!inherits(ConMix_fit,"try-error")) {
        cat("MR-ConMix: beta.hat=", ConMix_fit$Estimate, " pval=", ConMix_fit$Pvalue," CL=", ConMix_fit$CL, " CU=", ConMix_fit$CU,"\n")
        
        res = data.frame(exposure = exposure, outcome = outcome, method =  "MR-ConMix", Threshold = Threshold,
                         nsnp = nsnp, beta =  ConMix_fit$Estimate, pval = ConMix_fit$Pvalue, CL=ConMix_fit$CL, CU=ConMix_fit$CU)
      }
    }
  }
  
  return(res)
}