#' @importFrom stats dnorm lm median pchisq pnorm pt qchisq qnorm qt sd var weighted.mean
#' @importFrom dplyr %>%
#' @title A function for implementing XMR.
#' @description XMR: a cross-population Mendelian randomization method for
#'   causal inference using genome-wide summary statistics.
#'   XMR uses a variational EM algorithm for the estimation of parameters.
#'   XMR uses likelihood ratio test for inference.
#'
#' @param data data.frame that at least contains the following columns:
#'   b.exp.pop1, b.exp.pop2, b.out.pop2, se.exp.pop1, se.exp.pop2, se.out.pop2,
#'   L2.pop1, L12, L2.pop2. b: SNP-trait effect; se: the standard error of b
#'   (exp: exposure, out: outcome, pop1: the auxiliary population,
#'   pop2: the target population); L2.pop1: LD score in the auxiliary population;
#'   L2.pop2: LD score in the target population; L12: cross-population LD score
#'   between pop1 and pop2.
#' @param exposure exposure name.
#' @param outcome outcome name.
#' @param SigmaX initial value for SigmaX, default `NULL` will use the default
#'   initialize procedure.
#' @param tau.sq initial value for tau.sq, default `NULL` will use the default
#'   initialize procedure.
#' @param pi0 initial value for pi0, default `NULL` will use the default
#'   initialize procedure.
#' @param C the estimated C matrix capturing the effects of sample structure.
#'   Default `diag(3)`.
#' @param Omega0 the estimated variance-covariance matrix of polygenic effects.
#'   Default `matrix(0, 3, 3)`.
#' @param tol1 tolerance of the 1st fitting round fixing beta=0, default `1e-06`.
#' @param tol2 tolerance of the 2nd fitting round which fits beta, default `1e-06`.
#' @param Threshold modified IV selection threshold for correction of selection bias.
#' @param ELBO logical, whether to track the evidence lower bound or not.
#'   If `FALSE`, check the maximum likelihood instead. Default `FALSE`.
#' @param verbose logical, whether to print results. Default `TRUE`.
#' @param min_thres lower bound of pi0 to avoid log(0); pi0 smaller than this
#'   value during fitting will be truncated to it. Default `1e-4`.
#' @param max_iter maximum iteration rounds of the variational-EM algorithm.
#'   Default `5000`.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{exposure}{exposure of interest}
#'   \item{outcome}{outcome of interest}
#'   \item{beta}{causal effect estimate}
#'   \item{beta.se}{standard error}
#'   \item{beta.pvalue}{p-value}
#'   \item{tau.sq}{variance of the outcome effect in the causal inference module}
#'   \item{SigmaX}{covariance matrix of exposure effects in the causal inference module}
#'   \item{post}{posterior probability of the validity of each input IV}
#'   \item{nIV}{number of input IVs}
#'   \item{nvalid}{number of valid IVs detected by XMR}
#'   \item{pi0}{the probability of an SNP with foreground signals after selection}
#'   \item{Threshold}{input modified IV selection threshold for correction of selection bias}
#'   \item{fit1_elbos}{records of the ELBO changes in the 1st fitting round}
#'   \item{fit2_elbos}{records of the ELBO changes in the 2nd fitting round}
#'   \item{fit1_likelis}{records of the likelihood changes in the 1st fitting round}
#'   \item{fit2_likelis}{records of the likelihood changes in the 2nd fitting round}
#' }
#'
#' @examples
#' exposure <- "LDLC"
#' outcome  <- "MI"
#' threshold <- 5e-05
#' N1 <- 343621 # EUR sample size
#' N2 <- 72866  # EAS sample size
#' t0 <- abs(qnorm(threshold / 2))
#' dt <- 0.13 / (sqrt(N2 / N1))
#' modified_threshold <- 2 * (1 - pnorm(abs(t0 + dt))) # modified threshold
#' data(C)
#' data(Omega)
#' data(clumped_data)
#' XMR_res <- fit_XMR(
#'   data = clumped_data,
#'   C = C,
#'   Omega0 = Omega,
#'   Threshold = modified_threshold,
#'   tol1 = 1e-07,
#'   tol2 = 1e-07,
#'   min_thres = 1e-2
#' )
#' @export


fit_XMR <- function(data,
                    exposure="exposure",
                    outcome="outcome",
                    SigmaX = NULL,
                    tau.sq = NULL,
                    pi0 = NULL,
                    C = diag(3),
                    Omega0 = matrix(0, 3, 3),
                    tol1=1e-06,
                    tol2=1e-06,
                    Threshold=5e-05,
                    ELBO=F,
                    verbose=T,
                    min_thres=1e-4,
                    max_iter=5000){
  
  fit1 = XMR_EM_func(data = data,
                     beta = 0,
                     SigmaX = SigmaX,
                     tau.sq = tau.sq,
                     pi0 = pi0,
                     fix.beta = T,
                     C = C,
                     Omega0 = Omega0,
                     tol=tol1,
                     Threshold=Threshold,
                     ELBO=ELBO,
                     min_thres=min_thres,
                     max_iter=max_iter)
  
  fit2 = XMR_EM_func(data = data,
                     beta = 0,
                     SigmaX = fit1$SigmaX,
                     tau.sq = fit1$tau.sq,
                     pi0 = fit1$pi0,
                     fix.beta = F,
                     C = C,
                     Omega0 = Omega0,
                     tol=tol2,
                     Threshold=Threshold,
                     ELBO=ELBO,
                     min_thres=min_thres,
                     max_iter=max_iter)
  
  # Inference
  #LR1 = 2*(fit2$log_elbo - fit1$log_elbo)
  LR1 = 2*(fit2$likeli - fit1$likeli)
  pvalue = pchisq(LR1, 1, lower.tail = F)
  pvalue = formatC(pvalue, format = "e", digits = 4)
  beta.se = suppressWarnings(abs(fit2$beta/sqrt(LR1)))
  SigmaX = fit2$SigmaX
  
  if(verbose){
    cat("***********************************************************\n")
    cat("MR test results of ", exposure , " on ", outcome, ": \n")
    cat("beta = ", round(fit2$beta,4), "beta.se = ", round(beta.se, 4), "beta.pvalue = ", pvalue,  "\n")
    cat("Total No.of IVs:", nrow(data), "Effective NO. of IVs:", fit2$pi0 * nrow(data), "\n")
    cat("***********************************************************\n")
  }

  return( list(exposure = exposure,
               outcome = outcome,
               beta = round(fit2$beta,4),
               beta.se = beta.se,
               beta.pvalue = pvalue,
               tau.sq = fit2$tau.sq,
               SigmaX = fit2$SigmaX,
               post=fit2$post,
               nIV = nrow(data),
               nvalid = fit2$pi0 * nrow(data),
               pi0 = fit2$pi0,
               Threshold = Threshold,
               fit1_elbos = fit1$elbo_comps,
               fit2_elbos = fit2$elbo_comps,
               fit1_likelis = fit1$likelis,
               fit2_likelis = fit2$likelis))
}