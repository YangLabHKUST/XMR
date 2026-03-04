library(XMR)
exposure <- "LDLC"
outcome  <- "MI"
threshold <- 5e-05

N1 <- 343621 # EUR sample size
N2 <- 72866 # EAS sample size

t0 = abs(qnorm(threshold / 2))
dt = 0.13 / (sqrt(N2 / N1))
modified_threshold = 2 * (1 - pnorm(abs(t0 + dt))) # modified threshold to correct selectio bias

data(C)
data(Omega)
data(clumped_data)

XMR_res <- fit_XMR(
  data = clumped_data, 
  C = C, 
  Omega0 = Omega,
  Threshold = modified_threshold,
  tol1 = 1e-07, 
  tol2 = 1e-07, 
  min_thres = 1e-2
)