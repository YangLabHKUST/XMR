library(XMR)
exposure <- "LDLC"
outcome  <- "MI"
threshold <- 5e-05

# These bundled data contain the LDLC-MI background matrices and selected IVs.
# See the HTML tutorial and workflow.R for data preparation.

data(C)
data(Omega)
data(clumped_data)

XMR_res <- fit_XMR(
  data = clumped_data, 
  C = C, 
  Omega0 = Omega,
  Threshold = threshold,
  tol1 = 1e-07, 
  tol2 = 1e-07, 
  min_thres = 1e-2
)