#' Quality control of background parameter matrices
#'
#' Reshape the output of \code{estimate_bgparas()} into symmetric Omega and C
#' matrices, check for negative variances and positive semi-definiteness of C,
#' and truncate off-diagonal correlations whose absolute value exceeds 1.
#'
#' @param bgparas A \code{data.frame} returned by \code{estimate_bgparas()},
#'   containing columns \code{trait1}, \code{trait2}, \code{Omega12}, and \code{c12}.
#' @param exposure_code Character. Exposure trait code, e.g. \code{"30780"}.
#' @param outcome_code Character. Outcome trait code, e.g. \code{"411.2"}.
#' @param pop_ref Character. Reference (auxiliary) population label, e.g. \code{"EUR"}.
#' @param pop_target Character. Target population label, e.g. \code{"BBJ"}.
#' @param truncation_corr Numeric. Value used to cap off-diagonal correlations
#'   that exceed 1 in absolute value. Default \code{0.95}.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{Omega0}{A named symmetric matrix of background genetic covariance estimates.}
#'   \item{C}{A named symmetric matrix of LDSC intercept estimates.}
#' }
#'
#' @importFrom reshape2 acast
#' @export
qc_bgparas <- function(bgparas,
                        exposure_code,
                        outcome_code,
                        pop_ref,
                        pop_target,
                        truncation_corr = 0.95) {

  ## ---- reshape to symmetric form ----
  paras <- bgparas[, c("trait1", "trait2", "Omega12", "c12")]
  paras_opt <- bgparas[, c("trait2", "trait1", "Omega12", "c12")]
  colnames(paras_opt) <- c("trait1", "trait2", "Omega12", "c12")
  paras_opt <- subset(paras_opt, trait1 != trait2)
  paras_all <- rbind(paras, paras_opt)

  ## ---- cast to matrices ----
  Omega0 <- reshape2::acast(
    unique(paras_all[, c("trait1", "trait2", "Omega12")]),
    trait1 ~ trait2, value.var = "Omega12"
  )
  C <- reshape2::acast(
    unique(paras_all[, c("trait1", "trait2", "c12")]),
    trait1 ~ trait2, value.var = "c12"
  )

  ## ---- reorder rows / columns ----
  idx <- c(paste0(pop_ref, "_", exposure_code),
           paste0(pop_target, "_", exposure_code),
           paste0(pop_target, "_", outcome_code))
  Omega0 <- Omega0[idx, idx]
  C      <- C[idx, idx]

  ## ---- QC: negative diagonal in Omega0 ----
  if (Omega0[1, 1] < 0 | Omega0[2, 2] < 0 | Omega0[3, 3] < 0) {
    warning("Omega0 contains negative variance on the diagonal.")
  }

  ## ---- QC: positive semi-definiteness of C ----
  flag_psd <- TRUE
  tryCatch(
    { chol(C) },
    error = function(e) { flag_psd <<- FALSE }
  )
  if (!flag_psd) {
    warning("C is not positive semi-definite.")
  }

  ## ---- QC: truncate off-diagonal correlations exceeding 1 ----
  # (1, 3) pair
  if (abs(Omega0[1, 3] / sqrt(Omega0[1, 1] * Omega0[3, 3])) > 1) {
    Omega0[1, 3] <- sqrt(Omega0[1, 1] * Omega0[3, 3]) * truncation_corr * sign(Omega0[1, 3])
    Omega0[3, 1] <- Omega0[1, 3]
  }
  # (2, 3) pair
  if (abs(Omega0[2, 3] / sqrt(Omega0[2, 2] * Omega0[3, 3])) > 1) {
    Omega0[2, 3] <- sqrt(Omega0[2, 2] * Omega0[3, 3]) * truncation_corr * sign(Omega0[2, 3])
    Omega0[3, 2] <- Omega0[2, 3]
  }
  # (1, 2) pair
  if (abs(Omega0[1, 2] / sqrt(Omega0[1, 1] * Omega0[2, 2])) > 1) {
    Omega0[1, 2] <- sqrt(Omega0[1, 1] * Omega0[2, 2]) * truncation_corr * sign(Omega0[1, 2])
    Omega0[2, 1] <- Omega0[1, 2]
  }

  return(list(Omega0 = Omega0, C = C))
}