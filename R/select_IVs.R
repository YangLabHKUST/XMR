#' Select and clump candidate IVs for the XMR workflow
#'
#' Apply the auxiliary-exposure P-value threshold, rank eligible SNPs by signed
#' inverse-variance fixed-effect meta-analysis (FEMA) P values, and clump using
#' a local target-population LD reference. There is no additional FEMA P-value
#' threshold or heterogeneity filter. Effects must already be allele-aligned
#' and expressed in compatible exposure units across populations.
#'
#' @param data A harmonized data.frame containing SNP, b.exp.pop1, b.exp.pop2,
#'   se.exp.pop1, se.exp.pop2, and pval.exp.pop1. Other columns, including outcome
#'   associations and LD scores required by fit_XMR(), are preserved.
#' @param Threshold Numeric. Auxiliary-exposure P-value threshold. Default 5e-5.
#'   Pass this same unmodified threshold to fit_XMR().
#' @param bfile Character. Prefix of local target-population PLINK bed/bim/fam
#'   reference files. The population of this reference must be checked by the user.
#' @param plink_bin Character. Path to a local PLINK executable.
#' @param clump_kb Numeric. LD clumping window in kilobases. Default 10000 (10 Mb).
#' @param clump_r2 Numeric. LD clumping r-squared cutoff. Default 0.001.
#'
#' @return A data.frame of retained candidates, with the input columns plus
#'   b.exp.fema, se.exp.fema, pval.exp.fema and Threshold. An empty data.frame
#'   is returned if no candidate passes the auxiliary threshold or clumping.
#' @export
select_IVs <- function(data,
                       Threshold = 5e-5,
                       bfile,
                       plink_bin,
                       clump_kb = 10000,
                       clump_r2 = 0.001) {
    if (!is.data.frame(data)) stop("data must be a data.frame")
    required <- c("SNP", "b.exp.pop1", "b.exp.pop2", "se.exp.pop1",
                  "se.exp.pop2", "pval.exp.pop1")
    missing_columns <- setdiff(required, names(data))
    if (length(missing_columns)) stop("Missing columns: ", paste(missing_columns, collapse = ", "))
    if (!is.numeric(Threshold) || length(Threshold) != 1L ||
        !is.finite(Threshold) || Threshold <= 0 || Threshold > 1) {
        stop("Threshold must be in (0, 1]")
    }
    if (!is.numeric(clump_kb) || length(clump_kb) != 1L ||
        !is.finite(clump_kb) || clump_kb <= 0) stop("clump_kb must be positive")
    if (!is.numeric(clump_r2) || length(clump_r2) != 1L ||
        !is.finite(clump_r2) || clump_r2 < 0 || clump_r2 > 1) {
        stop("clump_r2 must be in [0, 1]")
    }
    dat <- as.data.frame(data)
    if (anyNA(dat$SNP) || any(!nzchar(as.character(dat$SNP))) || anyDuplicated(dat$SNP)) {
        stop("SNP identifiers must be non-missing and unique")
    }
    numeric_cols <- setdiff(required, "SNP")
    if (!all(vapply(dat[numeric_cols], is.numeric, logical(1))) ||
        any(!is.finite(as.matrix(dat[numeric_cols])))) {
        stop("Exposure effects, standard errors and P values must be finite numeric values")
    }
    if (any(dat$se.exp.pop1 <= 0 | dat$se.exp.pop2 <= 0)) {
        stop("Exposure standard errors must be positive")
    }
    if (any(dat$pval.exp.pop1 < 0 | dat$pval.exp.pop1 > 1)) {
        stop("Auxiliary-exposure P values must be in [0, 1]")
    }
    dat <- dat[dat$pval.exp.pop1 < Threshold, , drop = FALSE]
    w1 <- 1 / dat$se.exp.pop1^2
    w2 <- 1 / dat$se.exp.pop2^2
    dat$b.exp.fema <- (w1 * dat$b.exp.pop1 + w2 * dat$b.exp.pop2) / (w1 + w2)
    dat$se.exp.fema <- sqrt(1 / (w1 + w2))
    dat$pval.exp.fema <- pmax(
        2 * stats::pnorm(abs(dat$b.exp.fema / dat$se.exp.fema), lower.tail = FALSE),
        .Machine$double.xmin
    )
    dat$Threshold <- rep(Threshold, nrow(dat))
    if (!nrow(dat)) return(dat)
    if (any(!is.finite(dat$b.exp.fema)) || any(!is.finite(dat$se.exp.fema)) ||
        any(!is.finite(dat$pval.exp.fema)) || any(dat$se.exp.fema <= 0)) {
        stop("Non-finite FEMA calculation; check effect and standard-error scales")
    }
    if (missing(bfile) || !is.character(bfile) || length(bfile) != 1L ||
        is.na(bfile) || !nzchar(bfile)) stop("Supply a local target-population bfile prefix")
    if (missing(plink_bin) || !is.character(plink_bin) || length(plink_bin) != 1L ||
        is.na(plink_bin) || !nzchar(plink_bin)) stop("Supply a local plink_bin path")
    if (!all(file.exists(paste0(bfile, c(".bed", ".bim", ".fam"))))) {
        stop("Missing target-population PLINK reference files")
    }
    if (!file.exists(plink_bin)) stop("PLINK executable not found: ", plink_bin)
    if (!requireNamespace("MRAPSS", quietly = TRUE)) {
        stop("select_IVs() requires the MRAPSS package")
    }
    clumped <- MRAPSS::clump(
        dat, IV.Threshold = 1, SNP_col = "SNP", pval_col = "pval.exp.fema",
        clump_kb = clump_kb, clump_r2 = clump_r2, clump_p = 1,
        bfile = bfile, plink_bin = plink_bin
    )
    if (is.null(clumped) || !nrow(clumped)) return(dat[FALSE, , drop = FALSE])
    clumped <- as.data.frame(clumped)
    # MRAPSS::clump records its internal IV.Threshold (1); XMR needs the
    # auxiliary-exposure eligibility threshold.
    clumped$Threshold <- rep(Threshold, nrow(clumped))
    clumped
}
