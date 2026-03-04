#' Estimate XMAP background parameters via two-stage bivariate LDSC
#'
#' @param ldscore0 A data.frame of LD scores.
#' @param exposure_code Character. Exposure trait code, e.g. "30780".
#' @param outcome_code Character. Outcome trait code, e.g. "411.2".
#' @param pop_ref Character. Auxiliary population label, e.g. "EUR".
#' @param pop_target Character. Target population label, e.g. "BBJ".
#' @param gwas_dir Character. Directory containing formatted GWAS files.
#' @param snp_col Character. Column name for SNP rsID in ldscore0. Default "rsid".
#' @param a1_col Character. Column name for allele1 in ldscore0. Default "allele1".
#' @param z2_threshold Numeric. Z-squared threshold for filtering in stage 1. Default 30.
#'
#' @return A data.frame with 6 rows, one per trait-population pair,
#'   containing columns: trait1, trait2, pop1, pop2, omega11, omega12, omega22, c1, c12, c2.
#'
#' @export
estimate_bgparas <- function(ldscore0,
                             exposure_code,
                             outcome_code,
                             pop_ref = "EUR",
                             pop_target = "BBJ",
                             gwas_dir = "./formatted_data",
                             snp_col = "rsid",
                             a1_col = "allele1",
                             z2_threshold = 30) {

    # Build the 6 trait-population pairs:
    # all unique (unordered) pairs from {pop_ref_exposure, pop_target_exposure, pop_target_outcome}
    traits <- c(paste0(pop_ref,    "_", exposure_code),
                paste0(pop_target, "_", exposure_code),
                paste0(pop_target, "_", outcome_code))

    pairs <- data.frame(
        trait1 = c(traits[1], traits[1], traits[1], traits[2], traits[2], traits[3]),
        trait2 = c(traits[1], traits[2], traits[3], traits[2], traits[3], traits[3]),
        pop1   = c(pop_ref,   pop_ref,   pop_ref,   pop_target, pop_target, pop_target),
        pop2   = c(pop_ref,   pop_target, pop_target, pop_target, pop_target, pop_target),
        stringsAsFactors = FALSE
    )

    results_list <- vector("list", nrow(pairs))

    for (i in seq_len(nrow(pairs))) {

        t1 <- pairs$trait1[i];  t2 <- pairs$trait2[i]
        p1 <- pairs$pop1[i];    p2 <- pairs$pop2[i]

        message(sprintf("Pair %d/%d: %s vs %s (%s x %s)",
                        i, nrow(pairs), t1, t2, p1, p2))

        result <- try({

            dat1 <- data.table::fread(file.path(gwas_dir, t1))
            dat2 <- data.table::fread(file.path(gwas_dir, t2))

            # Overlap SNPs across LD scores and both GWAS datasets
            snps <- Reduce(intersect, list(ldscore0[[snp_col]], dat1$SNP, dat2$SNP))
            dat1_ldsc <- dat1[match(snps, dat1$SNP), ]
            dat2_ldsc <- dat2[match(snps, dat2$SNP), ]
            ldscore   <- ldscore0[match(snps, ldscore0[[snp_col]]), ]

            # Flip alleles to harmonize SNP effect directions
            z1 <- dat1_ldsc$Z
            z2 <- dat2_ldsc$Z

            flip1 <- which(dat1_ldsc$A1 != ldscore[[a1_col]] &
                           dat1_ldsc$A1 != comple(ldscore[[a1_col]]))
            z1[flip1] <- -z1[flip1]

            flip2 <- which(dat2_ldsc$A1 != ldscore[[a1_col]] &
                           dat2_ldsc$A1 != comple(ldscore[[a1_col]]))
            z2[flip2] <- -z2[flip2]

            # Filter extreme Z-scores for stage 1
            idx1 <- which(z1^2 < z2_threshold & z2^2 < z2_threshold)

            # LD-score-based regression weights
            w1 <- 1 / sapply(ldscore[, p1], function(x) max(x, 1))
            w2 <- 1 / sapply(ldscore[, p2], function(x) max(x, 1))
            cross_col <- paste0(p1, "_", p2)

            # Stage 1: estimate intercepts (unconstrained)
            fit1 <- XMAP::estimate_gc(
                data.frame(Z = z1[idx1], N = dat1_ldsc$N[idx1]),
                data.frame(Z = z2[idx1], N = dat2_ldsc$N[idx1]),
                ldscore[idx1, p1], ldscore[idx1, p2], ldscore[idx1, cross_col],
                reg_w1 = w1[idx1], reg_w2 = w2[idx1],
                reg_wx = sqrt(w1[idx1] * w2[idx1]),
                constrain_intercept = FALSE
            )

            # Stage 2: fix intercepts, estimate slopes
            fit2 <- XMAP::estimate_gc(
                data.frame(Z = z1, N = dat1_ldsc$N),
                data.frame(Z = z2, N = dat2_ldsc$N),
                ldscore[, p1], ldscore[, p2], ldscore[, cross_col],
                reg_w1 = w1, reg_w2 = w2,
                reg_wx = sqrt(w1 * w2),
                constrain_intercept = TRUE,
                fit1$tau1$coefs[1], fit1$tau2$coefs[1], fit1$theta$coefs[1]
            )

            data.frame(trait1 = t1, trait2 = t2, pop1 = p1, pop2 = p2,
                       omega11 = fit2$tau1$coefs[2],
                       omega12 = fit2$theta$coefs[2],
                       omega22 = fit2$tau2$coefs[2],
                       c1  = fit2$tau1$coefs[1],
                       c12 = fit2$theta$coefs[1],
                       c2  = fit2$tau2$coefs[1])

        }, silent = TRUE)

        if (inherits(result, "try-error")) {
            warning(sprintf("Pair %d failed, skipping. Error: %s", i, conditionMessage(attr(result, "condition"))))
            results_list[[i]] <- NULL
        } else {
            results_list[[i]] <- result
        }
    }

    do.call(rbind, results_list)
}