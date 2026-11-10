# Run after preparing full-genome GWAS files, LD scores and harmonized candidate data.
# GWAS filenames and LD-score columns follow ?XMR::estimate_bgparas.
# harmonized_data must contain the columns in ?XMR::select_IVs and ?XMR::fit_XMR.
# Supply a local TARGET-population PLINK reference; no reference is downloaded.
run_xmr_workflow <- function(ldscore0, harmonized_data,
                             exposure, outcome, pop_ref, pop_target,
                             gwas_dir, bfile, plink_bin, threshold = 5e-5) {
    background <- XMR::estimate_bgparas(
        ldscore0, exposure, outcome, pop_ref, pop_target, gwas_dir,
        cross_intercept_zero = TRUE
    )
    matrices <- XMR::qc_bgparas(background, exposure, outcome, pop_ref, pop_target)
    ivs <- XMR::select_IVs(harmonized_data, Threshold = threshold,
                          bfile = bfile, plink_bin = plink_bin)
    if (!nrow(ivs)) stop("No candidate IV remains")
    tolerance <- if (nrow(ivs) < 100L) 1e-8 else 1e-7
    XMR::fit_XMR(ivs, exposure = exposure, outcome = outcome,
                 C = matrices$C, Omega0 = matrices$Omega0,
                 Threshold = threshold, tol1 = tolerance, tol2 = tolerance,
                 min_thres = 1e-2)
}
