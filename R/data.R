#' C matrix
#'
#' Background intercept matrix for the EUR-to-BBJ LDL-C to MI example.
#' Cross-population intercepts are fixed at zero during stage-2 LDSC;
#' within-population intercepts remain estimated.
#'
#' @format A 3x3 numeric matrix ordered as EUR LDL-C, BBJ LDL-C and BBJ MI.
#' @source Rebuilt by inst/tutorial/XMR_tutorial_LDLC_MI.Rmd; input checksums
#'   and settings are recorded in inst/extdata/LDLC_MI_example_provenance.json.
"C"

#' Omega matrix
#'
#' Background genetic coefficients for the EUR-to-BBJ LDL-C to MI
#' example, estimated jointly with the bundled C and processed by qc_bgparas().
#'
#' @format A 3x3 numeric matrix ordered as EUR LDL-C, BBJ LDL-C and BBJ MI.
#' @source Rebuilt by inst/tutorial/XMR_tutorial_LDLC_MI.Rmd; input checksums
#'   and settings are recorded in inst/extdata/LDLC_MI_example_provenance.json.
"Omega"

#' Example clumped data
#'
#' Harmonized LDL-C and MI summary statistics for the EUR-to-BBJ example.
#' Candidates satisfy auxiliary-exposure P < 5e-5, are ranked by signed
#' inverse-variance FEMA, and are clumped using EAS LD (10 Mb, r-squared 0.001).
#' The Threshold column contains 5e-5. FEMA effects, standard errors and P
#' values are retained together with the SNP associations and LD scores.
#'
#' @format A data frame with 289 rows and 19 variables.
#' @source Rebuilt by inst/tutorial/XMR_tutorial_LDLC_MI.Rmd; input checksums
#'   and settings are recorded in inst/extdata/LDLC_MI_example_provenance.json.
"clumped_data"
