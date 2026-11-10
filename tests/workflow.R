library(XMR)

must_fail <- function(expr) stopifnot(inherits(tryCatch({ force(expr); NULL }, error = identity), "error"))

# Work on independent closures: do not modify installed package namespaces.
mock_calls <- list()
env <- new.env(parent = environment(select_IVs))
env$`::` <- function(pkg, name) {
    if (as.character(substitute(pkg)) == "MRAPSS" && as.character(substitute(name)) == "clump") {
        return(function(dat, ...) {
            mock_calls[[length(mock_calls) + 1L]] <<- list(data = dat, args = list(...))
            dat$Threshold <- list(...)$IV.Threshold
            dat[order(dat$pval.exp.fema), , drop = FALSE]
        })
    }
    getExportedValue(as.character(substitute(pkg)), as.character(substitute(name)))
}
env$requireNamespace <- function(...) TRUE
selector <- select_IVs
environment(selector) <- env

tmp <- tempfile("xmr-selection-"); dir.create(tmp)
ref <- file.path(tmp, "target")
plink <- file.path(tmp, "plink")
invisible(file.create(c(paste0(ref, c(".bed", ".bim", ".fam")), plink)))
dat <- data.frame(SNP = paste0("rs", 1:4),
                  b.exp.pop1 = c(1, 1, 1, 1), b.exp.pop2 = c(-0.9, 1, 1, 1),
                  se.exp.pop1 = 0.1, se.exp.pop2 = 0.1,
                  pval.exp.pop1 = c(1e-8, 1e-8, 5e-5, 0.5), b.out.pop2 = 0.2)
original <- dat
ivs <- selector(dat, bfile = ref, plink_bin = plink)
stopifnot(identical(dat, original), setequal(ivs$SNP, c("rs1", "rs2")))
# Opposite signs cancel instead of being meta-analysed as unsigned Z statistics.
stopifnot(ivs$pval.exp.fema[ivs$SNP == "rs1"] > 0.4,
          ivs$SNP[1] == "rs2", all(ivs$Threshold == 5e-5))
args <- mock_calls[[1]]$args
stopifnot(args$IV.Threshold == 1, args$clump_p == 1,
          args$pval_col == "pval.exp.fema", args$bfile == ref,
          args$clump_kb == 10000, args$clump_r2 == 0.001)
stopifnot(nrow(selector(dat[FALSE, ])) == 0,
          nrow(selector(transform(dat, pval.exp.pop1 = 1))) == 0)
must_fail(selector(dat, Threshold = 0))
must_fail(selector(transform(dat, se.exp.pop1 = 0)))
must_fail(selector(transform(dat, SNP = "duplicate")))
must_fail(selector(transform(dat, b.exp.pop1 = NA_real_)))
must_fail(selector(dat, bfile = "missing-reference", plink_bin = plink))
unlink(tmp, recursive = TRUE)

# Background stage 2 fixes only cross-population intercepts and supports estimated cross-population intercepts.
calls <- list()
env <- new.env(parent = environment(estimate_bgparas))
gwas <- data.frame(SNP = paste0("rs", 1:4), A1 = "A", Z = 1:4 / 2, N = 1000)
env$`::` <- function(pkg, name) {
    package <- as.character(substitute(pkg)); fun <- as.character(substitute(name))
    if (package == "data.table" && fun == "fread") return(function(...) gwas)
    if (package == "XMAP" && fun == "estimate_gc") return(function(...) {
        a <- list(...); calls[[length(calls) + 1L]] <<- a
        if (isTRUE(a$constrain_intercept)) {
            list(tau1 = list(coefs = c(a$int1, 0.001)),
                 tau2 = list(coefs = c(a$int2, 0.001)),
                 theta = list(coefs = c(a$intx, 0.0001)))
        } else {
            list(tau1 = list(coefs = c(1.1, 0.001)),
                 tau2 = list(coefs = c(1.2, 0.001)),
                 theta = list(coefs = c(0.2, 0.0001)))
        }
    })
    getExportedValue(package, fun)
}
estimator <- estimate_bgparas; environment(estimator) <- env
ld <- data.frame(rsid = gwas$SNP, allele1 = "A", EUR = 1:4, BBJ = 2:5,
                 EUR_EUR = 1:4, EUR_BBJ = 1:4, BBJ_BBJ = 2:5)
bg <- estimator(ld, "X", "Y")
cross <- bg$pop1 != bg$pop2
stopifnot(nrow(bg) == 6L, all(bg$c12[cross] == 0), all(bg$c12[!cross] == 0.2),
          all(bg$c1 == 1.1), all(bg$c2 == 1.2), length(calls) == 12L)
bg_old <- estimator(ld, "X", "Y", cross_intercept_zero = FALSE)
stopifnot(all(bg_old$c12 == 0.2))
must_fail(estimator(ld, "X", "Y", cross_intercept_zero = NA))

# estimate_bgparas output can be consumed directly; uppercase Omega12 is also accepted.
bg$omega12 <- c(0.01, 0.001, 0.001, 0.01, 0.001, 0.01)
bg$c12 <- c(1.1, 0, 0, 1.1, 0.1, 1.1)
uppercase_input <- bg; names(uppercase_input)[names(uppercase_input) == "omega12"] <- "Omega12"
stopifnot(identical(qc_bgparas(bg, "X", "Y", "EUR", "BBJ"),
                    qc_bgparas(uppercase_input, "X", "Y", "EUR", "BBJ")))
cat("PASS: IV selection, intercept constraints and accepted input column names\n")

# Bundled LDLC-MI data must match the XMR workflow.
example_data <- new.env()
data(list = c("C", "Omega", "clumped_data"), package = "XMR", envir = example_data)
stopifnot(all(example_data$C[1, 2:3] == 0),
          all(example_data$C[2:3, 1] == 0),
          identical(dimnames(example_data$C), dimnames(example_data$Omega)),
          !anyDuplicated(example_data$clumped_data$SNP),
          all(example_data$clumped_data$Threshold == 5e-5),
          all(example_data$clumped_data$pval.exp.pop1 < 5e-5),
          all(c("b.exp.fema", "se.exp.fema", "pval.exp.fema") %in%
              names(example_data$clumped_data)))
cat("PASS: bundled example uses matched C0 backgrounds and IV selection\n")
