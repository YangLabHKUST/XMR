# XMR 0.1.0.9000

* `estimate_bgparas()` fixes cross-population LDSC intercepts to zero during
  stage 2 by default. Set `cross_intercept_zero = FALSE` to estimate these
  intercepts without fixing them to zero.
* `select_IVs()` selects candidates using the auxiliary-exposure P value,
  ranks them by signed inverse-variance FEMA, and clumps them using
  target-population LD (10 Mb, r-squared 0.001 by default).
* `qc_bgparas()` accepts either `omega12` or `Omega12` as the column name.
* `fit_XMR()` uses the auxiliary-exposure selection threshold. Posterior
  foreground-membership probabilities quantify support for IV contribution;
  `nvalid` reports model-estimated effective IV counts.
* The package includes LDLC-MI background matrices, clumped IV data, a
  five-step HTML tutorial and its R Markdown source.
