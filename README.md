# XMR

**XMR** (Cross-Population Mendelian Randomization) is a probabilistic method for estimating causal effects between an exposure and an outcome using genome-wide summary statistics from multiple populations.

XMR improves the power and robustness of causal inference in underrepresented (small-sample) populations by leveraging information from a large-sample auxiliary population. 
Specifically, XMR decomposes the observed SNP–trait effects into true causal effects and confounding factors (e.g., pleiotropy, population structure) hidden in summary statistics. 
By explicitly modelling the genetic correlation between two populations, XMR effectively borrows strength from the large-sample group. 
XMR accounts for auxiliary-exposure IV selection in its likelihood.

<img width="2314" height="1484" alt="649dc9cb-1" src="https://github.com/user-attachments/assets/1c59e233-90c5-4ffc-b6dc-a6c8db2b2335" />

**Overview of the XMR method.** XMR estimates the causal effect $\beta$ between exposure $X\_2$ and outcome $Y\_2$ in a small-sample population by leveraging data on the same exposure $X\_1$ from a large-sample population. The method involves several key elements:
**(A)** IVs are selected from the large-sample population ($X\_1$) to improve power compared to the limited IVs available from the small-sample population ($X\_2$). The distributions of observed $-\log\_{10}(p)$ values for SNP–exposure associations across chromosomes are shown.
**(B)** The XMR model diagram. Arrowed lines represent directed effects. The blue dashed line indicates the correlation between $X\_1$ and $X\_2$.
**(C)** Selection bias and confounding factors contribute to the observed SNP–trait associations.
**(D)** An illustrative example of causal inference between SHBG (sex hormone-binding globulin) and T2D (type 2 diabetes) in an African population, using conventional two-sample MR methods (left) and XMR (right). The estimated causal effect is shown as a red line, with the 95% confidence interval shaded in transparent red. Triangles represent observed SNP effect sizes ($\hat{\gamma}\_{2,j}$ and $\hat{\Gamma}\_{2,j}$), colored by their posterior probability of foreground membership (IV contribution) ($Z\_j = 1$ in dark blue; $Z\_j = 0$ in light blue).

## Installation

```r
# install.packages("devtools")
devtools::install_github("YangLabHKUST/XMR")
```

Typical installation time: < 5 minutes. Required time on a standard desktop computer should be comparable.

## Usage

We illustrate how to perform cross-population MR analysis using XMR with a real-data example: LDL cholesterol (LDLC, exposure) and myocardial infarction (MI, outcome), with Europeans (EUR) as the auxiliary population and East Asians (EAS) as the target population.

The XMR analysis comprises two main steps:

- **Step 1:** Prepare data and estimate background parameters (the C matrix and Omega matrix via cross-population LD score regression).
- **Step 2:** Fit XMR for causal inference.

The [tutorial: causal effect of LDLC on MI](https://htmlpreview.github.io/?https://github.com/YangLabHKUST/XMR/blob/main/XMR_tutorial_LDLC_MI.html) ([download link](https://github.com/YangLabHKUST/XMR/blob/main/XMR_tutorial_LDLC_MI.html)) walks through data preparation, IV selection and model fitting, with code and example results. Its executable source is [`inst/tutorial/XMR_tutorial_LDLC_MI.Rmd`](inst/tutorial/XMR_tutorial_LDLC_MI.Rmd). A compact wrapper is available in [`inst/examples/workflow.R`](inst/examples/workflow.R).

For a quick start, you can skip Step 1 and proceed directly to Step 2 using the example data we have prepared.

```r
library(XMR)

exposure <- "LDLC"
outcome  <- "MI"

# Auxiliary-exposure IV selection threshold
threshold <- 5e-05

# Load example data
data(C)
data(Omega)
data(clumped_data) # after IV selection and LD clumping

# Fit XMR
XMR_res <- fit_XMR(
  data = clumped_data,
  C = C,
  Omega0 = Omega,
  Threshold = threshold,
  tol1 = 1e-07,
  tol2 = 1e-07,
  min_thres = 1e-2
)
```

### Data preparation

Use harmonized effects with aligned alleles and compatible exposure units in both populations. Estimate background parameters from genome-wide summary statistics, not only selected IVs. The default `estimate_bgparas(..., cross_intercept_zero = TRUE)` fixes cross-population intercepts to zero during the second LDSC stage while retaining the estimated marginal and within-population intercepts. Set `cross_intercept_zero = FALSE` to estimate cross-population intercepts without fixing them to zero.

```r
# ldscore0 and formatted GWAS files follow ?estimate_bgparas
bg <- estimate_bgparas(ldscore0, exposure, outcome,
                      pop_ref = "EUR", pop_target = "BBJ",
                      gwas_dir = "formatted_data")
matrices <- qc_bgparas(bg, exposure, outcome, "EUR", "BBJ")

# harmonized_data contains aligned effects/SEs, auxiliary P values and LD scores.
# bfile must be a TARGET-population PLINK reference prefix.
ivs <- select_IVs(harmonized_data, Threshold = threshold,
                 bfile = "reference/target", plink_bin = "/path/to/plink")
fit <- fit_XMR(ivs, C = matrices$C, Omega0 = matrices$Omega0,
              Threshold = threshold)
```

`select_IVs()` first selects `pval.exp.pop1 < Threshold`, then computes signed inverse-variance FEMA P values from the two exposure GWASs. FEMA P values prioritize eligible candidates during target-population LD clumping. It uses `MRAPSS::clump()` with a 10-Mb window and `r² = 0.001` by default, requiring a local target-population reference and PLINK executable. Use the auxiliary-exposure selection threshold as `Threshold` in `fit_XMR()`.

### Input data format

The input `data` should be a data.frame containing the following columns:

| Column | Description |
|--------|-------------|
| `b.exp.pop1` | SNP–exposure effect in the auxiliary population (pop1) |
| `b.exp.pop2` | SNP–exposure effect in the target population (pop2) |
| `b.out.pop2` | SNP–outcome effect in the target population (pop2) |
| `se.exp.pop1` | Standard error of `b.exp.pop1` |
| `se.exp.pop2` | Standard error of `b.exp.pop2` |
| `se.out.pop2` | Standard error of `b.out.pop2` |
| `L2.pop1` | LD score in the auxiliary population; defaults to all ones if not provided (i.e., no LD correction) |
| `L2.pop2` | LD score in the target population; defaults to all ones if not provided (i.e., no LD correction) |
| `L12` | Cross-population LD score between pop1 and pop2; defaults to all ones if not provided (i.e., no LD correction) |

### Parameters related to confounding factors

- **C matrix**: A 3×3 matrix capturing the effects of sample structure (population stratification, cryptic relatedness, sample overlap, etc.).
- **Omega matrix**: A 3×3 variance–covariance matrix of polygenic effects.

Both can be estimated using bivariate LD score regression.

### Expected output

`XMR_res` is a list containing:
- `exposure`: phecode of the exposure
- `outcome`: phecode of the outcome
- `beta`: estimated causal effect of the exposure on the outcome
- `beta.se`: standard error of the causal effect estimate
- `beta.pvalue`: p-value for the causal effect
- `tau.sq`: estimated variance of the pleiotropic effect
- `SigmaX`: estimated SNP–trait covariance matrix
- `post`: posterior moments; `post$Pi` contains foreground-membership probabilities (support for IV contribution)
- `nIV`: number of IVs used
- `nvalid`: model-estimated effective IV counts (`pi0 * nIV`)
- `pi0`: estimated foreground-mixture proportion among selected candidates
- `Threshold`: IV selection threshold
- `fit1_elbos`, `fit2_elbos`: evidence lower bounds in the two stages
- `fit1_likelis`, `fit2_likelis`: likelihoods in the two stages

The bundled LDLC → MI example contains matching background matrices and clumped IV data.

For the bundled example, the result is `nIV = 289`, `beta = 0.183`, `SE = 0.044908` and `P = 4.5902e-05`.
Input checksums, settings and software versions are recorded in [`inst/extdata/LDLC_MI_example_provenance.json`](inst/extdata/LDLC_MI_example_provenance.json). Runtime depends on the number of IVs and convergence.

## Reproducibility

Experiment scripts, data downloads and instructions are maintained in the companion repository [YangLabHKUST/XMR_reproduce](https://github.com/YangLabHKUST/XMR_reproduce).

GWAS data sources are listed in the following tables:

| Analysis | Data sources |
|----------|--------------|
| Negative-control studies | [nc_data_source.csv](nc_data_source.csv) |
| East Asian populations | [real_data_EAS_data_source.csv](real_data_EAS_data_source.csv) |
| African populations | [real_data_AFR_data_source.csv](real_data_AFR_data_source.csv) |
| Central/South Asian populations | [real_data_CSA_data_source.csv](real_data_CSA_data_source.csv) |

## System Requirements

- **Software dependencies**: R (>= 4.2.0); R packages: Matrix, expm, stats, mvtnorm, dplyr, reshape2, glmnet, MASS, data.table, readr, magrittr, doParallel, XMAP, MRAPSS
- **Tested on**: Intel® Xeon® CPU E5-2699 v4 @ 2.20 GHz, R 4.2.2
- **Hardware**: No non-standard hardware required. A standard desktop computer is sufficient.
  

## Reference

Xinrui Huang, Zitong Chao, Zhiwei Wang, Xianghong Hu, and Can Yang. *XMR: A cross-population Mendelian randomization method for causal inference using genome-wide summary statistics.* medRxiv, 2026. DOI: [10.64898/2026.03.10.26348003](https://www.medrxiv.org/content/10.64898/2026.03.10.26348003v1)


## Contact

Please feel free to contact Xinrui Huang (xhuangcn@connect.ust.hk), Prof. Xianghong Hu (huxh@szu.edu.cn), or Prof. Can Yang (macyang@ust.hk) if you have any questions.
