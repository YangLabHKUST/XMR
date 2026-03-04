# XMR

**XMR** (Cross-Population Mendelian Randomization) is a probabilistic method for estimating causal effects between an exposure and an outcome using genome-wide summary statistics from multiple populations.

XMR improves the power and robustness of causal inference in underrepresented (small-sample) populations by leveraging information from a large-sample auxiliary population. 
Specifically, XMR decomposes the observed SNP–trait effects into true causal effects and confounding factors (e.g., pleiotropy, population structure) hidden in summary statistics. 
By explicitly modelling the genetic correlation between two populations, XMR effectively borrows strength from the large-sample group. 
XMR further corrects bias introduced by IV selection and LD clumping to reduce false positive rates.

<img width="2314" height="1484" alt="649dc9cb-1" src="https://github.com/user-attachments/assets/1c59e233-90c5-4ffc-b6dc-a6c8db2b2335" />

**Overview of the XMR method.** XMR estimates the causal effect $\beta$ between exposure $X\_2$ and outcome $Y\_2$ in a small-sample population by leveraging data on the same exposure $X\_1$ from a large-sample population. The method involves several key elements:
**(A)** IVs are selected from the large-sample population ($X\_1$) to improve power compared to the limited IVs available from the small-sample population ($X\_2$). The distributions of observed $-\log\_{10}(p)$ values for SNP–exposure associations across chromosomes are shown.
**(B)** The XMR model diagram. Arrowed lines represent directed effects. The blue dashed line indicates the correlation between $X\_1$ and $X\_2$.
**(C)** Selection bias and confounding factors contribute to the observed SNP–trait associations.
**(D)** An illustrative example of causal inference between SHBG (sex hormone-binding globulin) and T2D (type 2 diabetes) in an African population, using conventional two-sample MR methods (left) and XMR (right). The estimated causal effect is shown as a red line, with the 95% confidence interval shaded in transparent red. Triangles represent observed SNP effect sizes ($\hat{\gamma}\_{2,j}$ and $\hat{\Gamma}\_{2,j}$), colored by their posterior probability of IV validity ($Z\_j = 1$ in dark blue; $Z\_j = 0$ in light blue).

## Installation

```r
# install.packages("devtools")
devtools::install_github("yudan121/XMR")
```

## Usage

We illustrate how to perform cross-population MR analysis using XMR with a real-data example: LDL cholesterol (LDLC, exposure) and myocardial infarction (MI, outcome), with EUR as the auxiliary population and EAS as the target population.

The XMR analysis comprises two main steps:

- **Step 1:** Prepare data and estimate background parameters (the C matrix and Omega matrix via cross-population LD score regression).
- **Step 2:** Fit XMR for causal inference.

For a quick start, you can skip Step 1 and proceed directly to Step 2 using the example data we have prepared.

```r
library(XMR)

exposure <- "LDLC"
outcome  <- "MI"

# Sample sizes
N1 <- 343621  # EUR (auxiliary population)
N2 <- 72866   # EAS (target population)

# Modified IV selection threshold for correction of selection bias
threshold <- 5e-05 # IV selection threshold
t0 <- abs(qnorm(threshold / 2))
dt <- 0.13 / (sqrt(N2 / N1))
modified_threshold <- 2 * (1 - pnorm(abs(t0 + dt)))

# Load example data
data(C)
data(Omega)
data(clumped_data) # after IV selection and LD clumping

# Fit XMR
XMR_res <- fit_XMR(
  data = clumped_data,
  C = C,
  Omega0 = Omega,
  Threshold = modified_threshold,
  tol1 = 1e-07,
  tol2 = 1e-07,
  min_thres = 1e-2
)
```

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

## Reproducibility



## Reference

Xinrui Huang, Zitong Chao, Zhiwei Wang, Xianghong Hu, and Can Yang. *XMR: A cross-population Mendelian randomization method for causal inference using genome-wide summary statistics.* 2026.

## Contact

Please feel free to contact Xinrui Huang (xhuangcn@connect.ust.hk), Prof. Xianghong Hu (huxh@szu.edu.cn), or Prof. Can Yang (macyang@ust.hk) if you have any questions.
