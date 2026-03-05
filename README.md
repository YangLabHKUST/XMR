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
devtools::install_github("YangLabHKUST/XMR")
```

## Usage

We illustrate how to perform cross-population MR analysis using XMR with a real-data example: LDL cholesterol (LDLC, exposure) and myocardial infarction (MI, outcome), with Europeans (EUR) as the auxiliary population and East Asians (EAS) as the target population.

The XMR analysis comprises two main steps:

- **Step 1:** Prepare data and estimate background parameters (the C matrix and Omega matrix via cross-population LD score regression).
- **Step 2:** Fit XMR for causal inference.

For a step-by-step walkthrough, see the [XMR tutorial: causal effect of LDLC on MI](https://htmlpreview.github.io/?https://github.com/YangLabHKUST/XMR/blob/main/XMR_tutorial_LDLC_MI.html) ([download link](https://github.com/YangLabHKUST/XMR/blob/main/XMR_tutorial_LDLC_MI.html)).

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

We applied XMR and 15 existing summary-level MR methods to
- (1) simulations;
- (2) test the causal effects of 35 traits on 2 negative control outcomes (Skin tanning ability, Natural hair color) in Africans (AFR) and Central/South Asians (CSA);
- (3) infer causal relationships in 3 underrepresented populations: East Asians, Central/South Asians and Africans.

We provide [source codes](https://github.com/YangLabHKUST/XMR_reproduce) for replicating the simulation and real data analysis results in the XMR paper.

**Simulation:**

[Experiments and visualization](https://github.com/YangLabHKUST/XMR_reproduce/blob/main/sim/code/simulation_reproduce_plot.ipynb).

**Negative-control studies:**
[Format data](https://github.com/YangLabHKUST/XMR_reproduce/blob/main/nc/code/format_data.ipynb);
[XMR in AFR](https://github.com/YangLabHKUST/XMR_reproduce/blob/main/nc/code/run_XMR_AFR.ipynb);
[XMR in CSA](http://github.com/YangLabHKUST/XMR_reproduce/blob/main/nc/code/run_XMR_CSA.ipynb);
[Other methods in AFR](https://github.com/YangLabHKUST/XMR_reproduce/blob/main/nc/code/run_other_AFR.ipynb);
[Other methods in CSA](https://github.com/YangLabHKUST/XMR_reproduce/blob/main/nc/code/run_other_CSA.ipynb);
[Visualization](https://github.com/YangLabHKUST/XMR_reproduce/blob/main/nc/code/nc-plot.ipynb).

**Real-data analysis for EAS:**
[Format data](https://github.com/YangLabHKUST/XMR_reproduce/blob/main/real_data_EAS/code/format_data.ipynb);
[XMR in BBJ](https://github.com/YangLabHKUST/XMR_reproduce/blob/main/real_data_EAS/code/run_XMR_BBJ.ipynb); [XMR in TPMI](https://github.com/YangLabHKUST/XMR_reproduce/blob/main/real_data_EAS/code/run_XMR_TPMI.ipynb);
[ ]();
[ ]();
[ ]();
[ ]().

**Real-data analysis for CSA:**
Coming soon
[ ]();
[ ]();
[ ]();
[ ]();
[ ]();
[ ]().

**Real-data analysis for AFR:**
Coming soon
[ ]();
[ ]();
[ ]();
[ ]();
[ ]();
[ ]().


We have also released data and results related to the abbove experiments. Follow the below setup for reproduction:
### Setup

#### 1. Clone this repository
```bash
git clone https://github.com/YangLabHKUST/XMR_reproduce.git
cd XMR_reproduce
```

#### Directory structure

```
XMR_reproduce/
├── nc/                  # Negative-control analysis (AFR & CSA; coming soon)
├── real_data_CSA_AFR/   # Real data analysis in CSA & AFR
├── real_data_EAS/       # Real data analysis in EAS
└── sim/                 # Simulations
```

#### 2. Download data

Download the following archives and place them in the repository root:

| File | Size | Link |
|------|------|------|
| `nc_data.tar.gz` | ~X GB | [Google Drive](link) |
| `real_data_CSA_AFR_data.tar.gz` | ~X GB | Coming soon |
| `real_data_EAS_data.tar.gz` | ~X GB | [Google Drive](link) |
| `sim_data.tar.gz` | ~23 MB | [Zenodo link](https://doi.org/10.5281/zenodo.18872085) |


#### 3. Extract

```bash
tar xzvf nc_data.tar.gz
tar xzvf real_data_CSA_AFR_data.tar.gz
tar xzvf real_data_EAS_data.tar.gz
tar xzvf sim_data.tar.gz
```

The data files will be automatically merged into the existing code directories.

#### 4. External resources (download separately)

The following large reference files may not be included in the archives due to size. Please download them manually:

- **1000 Genomes PLINK files**: download from [link] and place in `nc/1kg_pops/`
- **PLINK software**: download from [https://www.cog-genomics.org/plink2](https://www.cog-genomics.org/plink2)


#### 5. Run the analysis
All scripts assume the **working directory is the repository root** (`XMR_reproduce/`).

```r
# In R
setwd("/path/to/XMR_reproduce")  # set to your local path
source("nc/format_data.R")
```

To run the analysis for XMR and 15 compared methods, you need to install the corresponding packages first:

```r
# In R
#install.packages("devtools") #install.packages("remotes")
devtools::install_github("YangLabHKUST/XMAP")
devtools::install_github("hhoulei/TEMR")
devtools::install_github("YangLabHKUST/MR-APSS")
remotes::install_github("MRCIEU/TwoSampleMR")
devtools::install_github("jean997/cause@v1.2.0")
install_github("tye27/mr.divw")
devtools::install_github("gqi/MRMix")
devtools::install_github("xue-hr/MRcML")
devtools::install_github("rondolab/MR-PRESSO")
install.packages("MendelianRandomization")
install.packages(“robustbase”)
```


## Reference

Xinrui Huang, Zitong Chao, Zhiwei Wang, Xianghong Hu, and Can Yang. *XMR: A cross-population Mendelian randomization method for causal inference using genome-wide summary statistics.* 2026.

## Contact

Please feel free to contact Xinrui Huang (xhuangcn@connect.ust.hk), Prof. Xianghong Hu (huxh@szu.edu.cn), or Prof. Can Yang (macyang@ust.hk) if you have any questions.
