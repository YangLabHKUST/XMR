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

We applied XMR and 15 existing summary-level MR methods across three key domains:
- **Simulations**: evaluating method performance under various scenarios.
- **Negative-control studies**: testing the causal effects of 35 traits on 2 negative-control outcomes (skin tanning ability, natural hair color) in Africans (AFR) and Central/South Asians (CSA).
- **Real-data analysis**: inferring causal relationships in 3 underrepresented populations — East Asians (EAS), Central/South Asians (CSA), and Africans (AFR).

Source code and data for reproducing all results are available at [YangLabHKUST/XMR_reproduce](https://github.com/YangLabHKUST/XMR_reproduce). 
The XMR execution scripts provided below feature a parallelized framework designed to efficiently analyze multiple trait pairs simultaneously.

**Simulations:**

[Experiments and visualization](https://github.com/YangLabHKUST/XMR_reproduce/blob/main/sim/code/simulation_reproduce_plot.ipynb)

**Negative-control studies:**

[Format data](https://github.com/YangLabHKUST/XMR_reproduce/blob/main/nc/code/format_data.ipynb) |
[XMR in AFR](https://github.com/YangLabHKUST/XMR_reproduce/blob/main/nc/code/run_XMR_AFR.ipynb) |
[XMR in CSA](https://github.com/YangLabHKUST/XMR_reproduce/blob/main/nc/code/run_XMR_CSA.ipynb) |
[Other methods in AFR](https://github.com/YangLabHKUST/XMR_reproduce/blob/main/nc/code/run_other_AFR.ipynb) |
[Other methods in CSA](https://github.com/YangLabHKUST/XMR_reproduce/blob/main/nc/code/run_other_CSA.ipynb) |
[Visualization](https://github.com/YangLabHKUST/XMR_reproduce/blob/main/nc/code/nc-plot.ipynb)

**Real-data analysis for EAS:**

[Format data](https://github.com/YangLabHKUST/XMR_reproduce/blob/main/real_data_EAS/code/format_data.ipynb) |
[XMR in BBJ](https://github.com/YangLabHKUST/XMR_reproduce/blob/main/real_data_EAS/code/run_XMR_BBJ.ipynb) |
[XMR in TPMI](https://github.com/YangLabHKUST/XMR_reproduce/blob/main/real_data_EAS/code/run_XMR_TPMI.ipynb) |
[Other methods in BBJ](https://github.com/YangLabHKUST/XMR_reproduce/blob/main/real_data_EAS/code/run_other_BBJ.ipynb) |
[Other methods in TPMI](https://github.com/YangLabHKUST/XMR_reproduce/blob/main/real_data_EAS/code/run_other_TPMI.ipynb) |
[Visualization](https://github.com/YangLabHKUST/XMR_reproduce/blob/main/real_data_EAS/code/EAS-plot.ipynb)

**Real-data analysis for CSA and AFR:** 

[Format data](https://github.com/YangLabHKUST/XMR_reproduce/blob/main/real_data_CSA_AFR/code/format_data.ipynb) |
[XMR in CSA](https://github.com/YangLabHKUST/XMR_reproduce/blob/main/real_data_CSA_AFR/code/CSA/run_XMR_CSA.ipynb) |
[Other methods in CSA](http://github.com/YangLabHKUST/XMR_reproduce/blob/main/real_data_CSA_AFR/code/CSA/run_other_CSA.ipynb) |
[Visualization in CSA](https://github.com/YangLabHKUST/XMR_reproduce/blob/main/real_data_CSA_AFR/code/CSA-plot.ipynb)

[XMR in AFR](https://github.com/YangLabHKUST/XMR_reproduce/blob/main/real_data_CSA_AFR/code/AFR/run_XMR_AFR.ipynb) |
[Other methods in AFR](https://github.com/YangLabHKUST/XMR_reproduce/blob/main/real_data_CSA_AFR/code/AFR/run_other_AFR.ipynb) |
[Visualization in AFR](https://github.com/YangLabHKUST/XMR_reproduce/blob/main/real_data_CSA_AFR/code/AFR-plot.ipynb)

[Case study in EAS: SHBG to T2D](https://github.com/YangLabHKUST/XMR_reproduce/blob/main/real_data_CSA_AFR/code/case_study_in_EAS.ipynb) |
[Visualization: case study across 4 populations](https://github.com/YangLabHKUST/XMR_reproduce/blob/main/real_data_CSA_AFR/code/case_study_plot.ipynb)

### Setup

All data and results needed to reproduce the above experiments are publicly available.
See **Step 2** for download links.
Follow below steps for reproduction:

#### 1. Clone this repository

```bash
git clone https://github.com/YangLabHKUST/XMR_reproduce.git
cd XMR_reproduce
```

#### Directory structure

```
XMR_reproduce/
├── nc/                  # Negative-control analysis in AFR & CSA
├── real_data_CSA_AFR/   # Real data analysis in CSA & AFR (coming soon)
├── real_data_EAS/       # Real data analysis in EAS
└── sim/                 # Simulations
```

#### 2. Download data

We provide archived files containing formatted data, LD score files, analysis results, and other files needed for reproduction.

**Raw GWAS summary statistics** are not included due to their large size (~8–10 GB each).
Data sources are listed in the following tables — download the raw files, find the target folder in the above directory, place them in the corresponding `raw_data/` folder, and run `format_data.ipynb` in the target folder to format:

| Experiment | Data source table |
|------------|-------------------|
| Negative-control studies | [`nc_data_source.csv`](https://github.com/YangLabHKUST/XMR/blob/main/nc_data_source.csv) |
| Real-data analysis (EAS) | [`real_data_EAS_data_source.csv`](https://github.com/YangLabHKUST/XMR/blob/main/real_data_EAS_data_source.csv) |
| Real-data analysis (AFR) | [`real_data_AFR_data_source.csv`](https://github.com/YangLabHKUST/XMR/blob/main/real_data_AFR_data_source.csv) |
| Real-data analysis (CSA) | [`real_data_CSA_data_source.csv`](https://github.com/YangLabHKUST/XMR/blob/main/real_data_CSA_data_source.csv) |

Alternatively, you can **skip the raw data step** and start directly from our pre-formatted data by downloading the archives below. Then place them in the repository root and extract:

| File | Size | Link |
|------|------|------|
| `sim_data.tar.gz` | ~28 MB | [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18872086.svg)](https://doi.org/10.5281/zenodo.18872086)|
| `nc_data.tar.gz` | ~6.8 GB | [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18872963.svg)](https://doi.org/10.5281/zenodo.18872963) |
| `real_data_EAS.tar.gz` | ~5.8 GB | [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18874873.svg)](https://doi.org/10.5281/zenodo.18874873) |
| `real_data_CSA_AFR.tar.gz` | ~X GB | Coming soon |

```bash
tar xzvf sim_data.tar.gz
tar xzvf nc_data.tar.gz
tar xzvf real_data_EAS.tar.gz
tar xzvf real_data_CSA_AFR.tar.gz
```

Each archive preserves the directory structure and will merge into existing directories automatically.

#### 3. External resources (download separately)

The following large reference files may not be included in the archives due to size. Please download them manually:

- **1000 Genomes PLINK files**: download and place in `nc/1kg_pops/`; refer to [prepare_1kg_reference.sh](https://github.com/YangLabHKUST/XMR_reproduce/blob/main/prepare_1kg_reference.sh)
- **PLINK software**: download from [https://www.cog-genomics.org/plink2](https://www.cog-genomics.org/plink2)


#### 4. Run the analysis
All scripts assume the **working directory is the repository root** (`XMR_reproduce/`).

```r
# In R
setwd("/path/to/XMR_reproduce")  # set to your local path
source("nc/code/run_XMR_AFR.R")
```

To run XMR and the 15 compared methods, install the required R packages first:

```r
# In R
#install.packages("devtools") #install.packages("remotes")
devtools::install_github("YangLabHKUST/XMAP")
devtools::install_github("hhoulei/TEMR")
devtools::install_github("YangLabHKUST/MR-APSS")
remotes::install_github("MRCIEU/TwoSampleMR")
devtools::install_github("jean997/cause@v1.2.0")
devtools::install_github("tye27/mr.divw")
devtools::install_github("gqi/MRMix")
devtools::install_github("xue-hr/MRcML")
devtools::install_github("rondolab/MR-PRESSO")
install.packages("MendelianRandomization")
install.packages("robustbase")
```


## Reference

Xinrui Huang, Zitong Chao, Zhiwei Wang, Xianghong Hu, and Can Yang. *XMR: A cross-population Mendelian randomization method for causal inference using genome-wide summary statistics.* 2026.

## Contact

Please feel free to contact Xinrui Huang (xhuangcn@connect.ust.hk), Prof. Xianghong Hu (huxh@szu.edu.cn), or Prof. Can Yang (macyang@ust.hk) if you have any questions.
