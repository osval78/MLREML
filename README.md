# MLREML <img src="https://img.shields.io/badge/version-1.0.0-blue" alt="version 1.0.0"> <img src="https://img.shields.io/badge/license-GPL--3-green" alt="GPL-3"> <img src="https://img.shields.io/badge/R-%E2%89%A5%204.0.0-276DC3" alt="R >= 4.0.0">

**Linear mixed models by Average-Information REML for plant and animal breeding**

MLREML is an open-source R package that fits linear mixed models by restricted maximum likelihood (REML) with the Average Information (AI) algorithm. Its syntax is inspired by ASReml: every random term and the residual are written as **direct (Kronecker) products of variance structures**. For example, `fa(env, 2):vm(gen, G)` declares genomic effects whose covariance among environments follows a two-factor factor analytic model. MLREML is independent software: it contains no ASReml code, requires no licence and is distributed under GPL-3.

The package is meant for breeders, quantitative geneticists, statisticians, lecturers and students who need ASReml-type analyses (variance components, heritabilities, BLUP/GEBV, multi-environment trials, spatial analysis, pedigree and genomic models) without the cost of a commercial licence.

> **Guides in four languages** · **Guías en cuatro idiomas** · **Guias em quatro idiomas** · **四种语言的用户指南**
> English · Español · Português · 中文 (see [Documentation](#documentation)).

---

## Contents of this repository

| File | Description |
|---|---|
| [`MLREML_1.0.0.tar.gz`](MLREML_1.0.0.tar.gz) | Source package of MLREML 1.0.0 (R and C++ code, help pages, automatic tests, citation information and the English and Spanish user guides). Install it as explained in [Installation](#installation). |
| [`Paper_data_and_code.zip`](Paper_data_and_code.zip) | Data and R code needed to reproduce the article. After unzipping, the folder `Paper_data_and_code` contains `code/` (scripts `00_setup.R` to `07_example6_timing.R`) and `data/` (genomic relationship matrices `G_<year>.rds`, environmental kernels `W_<year>.rds` and phenotypes). |
| [User guide – **English** (PDF)](MLREML_documentation_Espa%C3%B1ol_v1.0.0.pdf) | MLREML user guide, version 1.0.0, in **English**. |
| [Guía del usuario – **español** (PDF)](MLREML_documentacion_ENGLISH_v1.0.0.pdf) | Guía del usuario de MLREML, versión 1.0.0, en **español**. |
| [Guia do usuário – **português** (PDF)](MLREML_documentacao_Portugues_v1.0.0.pdf) | Guia do usuário do MLREML, versão 1.0.0, em **português**. |
| [用户指南 – **中文** (PDF)](MLREML_documentation_Chino_v1.0.0.pdf) | MLREML 用户指南，版本 1.0.0（**中文**）。 |

Each guide (about 100 pages) contains the theory, the model syntax, nine complete examples with public data, advanced use and troubleshooting, a simulation study, the numerical validation of the package, and the reference manual of every function.

**Inside the source package (`MLREML_1.0.0.tar.gz`)**

```
MLREML/
├── DESCRIPTION, NAMESPACE, NEWS.md (DESCRIPTION lists the six authors)
├── R/            mlreml.R (main function), structures.R (variance structures), sparse.R (sparse engine),
│                 kr.R (Kenward–Roger), methods.R (summary, coef, predict, wald, ...),
│                 relationships.R (A, A⁻¹, G, H, QC), gblup.R (spectral GBLUP and cross-validation)
├── src/          sparse_tools.cpp (sparse inverse and inbreeding coefficients in C++, via Rcpp)
├── man/          help pages of all exported functions
├── tests/        automatic tests (testthat): agreement with lme4, nlme, pbkrtest and between engines
└── inst/
    ├── CITATION  reference returned by citation("MLREML")
    ├── doc/      user guides in English and Spanish (PDF, version 1.0.0) installed with the package
    └── scale/    sparse_scale.R: timing of the sparse engine with simulated animal models
```

---

## Installation

MLREML contains C++ code that calls the BLAS and LAPACK libraries of R, so it is installed **from source** and needs a working build toolchain. Install the toolchain for your operating system first (step 1) and then the package (step 2). The same two installation options are given in Box 1 of the article.

### Step 1. Build tools required on each operating system

| Operating system | What you need | How to get it |
|---|---|---|
| **Windows** | **Rtools** (C/C++ and Fortran compilers) for your version of R | Download the Rtools version that matches your R version from <https://cran.r-project.org/bin/windows/Rtools/> and install it with the default options. |
| **macOS** (Apple silicon and Intel) | (a) **Xcode command line tools** (C/C++ compiler) **and** (b) the **GNU Fortran (gfortran)** build used by CRAN | (a) `xcode-select --install` in the Terminal. (b) The CRAN installer from <https://mac.r-project.org/tools/> (see the detailed instructions below). |
| **Linux** | C/C++ and Fortran compilers and the R headers | Debian/Ubuntu: `sudo apt install r-base-dev` · Fedora/RHEL: `sudo dnf install R-devel` (both include `gcc`, `g++` and `gfortran`). |

#### Detailed instructions for macOS

MLREML itself contains no Fortran code, but it links to the BLAS/LAPACK libraries shipped with R, and R for macOS adds the **GNU Fortran run-time libraries** (`-lgfortran`, `-lquadmath`, located in `/opt/gfortran`) to the link step of every package that uses BLAS/LAPACK. The Xcode command line tools do **not** include a Fortran compiler. Without GNU Fortran, the installation stops at the link step with an error such as

```
ld: warning: search path '/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/14.2.0' not found
ld: library 'gfortran' not found
clang++: error: linker command failed with exit code 1
ERROR: compilation failed for package 'MLREML'
```

(older versions of Xcode write `ld: library not found for -lgfortran`). To avoid it:

1. **Install the Xcode command line tools** (once). In the Terminal:

   ```bash
   xcode-select --install
   ```

   If they are already installed, the command says so; you can continue.

2. **Check which GNU Fortran your R expects.** In the Terminal:

   ```bash
   R --version | head -1      # your R version
   R CMD config FLIBS         # Fortran libraries used by your R
   ```

   With R ≥ 4.5.0 the output of `R CMD config FLIBS` contains `/opt/gfortran/lib/gcc/.../14.2.0`, that is, **GNU Fortran 14.2**; with R 4.3.0–4.4.3 it contains `12.2.0`, that is, **GNU Fortran 12.2**.

3. **Install the matching GNU Fortran from CRAN** (it is the compiler used to build R itself; one installer works on both Apple silicon and Intel Macs):

   | Your R version | Installer to download and double-click |
   |---|---|
   | R ≥ 4.5.0 (including R 4.6.x) | [`gfortran-14.2-universal.pkg`](https://mac.r-project.org/tools/gfortran-14.2-universal.pkg) |
   | R 4.3.0 – 4.4.3 | [`gfortran-12.2-universal.pkg`](https://mac.r-project.org/tools/gfortran-12.2-universal.pkg) |

   The installer puts the compiler in `/opt/gfortran`. Please use these CRAN installers rather than the `gcc`/`gfortran` of Homebrew or MacPorts, which are installed in other folders and are not the ones that R for macOS expects.

4. **Check the installation** (in a new Terminal window):

   ```bash
   /opt/gfortran/bin/gfortran --version     # should print GNU Fortran 14.2 (or 12.2)
   ls /opt/gfortran/lib | grep gfortran     # libgfortran.dylib, the library that the linker was missing
   ```

5. **Restart R or RStudio** and continue with Step 2.

Optionally, the toolchain can be checked from R with `install.packages("pkgbuild"); pkgbuild::check_build_tools(debug = TRUE)`, which should end with *"Your system is ready to build packages!"*.

### Step 2. Install MLREML

```r
install.packages(c("Matrix", "Rcpp"))   # dependencies (Matrix >= 1.6-0)

# Option 1: install directly from this repository (internet connection needed)
install.packages("https://github.com/osval78/MLREML/raw/main/MLREML_1.0.0.tar.gz",
                 repos = NULL, type = "source")

# Option 2: download MLREML_1.0.0.tar.gz (button "Download raw file" on its GitHub page),
#           set the working directory to the folder that contains it and install the local file
install.packages("MLREML_1.0.0.tar.gz", repos = NULL, type = "source")

library(MLREML)
?mlreml                                                          # help of the main function
system.file("doc", "MLREML-guide-en.pdf", package = "MLREML")   # installed user guide
```

> **Note.** This repository distributes the source package as the file `MLREML_1.0.0.tar.gz`; the package sources are not stored unpacked at the top of the repository. For this reason `remotes::install_github("osval78/MLREML")` and `devtools::install_github("osval78/MLREML")` do **not** work (they stop with *"Does not appear to be an R package (no DESCRIPTION)"*). Please use Option 1 or Option 2.

To check that the installation works, fit a small model:

```r
library(MLREML)
set.seed(1)
d <- data.frame(g = factor(rep(1:20, each = 3)))
d$y <- rnorm(20)[d$g] + rnorm(60)
fit <- mlreml(y ~ 1, random = ~ g, data = d, trace = FALSE)
summary(fit)$varcomp          # genotype and residual variance components
```

The examples of the guides also use the packages `agridat`, `BGLR`, `pedigreemm`, `lme4`, `nlme` and `pbkrtest`. None of them is needed to use MLREML:

```r
install.packages(c("agridat", "BGLR", "pedigreemm", "lme4", "nlme", "pbkrtest"))
```

For speed, use R with an optimised BLAS library (OpenBLAS, Intel MKL or Accelerate on macOS).

### Troubleshooting the installation

| Message | Cause and solution |
|---|---|
| `ld: library 'gfortran' not found`, `library not found for -lgfortran`, or `search path '/opt/gfortran/...' not found` (macOS) | GNU Fortran is missing or is not the version expected by R. Install the CRAN gfortran that matches your R version (Step 1, macOS, points 2–4) and restart R. |
| `xcrun: error: invalid active developer path` (macOS) | The Xcode command line tools are missing or were removed by a macOS update. Run `xcode-select --install`. |
| `Does not appear to be an R package (no DESCRIPTION)` | `remotes::install_github()` or `devtools::install_github()` was used. Use Option 1 or Option 2 of Step 2. |
| `WARNING: Rtools is required to build R packages` or `make: not found` (Windows) | Install the Rtools version that matches your R version and restart R. |
| `cannot open URL` / `download failed` | No internet access or a proxy blocks GitHub. Download `MLREML_1.0.0.tar.gz` with the browser and use Option 2. |
| `there is no package called 'Rcpp'` or `'Matrix'` | Run `install.packages(c("Matrix", "Rcpp"))` first. |

If the installation still fails, please open an [Issue](https://github.com/osval78/MLREML/issues) and include the output of `sessionInfo()` and the full error message.

---

## Quick start

### 1. Alpha-lattice trial: variance components, heritability and BLUPs

```r
library(MLREML)
data(john.alpha, package = "agridat")
d <- john.alpha

fit <- mlreml(yield ~ rep, random = ~ gen + rep:block, data = d)
summary(fit)                                  # variance components, logLik, AIC, BIC
vpredict(fit, H2 ~ V1 / (V1 + V3 / 3))        # heritability with delta-method SE
head(coef(fit, "random", "gen"))              # BLUP, PEV, SE and reliability

fit_f <- mlreml(yield ~ rep + gen, random = ~ rep:block, data = d)
wald(fit_f)                                   # Wald F tests, Kenward–Roger df
predict(fit_f, classify = "gen")$avsed        # adjusted means and average SED
```

### 2. Genomic prediction (GBLUP)

```r
data(wheat, package = "BGLR")
ids <- paste0("L", 1:599)
M <- wheat.X; rownames(M) <- ids
qc <- snp_qc(2 * M, maf = 0.05, callrate = 0.9)      # marker quality control
G  <- Gmatrix(qc$M, method = "VanRaden") + diag(1e-4, 599)
d1 <- data.frame(gen = ids, y = wheat.Y[, 1])

f1 <- mlreml(y ~ 1, random = ~ vm(gen, G), data = d1)   # general engine
g1 <- gblup(setNames(d1$y, ids), G)                      # exact spectral GBLUP
```

### 3. Multi-environment trial with a factor analytic model

```r
data(besag.met, package = "agridat")
fa1 <- mlreml(yield ~ county,
              random   = ~ county:rep + county:rep:block + fa(county, 1):id(gen),
              residual = ~ diag(county):units, data = besag.met, maxit = 100)
vc_matrix(fa1, "fa(county, 1):id(gen)")    # genetic covariance/correlation, loadings
```

### 4. Large animal models with the sparse engine

```r
# ped: pedigree (id, sire, dam); d: data with columns y, herd, sex and animal
Ai  <- Ainverse(ped)                                   # sparse inverse of A
fit <- mlreml(y ~ herd + sex,
              random = ~ vm(animal, Ai) + id(animal),  # A given through its inverse
              data = d, engine = "sparse")             # mixed model equations
```

### 5. Direct–maternal model

```r
# A <- Amatrix(ped); d: data with columns y, cg (contemporary group), animal and dam
fit <- mlreml(y ~ cg, random = ~ str(animal, dam, K = A) + id(dam), data = d)
vc_matrix(fit, "str(animal, dam, K = A)")$cor          # direct–maternal genetic correlation
```

---

## Model syntax at a glance

```r
mlreml(fixed, random = NULL, residual = NULL, data, start = NULL, fix = NULL,
       maxit = 50, tol = 1e-6, pev = TRUE, trace = TRUE, store = TRUE,
       engine = c("auto", "dense", "sparse"))
```

| Structure | Meaning |
|---|---|
| `id(f)` or `f` | independent effects |
| `vm(f, K)` | effects correlated by a known matrix: pedigree **A**, genomic **G**, single-step **H**, any kernel |
| `diag(f)` | a different variance for each level |
| `us(f)` | unstructured covariance matrix |
| `fa(f, k)` | factor analytic of order *k* |
| `cs(f)` | compound symmetry |
| `ar1(f)` | first-order autoregressive (actual distances if levels are numeric) |
| `expc(x)`, `sexp(x, y)` | exponential correlation for irregular covariates / 2-D coordinates |
| `pol(x, k)`, `leg(x, k)` | polynomial and Legendre bases (random regression) |
| `mm(f1, f2)` | multi-membership factor (e.g., GCA in diallels) |
| `str(f1, f2, K = K)` | correlated effects of several factors sharing **K** (direct–maternal models) |
| `units` | observation units |

Typical models:

| Analysis | MLREML |
|---|---|
| Alpha-lattice, random genotypes | `mlreml(y ~ rep, random = ~ gen + rep:block, data = d)` |
| Row × column spatial analysis | `mlreml(y ~ rep, random = ~ gen + units, residual = ~ ar1(row):ar1(col), data = d)` |
| Heterogeneous residuals by environment | `residual = ~ diag(env):units` |
| Genomic MET | `random = ~ us(env):vm(gen, G)` or `fa(env, 2):vm(gen, G)` |
| Animal model with repeatability | `random = ~ vm(animal, A) + id(animal)` |
| Bivariate model | `mlreml(cbind(y1, y2) ~ trait, random = ~ us(trait):vm(animal, A), residual = ~ us(trait):units, data = d)` |
| Random regression | `random = ~ us(leg(age, 2)):vm(animal, A)` |
| Diallel | `random = ~ mm(female, male) + cross` |

---

## Main functions

| Function | Purpose |
|---|---|
| `mlreml()` | Fit a linear mixed model by AI-REML (dense or sparse engine) |
| `summary()`, `coef()`, `fitted()`, `residuals()`, `vcov()`, `logLik()`, `plot()` | Standard methods; `coef(fit, "random")` gives BLUP, PEV, SE and reliability |
| `wald()` | Incremental and conditional Wald F tests (Kenward–Roger or residual df) |
| `predict()` | Predicted means, standard errors and SED matrix |
| `vpredict()` | Functions of variance components (heritability, correlations) with delta-method SE |
| `lrt()` | Likelihood ratio test, with boundary correction |
| `vc_matrix()` | Covariance/correlation matrices of `diag`, `us`, `fa` and `str()` terms; FA loadings |
| `variogram()` | Sample variogram of residuals of field trials |
| `pedigree_prepare()`, `inbreeding_coef()`, `Amatrix()`, `Ainverse()` | Pedigree tools |
| `snp_qc()`, `Gmatrix()`, `align_G()`, `Hmatrix()`, `bend()` | Genomic and single-step relationship matrices |
| `gblup()`, `cv_gblup()` | Fast spectral GBLUP and k-fold cross-validation |

---

## Two computing engines

* **Dense engine** – works with the marginal variance matrix **V**; supports every structure; comfortable up to about 5,000–6,000 records with 8 GB of memory.
* **Sparse engine** – works with Henderson's mixed model equations, a supernodal sparse Cholesky factorisation (CHOLMOD, through `Matrix`) and the Takahashi sparse inverse; supports `id`, `vm`, `diag`, `us` and `str()` terms with `units`, `diag(f):units` or `us(trait):units` residuals; fits animal models with hundreds of thousands of records and nested models with millions of records on a personal computer.

With `engine = "auto"` (default) the sparse engine is used when the data have more than 5,000 records and the model is supported. Both engines give the same likelihood, estimates, BLUPs and prediction error variances up to rounding error.

---

## Validation

The automatic tests of the package (`tests/testthat`; 18 tests with 65 expectations in version 1.0.0, all passed) compare MLREML with established software every time the package is checked:

* variance components, BLUPs, predicted means and SED against **lme4**;
* heterogeneous and autoregressive residuals and exponential correlations against **nlme**;
* Kenward–Roger F statistics and denominator df against **pbkrtest**;
* the general engine against an exact spectral **GBLUP**, and the sparse engine against the dense engine;
* the direct–maternal term `str()`, boundary handling, pedigree tools (A, A⁻¹, inbreeding) and a numerical check of the score at the optimum.

Chapters 14 and 15 of the user guides present a simulation study and the complete numerical validation.

---

## Documentation

| Language | User guide |
|---|---|
| English | [*MLREML: linear mixed models by average-information REML for plant and animal breeding. User guide, theoretical background and worked examples with real data — version 1.0.0*](MLREML_documentation_Espa%C3%B1ol_v1.0.0.pdf) |
| Español | [*MLREML: modelos lineales mixtos por REML con información promedio para el mejoramiento genético vegetal y animal. Guía de usuario, fundamentos teóricos y ejemplos con datos reales — versión 1.0.0*](MLREML_documentacion_ENGLISH_v1.0.0.pdf) |
| Português | [*MLREML: modelos lineares mistos por REML de informação média para o melhoramento de plantas e animais. Guia do usuário, fundamentos teóricos e exemplos resolvidos com dados reais — versão 1.0.0*](MLREML_documentacao_Portugues_v1.0.0.pdf) |
| 中文 | [*MLREML：用于植物和动物育种的平均信息 REML 线性混合模型——用户指南、理论基础及真实数据实例，版本 1.0.0*](MLREML_documentation_Chino_v1.0.0.pdf) |

The English and Spanish guides are also installed with the package:

```r
system.file("doc", "MLREML-guide-en.pdf", package = "MLREML")   # English
system.file("doc", "MLREML-guide-es.pdf", package = "MLREML")   # Spanish
```

Contents of the guides: (1) introduction and installation; (2) theory (mixed model, REML, AI algorithm, sparse engine, Kenward–Roger); (3) model syntax; (4–12) worked examples with real data (alpha-lattice oat trial, spatial wheat trial, maize MET with factor analytic models, GBLUP with 599 wheat lines, dairy animal model with pedigree, single-step GBLUP, repeated measures and random regression, diallel, unreplicated early-generation trial); (13) advanced use and troubleshooting; (14) simulation study; (15) validation, performance and limitations; (16) function summary and ASReml-R syntax equivalences; appendix with the reference manual.

---

## Limitations

* Only Gaussian responses (no generalised linear mixed models).
* Structures with correlation parameters (`cs`, `ar1`, `expc`, `sexp`) and `fa` are available only in the dense engine.
* Two random terms can be correlated only through `str()` (shared relationship matrix).
* Kenward–Roger df are available for the dense engine only (the sparse engine uses n − p); the Satterthwaite approximation is not implemented.
* Likelihood ratio tests between factor analytic models of different order have no standard reference distribution; `lrt()` reports the statistic, and calibration requires a parametric bootstrap.
* There is no built-in `weights` argument for two-stage analyses (propagating the uncertainty of first-stage BLUEs).
* A few ASReml-specific structures, such as ARMA and random cubic splines, are not implemented.
* MLREML is not a replacement for specialised programs (ASReml, BLUPF90, DMU, WOMBAT) used in national evaluations with tens of millions of records.

MLREML and ASReml are not interchangeable; only the features compared explicitly in the guides should be considered equivalent.

---

## How to cite

If you use MLREML, please cite the user guide (and the article, once published). The same reference is returned in R by `citation("MLREML")`.

> Montesinos-López, O. A., Montesinos-López, A., Crossa, J., Montesinos-López, J. C., Vitale, P., and Ortiz, R. (2026). *MLREML: Linear Mixed Models by Average-Information REML for Plant and Animal Breeding. User Guide, Theoretical Background and Worked Examples with Real Data, Version 1.0.0.* Available at: https://github.com/osval78/MLREML

The Spanish, Portuguese and Chinese guides have the same six authors and can be cited in the same way with their own titles (see [Documentation](#documentation)).

> Montesinos-López, O. A., Crossa, J., Montesinos-López, J. C., Montesinos-López, A., Vitale, P., and Ortiz, R. Linear mixed models with the MLREML R package for plant and animal breeding. *Manuscript under review.*

---

## Authors

* **Osval A. Montesinos-López** – Facultad de Telemática, Universidad de Colima, México (maintainer)
* **Abelardo Montesinos-López** – CUCEI, Universidad de Guadalajara, México
* **José Crossa** – Colegio de Postgraduados, México
* **J. Cricelio Montesinos-López** – Independent consultant, México
* **P. Vitale** – International Maize and Wheat Improvement Center (CIMMYT), México
* **Rodomiro Ortiz** – Swedish University of Agricultural Sciences (SLU), Sweden

Questions, bug reports and suggestions are welcome through the [Issues](https://github.com/osval78/MLREML/issues) page of this repository.

## License

MLREML is free software distributed under the **GNU General Public License v3.0 (GPL-3)**.
