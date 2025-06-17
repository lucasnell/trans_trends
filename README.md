
<!-- README.md is generated from README.Rmd. Please edit that file -->

# trans_trends

[![DOI](https://zenodo.org/badge/227422889.svg)](https://zenodo.org/badge/latestdoi/227422889)

Code and data for the paper, “Heterogeneous responses of tundra
arthropods to spatiotemporal variation in resource subsidies.” The R
package that fits the models can be found at
<https://github.com/lucasnell/trans_trends_pkg>.

This repository was created and is maintained by [Lucas A.
Nell](https://github.com/lucasnell).

# Organization

These are the main top-level files/folders:

``` bash
.
├── README.md
├── trans_trends.Rproj
├── data
├── figures
├── rds_files
└── scripts
```

The following files/folders should be present:

- `README.md`: this file
- `trans_trends.Rproj`: file saving this RStudio Project’s preferences
- `data`: Folder containing all datasets associated with this paper.
- `figures`: Folder containing figures created from this paper. Note
  that this folder should be empty and is populated by running the R
  scripts.
- `rds_files`: Folder containing RDS files created for model fits and
  cleaned data. Note that this folder should be empty and is populated
  by running the R scripts.
- `scripts`: Folder containing R scripts for data cleanup, model
  fitting, and creating figures / tables.

Note that all folders have their own `README.md` files that should be
referenced for specific files.

# Replicating the R environment

I used R version 4.5.0 (2025-04-11) (platform: aarch64-apple-darwin20)
for all my scripts.

This project uses the `renv` package, so if you want to use this, you
must first install it:

``` r
install.packages("renv")
```

Then to install all the packages I used for these analyses, you can
simply run the following while having this project’s main directory as
your working directory:

``` r
renv::restore()
```

If you’d rather avoid `renv`, then you can install all the packages (in
the versions I used) this way:

``` r
pkgs <- c("abind@1.4-8", "backports@1.5.0",
          "bayesplot@1.12.0", "BH@1.87.0-1",
          "callr@3.7.6", "checkmate@2.3.2",
          "cli@3.6.5", "cpp11@0.5.2",
          "desc@1.4.3", "distributional@0.5.0",
          "dplyr@1.1.4", "farver@2.1.2",
          "generics@0.1.4", "ggplot2@3.5.2",
          "ggridges@0.5.6", "ggspatial@1.1.9",
          "glue@1.8.0", "gridExtra@2.3",
          "gtable@0.3.6", "inline@0.3.21",
          "isoband@0.2.7", "labeling@0.4.3",
          "lifecycle@1.0.4", "lme4@1.1-37",
          "loo@2.8.0", "magrittr@2.0.3",
          "matrixStats@1.5.0", "minqa@1.2.8",
          "nloptr@2.2.1", "numDeriv@2016.8-1.1",
          "patchwork@1.3.0", "pillar@1.10.2",
          "pkgbuild@1.4.8", "pkgconfig@2.0.3",
          "plyr@1.8.9", "posterior@1.6.1",
          "processx@3.8.6", "ps@1.9.1",
          "purrr@1.0.4", "QuickJSR@1.8.0",
          "R6@2.6.1", "rbibutils@2.3",
          "RColorBrewer@1.1-3", "RColorBrewer@1.1.3",
          "Rcpp@1.0.14", "RcppEigen@0.3.4.0.2",
          "RcppParallel@5.1.10", "Rdpack@2.6.4",
          "reformulas@0.4.1", "remotes@2.5.0",
          "reshape2@1.4.4", "rlang@1.1.6",
          "rstan@2.32.7", "rstantools@2.4.0",
          "scales@1.4.0", "sf@1.0.21",
          "StanHeaders@2.32.10", "stringi@1.8.7",
          "stringr@1.5.1", "tensorA@0.36.2.1",
          "tibble@3.3.0", "tidyr@1.3.1",
          "tidyselect@1.2.1", "tidyverse@2.0.0",
          "utf8@1.2.6", "vctrs@0.6.5",
          "viridisLite@0.4.2", "withr@3.0.2")
install.packages(pkgs)
remotes::install_github("lucasnell/trans_trends_pkg@1.0.2")
```

Note that these only install the proper versions of the packages I
manually installed (or were dependencies of `TransTrendsPkg`), so other
dependencies might vary from what I used.

Similarly, but without version numbers at all:

``` r
pkgs <- c("abind", "backports", "bayesplot", "BH",
          "callr", "checkmate", "cli", "cpp11",
          "desc", "distributional", "dplyr", "farver",
          "generics", "ggplot2", "ggridges", "ggspatial",
          "glue", "gridExtra", "gtable", "inline",
          "isoband", "labeling", "lifecycle", "lme4",
          "loo", "magrittr", "matrixStats", "minqa",
          "nloptr", "numDeriv", "patchwork", "pillar",
          "pkgbuild", "pkgconfig", "plyr", "posterior",
          "processx", "ps", "purrr", "QuickJSR",
          "R6", "rbibutils", "RColorBrewer", "RColorBrewer",
          "Rcpp", "RcppEigen", "RcppParallel", "Rdpack",
          "reformulas", "remotes", "reshape2", "rlang",
          "rstan", "rstantools", "scales", "sf",
          "StanHeaders", "stringi", "stringr", "tensorA",
          "tibble", "tidyr", "tidyselect", "tidyverse",
          "utf8", "vctrs", "viridisLite", "withr")
install.packages(pkgs)
remotes::install_github("lucasnell/trans_trends_pkg")
```
