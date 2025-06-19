# `scripts`

Folder containing R scripts that run analyses and create figures and tables
for the manuscript.

Folder contents:

```
.
├── 00-preamble.R
├── 01-sample-map.R
├── 02-clean-data.R
├── 03-catch-rates-plots.R
├── 04-fit_models.R
├── 05-coefficients.R
├── 06-pca.R
└── README.md
```


File descriptions:


- `00-preamble.R`: preamble that loads packages, defines objects, and does other
  miscellaneous tasks that are used in other scripts. This file is sourced 
  from all other scripts in here.
- `01-sample-map.R`: creates maps showing where samples were taken at Lake Myvatn
- `02-clean-data.R`: cleans the raw data and preps it for model fitting and plotting
- `03-catch-rates-plots.R`: creates figure showing catch rates through time and 
  across distance from Lake Myvatn
- `04-fit_models.R`: fits the model using the `TransTrendsPkg` package
- `05-coefficients.R`: creates plots of coefficient estimates from the model fit.
  It also contains code that shows, for each taxon, the 
  autoregressive parameter and the effect of a 10-fold midge increase.
- `06-pca.R`: conducts PC analyses and creates resulting plots
- `README.md`: this file
