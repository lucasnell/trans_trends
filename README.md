

Code for analyses related to the paper, "Quantifying community responses to 
environmental variation from replicate time series."
The R package that fits the models can be found at 
<https://github.com/lucasnell/trans_trends_pkg>.


# File organization:

- `./analysis/`
    - `analysis.R` creates figures 1-3 and table 1
    - `data_fit.csv` is a cleaned version of the data used for model fitting
    - `model_fit.R` cleans the data, fits the model using the 
      `TransTrendsPkg` package, and summarizes the model output to 
      some temporary files
    - `output/`
        - `fit_sum.csv` contains summaries of the parameter-estimate 
          distributions from the model, without any cleaning of parameter names
    - `pca_funs.R` contains a bunch of helper functions for the PCA
- `./data/`
    - `clean_data.R` cleans the raw data
    - `myv_arth.csv` is the cleaned version of data
    - `myvatn_infalls.csv` is raw data from infall traps that estimate
      midge deposition
    - `myvatn_pitfalls.csv` is raw data from predatory arthropod pitfall traps

