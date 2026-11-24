# Code and statistical review: issues

Review of `scripts/` and of the package `trans_trends_pkg` (Stan model and
`armm()`), using the committed `data/myv_arth.csv`,
`rds_files/data-fit.rds`, and `rds_files/model-fit.rds`.

Every number below comes from a script in this folder, run from the project
root. Each script's output is saved next to it:

| Script | Output | Covers |
|---|---|---|
| `01-check-data.R` | `01-check-data-output.txt` | E1–E3, S4 (window by year) |
| `02-check-corrected-cleaning.R` | `02-check-corrected-cleaning-output.txt` | E2–E4, S4 (sensitivity) |
| `03-check-model.R` | `03-check-model-output.txt` | Sampler health, S1, S2, S5, S6 |
| `04-check-year-shocks.R` | `04-check-year-shocks-output.txt` | S3 |
| `05-overlap-checks.R` | `05-overlap-checks-output.txt` | E2 status, E5, replacement overlap checks |
| `06-overlaps-2012.R` | `06-overlaps-2012-output.txt` | E4 overlap groups and their counts |
| `07-e5-infall-duplicates.R` | `07-e5-infall-duplicates-output.txt` | E5 duplicate rows, curation notes, missing 2012 samples |
| `08-e4-sum-check.R` | `08-e4-sum-check-output.txt` | E4: is any sample a sum or average of others |
| `09-test-overlap-checks.R` | `09-test-overlap-checks-output.txt` | Tests of the new checks in `00-clean-data-round1.R` |
| `10-date-overlaps-all.R` | `10-date-overlaps-all-output.txt` | E6, E7: overlapping date ranges, all years, both trap types |
| `11-check-e6-fix.R` | `11-check-e6-fix-output.txt` | E6: verifies the btl 5 m 2013 row removal |
| `12-e5-assign-pairs.R` | `12-e5-assign-pairs-output.txt` | E5 pairs 1–2: assignment from other distances, calibration, sensitivity |
| `13-e5-random-forest.R` | `13-e5-random-forest-output.txt` | E5 pairs 1–2: random-forest assignment and calibration |
| `14-e5-method-comparison.R` | `14-e5-method-comparison-output.txt` | E5 pairs 1–2: random forest vs. quantile forest vs. GAM on 20,000 test pairs |
| `15-e5-drop-cost.R` | `15-e5-drop-cost-output.txt` | E5: cost of dropping both rows of pairs 1–2 |
| `16-e5-method-comparison-2.R` | `16-e5-method-comparison-2-output.txt` | E5 pairs 1–2: variance-modelling GAM, joint scoring, stacked classifier |

Status key: **Open**, **Resolved**, **Partly resolved**.

## Sampler health (no issue)

No divergences, no transitions at maximum tree depth, energy diagnostic
(E-BFMI) 0.34–0.43 per chain, maximum R-hat 1.009, and minimum effective
sample size 391 (`sig_proc`). The per-taxon slopes rebuilt in
`05-model-summaries.R` (`re_draws`) match Stan's own `beta` exactly, and
`rstan::extract(pars = c("sig_proc", "sig_obs", "sig_beta"))` returns them in
the requested order, so the `set_names()` labels in `extract_effect_info()`
are correct. (03)

---

## Code errors

### E1. Raw file names in `00-clean-data-round1.R` did not exist — Resolved

The script read `pitfalls-08-17.csv`, `pitfalls-11-19.csv`,
`infalls-08-17.csv`, and `infalls-13-19.csv`, but the files were named
`myvatn_pitfalls.csv`, `pitfalls.csv`, `myvatn_infalls.csv`, and `infalls.csv`.

- **Status:** The raw files have been renamed to match the script
  (committed in `186763b`). Three are byte-identical to the old versions;
  `infalls-08-17.csv` differs from `myvatn_infalls.csv` only by a removed
  byte-order mark and an added final newline. All four files now exist (01).
- `data/README.md` now names the metadata files as they exist
  (`metadata-infall.csv`, `metadata-pitfall.csv`); every file in its tree
  matches `data/` and `data/raw/`.

### E2. Pitfall samples from 2013–2017 are counted twice — Resolved in script (uncommitted; `myv_arth.csv` not regenerated)

`full_join(p1, p2, by = ...)` (`00-clean-data-round1.R` lines 56–59) matches
on every column, including `acar_total`. In `pitfalls-11-19.csv`,
`acar_total = trom_adult + trom_juv + acar` is NA for every 2011–2017 sample
because `trom_adult` and `trom_juv` are NA, while `pitfalls-08-17.csv` has
values. Overlapping samples therefore fail to match and appear twice: 377
Mývatn samples are duplicated (01).

The overlap check (lines 37–53) missed this for two reasons: it compares only
columns 7–14 one at a time, and `x[[2]] != x[[3]]` is NA when either side is
NA, which `filter()` drops. Counting NA-versus-value as a mismatch gives 549
mismatches in `acar_total` and none elsewhere (01).

Effect on `data/myv_arth.csv` (02):

| Year | Original / corrected `season_days` | Original / corrected count |
|---|---|---|
| 2013–2016 | 2 | 2 |
| 2017 | 1.73 | 1.91 |
| Other years | 1 | 1 |

- 2013–2016 catch rates (count / days) are unchanged, but the Poisson part of
  the model sees twice the individuals actually caught.
- 2017 is only partly duplicated, so its catch rates are shifted.
- **Fix:** join on the sample key only, for example
  `bind_rows(p1, anti_join(p2, p1, by = c("lakeid", "dist", "trans", "setdate", "coldate")))`.
  Mites (`acar`) are dropped from the analysis, so the `acar_total`
  disagreement does not need resolving for the analysed taxa.
- **Status:** `00-clean-data-round1.R` lines 56–57 now use this join. Run
  against the current raw files, it gives 1,170 pitfall rows (933 from
  `pitfalls-08-17.csv` plus 237 only in `pitfalls-11-19.csv`) with no
  duplicated sample among rows that survive cleaning. The resulting
  `season_days` is half the committed value in 2013–2016 and 1/1.73 in 2017,
  and equal in other years. The only plot-years whose `season_days` still
  exceeds their sampling window are in 2012 (E4). Where both files have a
  sample, they agree on every analysed column (05).
- **Still to do:** `write_csv()` (line 182) is commented out, so
  `data/myv_arth.csv`, `rds_files/data-fit.rds`, and the model fit still
  contain the duplicates.
- **Overlap checks (lines 37–53 and 130–146):** in addition to the NA problem
  above, they filter on `!is.na(x[[7]])`, which is a column from the first
  file and so never removes unmatched rows; they join on `daysout`, so a
  sample whose `daysout` differs between files is never compared; they don't
  check that each sample appears once within each file; and they only print.
  These are now replaced in `00-clean-data-round1.R` by `check_unique_key()`
  and `compare_overlap()`, which stop the script on an unexpected duplicate or
  mismatch. The known `acar_total` disagreement and the four E5 infall
  duplicates are explicitly allowed (uncommitted). Tests are in (09).

### E3. Samples with NA counts still add their days to sampling effort — Open

`summarize(count = sum(round(count), na.rm = TRUE), season_days = sum(daysout))`
drops an NA count but keeps that sample's days, which is the same as treating
the sample as zero individuals caught. The infall summary (lines 170–177) does
the same for midges.

**Pitfalls** (01, 02):
- 19 samples have NA counts for five taxa and 12 for wolf spiders, totalling
  321 and 219 trap-days across 15 and 9 plot-years.
- 36 of 1,152 analysed rows shift by more than 0.05 on the
  `log1p(count / season_days)` scale.
- Worst case: hag, 5 m, 2018 has 11 counted days, not 74, so its catch rates
  are understated about 7-fold (sheet weavers: 14 / 74 versus 14 / 11).

**Midges** (02):
- 78 infall samples have NA midge totals, affecting 41 of 192 analysed
  plot-years.
- In 2014 every plot's midge rate is too low by about 0.55 log units (about
  42%). Other affected plot-years shift by up to 0.49 log units.
- Midges are the key predictor, so this is error concentrated in one year,
  which can bias both the midge and time coefficients.

**Fix:** drop NA samples before summing, or sum `daysout[!is.na(count)]`.

### E4. Overlapping pitfall deployments in 2012 — Open

31 samples in 2012 have a `setdate` earlier than the previous sample's
`coldate` at the same plot (02). This is why 2012 has 85 sampling days inside
an 82-day window (day-of-year 140–222). It is either a date-entry error or two
traps per plot; it needs checking against the field sheets.

Details are in `E4-2012-overlaps.md`. Briefly, 54 samples fall into 23 overlap
groups, all from `pitfalls-08-17.csv`. Counts differ within every group, so
the samples are distinct catches, not duplicate entries (06). No 2012 sample equals the sum of two or more other samples from its plot, or the average of two, on any count column (08).

### E5. Conflicting duplicate samples within `infalls-08-17.csv` — Open

Four Mývatn infall samples appear twice within `infalls-08-17.csv` with
different values (05, 07). The current cleaning keeps both copies, so their
midges are summed and their days counted twice.

Details, analyses, limitations, and references are in
`E5-pair-assignment.md`. In brief:

- **Pairs 3–4** (non 150 m and vin 5 m, 2008): one copy of each is all NA (a
  "missing data" placeholder and a possibly fallen trap). Drop those rows.
- **Pairs 1–2** (hag 5 m and non 50 m, 2012): curation notes say one copy of
  each is probably a mislabelled sample from btl 50 m 2012/06/19–07/03 and
  kal 5 m 2012/07/17–07/31, which are both missing. Field sheets aren't
  available. Prediction from other distances, neighbouring samples, and
  season gives these assignments:
  - **Pair 2:** line 942 = non 50 m and line 944 = kal 5 m, at 92–95%
    accuracy on known test pairs at that evidence strength. The choice
    matters (0.5–0.7 log units).
  - **Pair 1:** line 864 = hag 5 m, at 81–87%. The premise is doubtful,
    because neither row looks like hag 5 m in early June, but the choice
    changes rates by under 0.05 log units.
- **Decision pending:** drop pair 1 and assign pair 2 (the middle option),
  assign both, or drop both. Dropping pair 2 leaves kal 5 m with no samples
  after mid-July and probably halves its 2012 midge rate.

The infall `full_join` on all columns (lines 149–151) has the same weakness
as the old pitfall join. It currently creates no cross-file duplicates,
because the two files agree on every overlapping sample (05), but a key-based
join would protect against future edits.

### E6. Overlapping infall deployment periods, including a mis-dated copy at btl 5 m in 2013 — Partly resolved

The duplicate and mismatch checks compare samples by exact set and collection
dates, so they can't see samples whose date ranges overlap. A within-plot
check (10) finds 16 infall samples that start before an earlier sample at the
same plot ended: 5 at btl 5 m in 2013 (all caused by one mis-dated row), 4
from the E5 pairs, 6 at Kálfaströnd (now E7), and 1 at flg 500 m.

- **btl 5 m, 2013 (the important one):** `infalls-13-19.csv` line 43 has a
  sample dated 2013/06/05 to 2013/08/14 with `daysout` 14 and counts
  identical to the 2013/07/31 to 08/14 sample (big 0, small 1,073, simuliids
  1). It's almost certainly a copy of that sample with the wrong set date.
  `infalls-08-17.csv` doesn't have it; its curation note mentions "two
  entries with the 2013-06-05 set date". Because the dates differ, it isn't
  matched and is kept. Days are recomputed from the dates, so it adds 1,073
  midges and 70 days. `data/myv_arth.csv` has 10,190 midges over 140 days
  (72.8/day) for this plot-year; without the copy it is 9,117 over 70 days
  (130.2/day), 0.58 log units higher.
- **E5 pairs:** 4 of the 16 are the E5 duplicates (same dates).
- **Kálfaströnd:** see E7.
- **flg 500 m, 2008:** a 1-day overlap. This plot isn't analysed, so no action
  is needed.

Pitfalls have no overlaps outside 2012 (E4).

**Fix:** drop line 43 of `infalls-13-19.csv` (or exclude it in
`00-clean-data-round1.R` with a comment), and add a within-plot date-overlap
check to `00-clean-data-round1.R` with the known E4, E5, and E7 cases
listed as allowed.

- **Status:** `00-clean-data-round1.R` lines 168–171 now drop the mis-dated
  btl 5 m row from `inf2` (uncommitted). Exactly one row is removed, and the
  rebuilt btl 5 m 2013 values are 9,117 midges over 70 days, as expected (11).
  `data/myv_arth.csv` isn't regenerated yet.
- **Still open:** the within-plot date-overlap check.

### E7. Overlapping Kálfaströnd infall dates in 2008 and 2011 — Open

At kal 5, 50, and 150 m, one infall sample starts before the previous one
ends, in both 2008 and 2011. All rows are in `infalls-08-17.csv` (10):

| Year | Earlier sample (lines) | Later sample (lines) | Overlap |
|---|---|---|---|
| 2008 | 06/17–06/26 (79–81) | 06/23–07/08 (110–112) | 3 days |
| 2011 | 07/05–07/20 (779–781) | 07/19–08/02 (803–805) | 1 day |

- **Recorded `daysout`:** in 2008, the `daysout` values (9 and 15) match the
  dates as written, so the overlap is internally consistent. In 2011, kal 5 m
  has `daysout` 7 and 12, where the dates give 15 and 14. kal 50 m and 150 m
  match their dates.
- **Pitfall visit dates at the same plots:** pitfalls were collected and
  reset on 2008/06/26 (next set 07/08) and collected on 2011/07/20. That
  suggests the infall set dates should be 06/26 and 07/20, not 06/23 and
  07/19. It's suggestive only, since infall and pitfall traps weren't
  necessarily serviced on the same day.
- **Comment:** kal 50 m 2011/07/05–07/20 (line 780) has "sample missing??"
  and "Refers to archived sample?", but it has counts.
- **Effect:** days are recomputed from the dates, so these plot-years each
  have 3 (2008) or 1 (2011) extra midge-trap-days, out of 75 and 78. That's a
  small bias in midge catch rate.

**Fix:** confirm the set dates against field sheets. If 06/26 and 07/20 are
right, correct them in the raw file or in `00-clean-data-round1.R` with a
comment.

---

## Statistical issues

### S1. `05` reads coefficients as `change = FALSE` effects, but the model uses `change = TRUE` — Open

With `change = TRUE`, the Stan model is

1. μ₁ = x₁ · β (first year of each series)
2. μₜ = xₜ · β + φ μₜ₋₁ + σ_proc εₜ (later years)

Here μ is the latent log catch rate, x the predictors, β the coefficients, φ
the taxon's autoregressive parameter, σ_proc the process-error standard
deviation, and ε a standard-normal deviate. A sustained one-unit change in a
predictor shifts the long-run level of μ by β / (1 − φ), not β. The
derivation commented in `05` (log y = β₀ + β₁ T̂) describes the
`change = FALSE` model.

Posterior medians (03):

| Taxon | φ | `05` total % change | Long-run total % change | `05` 10× midges | Long-run 10× midges |
|---|---|---|---|---|---|
| Ground beetles | 0.35 | 58.5 | 103.5 | 1.071 | 1.111 |
| Rove beetles | 0.35 | 27.7 | 46.0 | 1.044 | 1.069 |
| Harvestmen | 0.39 | −12.4 | −19.5 | 1.058 | 1.096 |
| Ground spiders | 0.04 | −23.3 | −24.2 | 1.590 | 1.625 |
| Sheet weavers | 0.16 | −58.9 | −65.9 | 0.850 | 0.821 |
| Wolf spiders | 0.48 | −74.2 | −92.7 | 0.985 | 0.972 |

For midges, both values are legitimate but answer different questions: β is
the effect of a one-year pulse in midges, and β / (1 − φ) is the effect of a
sustained change. The paper needs to state which one it reports. For a trend,
the long-run value is the relevant one. The one-year percent changes in `05`
have the same issue.

### S2. The first year of each series is inconsistent with the model's dynamics — Open

The first year is set to μ₁ = x₁ · β, with no autoregressive term and no
process error (`armm_ss_lnp.stan` lines 60 and 68). Under `change = TRUE` the
series then drifts toward x · β / (1 − φ), so a series with no trend still
changes over its first few years by about x · β · φ / (1 − φ) log units.

Running the model's own recursion with zero trend, taxon-level intercepts,
and midges and distance at their means gives first-to-last-year changes of
(03):

| Taxon | Change from start-up drift alone, median [95%] |
|---|---|
| Ground beetles | −62.9% [−81.7, −35.0] |
| Rove beetles | −61.1% [−79.3, −37.0] |
| Harvestmen | −46.6% [−69.1, −23.8] |
| Ground spiders | −8.9% [−25.7, −1.4] |
| Sheet weavers | −4.5% [−20.5, 2.0] |
| Wolf spiders | −10.9% [−42.6, 25.5] |

Plot and site intercepts change these magnitudes, so they show the size and
direction of the problem, not the bias in the fitted time coefficients.
Because time is 0 in the first year, the start-up drift competes with the
time trend; the large observation error (σ_obs ≈ 0.84) absorbs part of the
mismatch.

**Fix:** start each series from its long-run distribution (mean
x₁ · β / (1 − φ), standard deviation σ_proc / √(1 − φ²)), or refit with
`change = FALSE`. Only a refit shows how much the trends move. With
`change = FALSE` the first-year mean is correct, but its variance should
still be the long-run variance.

### S3. No shared year-level variation, so time and midge intervals are probably too narrow — Open

Process errors are independent across the 96 time series (taxon × site ×
distance), so the 16 plots of a taxon are treated as independent within each
year. The residuals show they are not (04):

| Residual | R² of year × taxon | Expected if independent | p |
|---|---|---|---|
| Process-error deviates (posterior means) | 0.224 | 0.062 | 5 × 10⁻²⁵ |
| Observation residuals (log scale) | 0.262 | 0.062 | 2 × 10⁻³⁶ |

Averaged over 16 plots, taxon-by-year residuals reach ±1 log unit (for
example, wolf spiders +1.02 in 2013 and −0.93 in 2017).

58% of the variance in `midges_z` is among years (04). So the midge and time
effects rest mostly on 12 years, not 192 plot-years, and their intervals are
probably too narrow.

`armm()` cannot currently fit a taxon × year random intercept: `get_ts_info()`
requires every random-effect grouping to be constant within a time series.
Adding one means changing the Stan code.

### S4. The within-season sampling window changed over the years — Open

Median end of pitfall sampling ranges from day-of-year 204 (2009) to 230
(2017). The correlation between end day-of-year and year is 0.65; for start
day-of-year it is −0.01 (01). Dividing by sampling days (or using
`offset(log(season_days))`) assumes catch per day is constant across the
season.

Rough sensitivity check: least-squares slope of log catch rate on year within
plots, using all samples versus only samples set on or after day-of-year 160
and collected on or before day 205. The restricted set keeps 179 of 192
plot-years. Standard errors ignore autocorrelation; this is not the paper's
model. (02)

| Taxon | Slope, all samples | Slope, common window |
|---|---|---|
| Ground beetles | 0.018 | −0.002 |
| Rove beetles | −0.010 | 0.002 |
| Harvestmen | −0.043 | −0.016 |
| Ground spiders | −0.033 | −0.036 |
| Sheet weavers | −0.080 | −0.046 |
| Wolf spiders | −0.178 | −0.182 |

Options: restrict to a common window, or model individual samples with
day-of-year as a covariate.

### S5. Effect sizes are reported as posterior medians without intervals — Open

`05` prints medians only. The 95% intervals for total percent change (03)
include zero for rove beetles [−17.6, 92.5], harvestmen [−43.6, 27.1], and
ground spiders [−52.8, 24.6]. The posterior probability of a positive time
effect is 0.87 (rove beetles), 0.25 (harvestmen), and 0.14 (ground spiders).
These shouldn't be presented as findings without their intervals. Sheet
weavers and wolf spiders (probability 0 of a positive time effect) and ground
beetles (0.985) are the clear trends. For the midge effect, the probability
of a positive coefficient is 1.00 for ground spiders, 0.075 for sheet
weavers, and 0.44–0.73 for the other four taxa.

### S6. Priors drive some among-taxa parameters — Open

Contraction (1 − posterior variance / prior variance; 1 means the data
determine the value, 0 means the prior does) (03):

| Parameter | Contraction |
|---|---|
| Among-taxa SD of intercepts | 0.585 |
| Among-plot SD of intercepts | 0.978 |
| Among-site SD of intercepts | 0.945 |
| Among-taxa SD of midge slopes | 0.935 |
| Among-taxa SD of time slopes | 0.916 |
| Among-taxa SD of distance slopes | 0.806 |
| φ (six taxa) | 0.763–0.981 |
| σ_proc, σ_obs | 0.916, 0.985 |

The `Gamma(1.5, 3)` prior (mean 0.5) shapes the among-taxa intercept SD. With
only six taxa, the among-taxa SDs in the `coefs-mean-sd` figure should come
with a prior-sensitivity check.

---

## Package (`trans_trends_pkg`)

### P1. `log_lik` leaves out the offset — Open

In `inst/stan/armm_ss_lnp.stan`, the model uses
`lognormal(ly_pred + offset, sig_obs)` but `log_lik` uses
`lognormal_lpdf(y_pred[i] | ly_pred[i], sig_obs)`. The LOO and WAIC printed by
`summary()` (called in `04-fit_models.R`) are therefore wrong. They are also
computed given the per-observation latent rates `y_pred`, which makes
leave-one-out unreliable even after the offset is fixed.

### P2. Documentation doesn't match the Stan code — Open

- The `priors` argument is accepted but never passed to Stan.
- The docs describe Normal priors on φ (μ = 0, σ = 0.5) and say
  `ar_bound = FALSE` allows φ > 1. Stan fixes φ to [0, 1] with a `Beta(2, 2)`
  prior, and `p_bound` is unused. Negative autocorrelation is impossible
  (ground spiders' φ sits near the bound: median 0.04).
- The intercept and fixed-effect prior is `Normal(0, 5)` in `armm_ss_lnp.stan`
  and `Normal(0, 1)` in `armm_ss.stan`, also undocumented.

### P3. The missing-levels check is backwards — Open (does not affect this fit)

In `check_len_sort_data()`, `length(levels(x)) < length(unique(x))` can never
be true; it should be `>`. Combined with `form_info_w_rand()` offsetting
`g_mat` columns by `max()` of the previous column rather than its number of
levels, an unused factor level would misalign random-effect indices. All
levels are present in this dataset.

---

## Minor

### M1. Convoluted `time` construction in `02-clean-data.R`

`factor(year, levels = c(1:max(year))) |> as.numeric()` creates a factor with
2,019 levels and returns `year`. `time = year - min(year)` is equivalent.

### M2. Results recorded only as comments in `05`

The effect-size tables in `05-model-summaries.R` exist only as pasted
comments. Writing them to a file would tie each number in the paper to a
saved output.
