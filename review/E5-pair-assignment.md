# E5: conflicting duplicate infall samples

*Note: dates are written year/month/day (e.g., 2012/06/05), or month/day when
the year is clear from context (e.g., 06/05).*

Scripts and saved outputs are in `review/` and are cited by number below, for
example (12) is `12-e5-assign-pairs.R` and `12-e5-assign-pairs-output.txt`.

## The problem

Four Mývatn infall samples appear twice within `infalls-08-17.csv` with the
same transect, distance, set date, and collection date, but different values
(05, 07). `unique()` doesn't remove them, and none of these dates is in
`infalls-13-19.csv`, so they are in the original data rather than created by
the merge. The current cleaning keeps both copies, so their midges are summed
and their days counted twice; the NA copies also feed E3.

| Pair | Plot | Dates | Line | Big | Small | `fract_smch` | `simu` |
|---|---|---|---|---|---|---|---|
| 1 | hag 5 m | 2012/06/05–06/19 | 863 | 4 | 291 | 1 | NA |
| | | | 864 | 7 | 468 | 1 | 5 |
| 2 | non 50 m | 2012/07/17–07/31 | 942 | 0 | 360 | 1 | 0 |
| | | | 944 | 0 | 601 | 0.229 | 0 |
| 3 | non 150 m | 2008/06/04–06/10 | 26 | NA | NA | NA | NA |
| | | | 27 | 26 | 173 | 1 | 5 |
| 4 | vin 5 m | 2008/05/30–06/09 | 17 | 1,842.5 | 28,482.4 | 1 | 0 |
| | | | 18 | NA | NA | 1 | NA |

Line numbers count the header as line 1. "Big" and "small" are raw counts
(`bgch`, `smch`) before dividing by the subsample fraction.

The raw file's `comments.2` column already flags all four pairs.

- **Pair 1:** *"One may be BTL 50 from the next sampling event, which is
  missing. No way to tell."* btl 50 m has no 2012/06/19–07/03 sample (5
  samples in 2012), and hag 5 m has 7, where most plots have 6.
- **Pair 2:** *"One of these is possibly kal 5. No way to tell."* kal 5 m has
  no 2012/07/17–07/31 sample (5 samples in 2012), and non 50 m has 7.
- **Pair 3:** line 26 is a "missing data" placeholder. Line 27's collection
  date was edited from 06/11 to 06/10, but its `daysout` is still 7, while
  the dates give 6.
- **Pair 4:** line 18 is a possibly fallen trap ("rubber band failure"); its
  zeros were changed to NA.

Pairs 3 and 4 are straightforward: drop the NA rows. Field sheets aren't
available for pairs 1 and 2, so the rest of this note is about assigning
those rows from the data. Questions sent to a collaborator are in
`E5-infall-duplicates.md`.

## What dropping both rows of pairs 1–2 would cost

I compared dropping both rows with the best-guess assignment below (15).
Rates assume the E3 fix, so samples with NA counts don't add trap-days.

| Plot (2012) | Samples, best guess → drop | Midges/day, best guess | Midges/day, drop | Difference (log) | Difference (SDs of the midge predictor) |
|---|---|---|---|---|---|
| hag 5 m | 6 → 5 | 155.5 | 181.7 | 0.16 | 0.09 |
| btl 50 m | 6 → 5 | 46.9 | 52.5 | 0.11 | 0.07 |
| non 50 m | 6 → 5 | 40.6 | 43.8 | 0.08 | 0.05 |
| kal 5 m | 4 → 3 | 77.8 | 41.3 | −0.63 | −0.39 |

The SD of log midge rate across all analysed plot-years (the `midges_z`
scale) is 1.65.

- **kal 5 m** already lacks its 06/19 and 07/31 samples. Dropping the pair
  leaves it with three samples (05/22, 06/05, and 07/03; 42 trap-days) and
  nothing after mid-July. Small midges were abundant at kal's other
  distances in that period (729 at 50 m, 578 at 150 m), so dropping probably
  understates its 2012 rate.
- **The other three plots** lose one of six periods each, which shifts their
  rates by under a tenth of an SD.
- **Either choice:** no plot-year rows are lost from the model; only one
  midge value per plot-year changes.

## Assigning rows from the data

Every method predicts the catch at each candidate slot from information other
than the ambiguous rows. It then compares how well the two rows fit
"row A → slot 1, row B → slot 2" against the swap. The evidence is a
log-likelihood ratio (LR), summed over three responses: big midges, small
midges, and simuliids, each as log1p(count / subsample fraction / days). The
slots are:

- **Pair 1:** hag 5 m 2012/06/05–06/19 and btl 50 m 2012/06/19–07/03.
- **Pair 2:** non 50 m 2012/07/17–07/31 and kal 5 m 2012/07/17–07/31.

Each method's stated probabilities were checked by applying the same
procedure to pairs of known samples hidden from the model. That gives the
accuracy at each strength of evidence, which is how the real pairs' results
are judged below.

### 1. Per-site multiplicative model (12)

This is a negative binomial model for 2012 only, fit per site:
count ~ distance + sampling period, with the subsample fraction as an
offset. Each distance is a constant multiple of the others across the
season.

- **Calibration was poor.** It was confidently wrong (P > 0.95 for the wrong
  assignment) for 6% of hag + btl test pairs and 24% of non + kal pairs.
- **At the real pairs' evidence strength,** accuracy on test pairs was 50%
  (pair 1) and 59% (pair 2).
- **Fit at hag 5 m:** it predicted about 2,400 small midges there for
  06/05–06/19, far above both rows, so the multiplicative assumption fits
  hag 5 m poorly in early June.

### 2. Random forest (13)

One forest per response (`ranger`), trained on all Mývatn infall samples
2008–2019. The predictors are:
- site, log distance, year, mid-sample day of year, and days out;
- the mean at the site's other distances on the same dates, and how many
  there are;
- the same plot's previous and next samples.

Predictors were rebuilt after removing each held-out sample, so its own value
never leaked into training. The predictive density is Gaussian, with SD equal
to the out-of-bag residual SD.

- On known 2012 pairs it was right 92% of the time, and confidently wrong in
  0% (hag + btl) and 0.5% (non + kal) of test pairs.
- Results were stable across 10 forest seeds.
- At hag 5 m, both pair 1 rows have far fewer big midges than expected (3.6
  and 4.0 out-of-bag SDs below the prediction), although hag 5 m's other
  2012 samples are predicted well (z within ±1.2 for big midges).

### 3. Method comparison on a larger test (14, 16)

These use 10-fold cross-validation on the five analysed sites, all years
(1,071 samples), with 20,000 random test pairs of held-out samples from
different sites in the same year.

The methods are:
- **Random forest:** as above, but trained on the five sites only.
- **Quantile regression forest:** the predictive density comes from
  per-tree draws.
- **GAM:** y ~ plot + a smooth of day of year per plot + random intercepts
  per site-date and per plot-year (`mgcv::bam`).
- **Variance-modelling GAM:** the GAM's mean, with log(squared residual)
  modelled by plot + a smooth of day of year, rescaled to the GAM's residual
  variance. A one-step `gaulss` fit failed, because the mean model has about
  as many coefficients as rows.
- **Joint scoring:** the three responses scored with a multivariate normal
  using the correlation of their cross-validation residuals.
- **Stacked classifier:** a logistic regression with no intercept on the
  per-response LRs of the random forest, GAM, and variance-modelling GAM.
  With no intercept, swapping the pair flips the prediction. It was trained
  on true-versus-swapped test pairs and evaluated leaving one year out.

Results from (16), all on the same folds and test pairs:

| Approach | Accuracy | Wrong at P > 0.95 | Mean log P(true) | Pair 1 LR (accuracy at that \|LR\|) | Pair 2 LR (accuracy at that \|LR\|) |
|---|---|---|---|---|---|
| Stacked classifier | 92.6% | 0.56% | −0.196 | −1.20 (87%) | 2.76 (95%) |
| GAM, joint | 92.1% | 0.86% | −0.211 | −1.66 (81%) | 2.60 (92%) |
| GAM, independent | 92.0% | 0.88% | −0.214 | −1.84 (81%) | 2.60 (92%) |
| Variance-modelling GAM, joint | 92.6% | 1.14% | −0.217 | −0.18 (63%) | 2.63 (92%) |
| Variance-modelling GAM, independent | 92.5% | 1.14% | −0.219 | −0.27 (63%) | 2.69 (92%) |
| Random forest, joint | 91.4% | 0.88% | −0.233 | −1.65 (85%) | 2.86 (93%) |
| Random forest, independent | 91.4% | 0.88% | −0.235 | −1.88 (85%) | 2.82 (93%) |

A negative pair 1 LR favours line 864 = hag 5 m; a positive pair 2 LR
favours line 942 = non 50 m. On a separate set of folds and pairs (14), the
quantile regression forest scored 89.5% accuracy, 1.5% confidently wrong,
and a mean log P(true) of −0.35. It was overconfident (stated 0.85,
observed 0.71), probably because of how its draws were turned into a
density, so it was dropped.

- **The GAM** is well calibrated: stated probabilities are within about 2
  percentage points of observed accuracy in every bin.
- **Joint scoring** barely matters, because residual correlations among
  responses are small (0.03–0.16).
- **Modelling the variance** doesn't help. It is slightly more accurate but
  more often confidently wrong (stated 0.97, observed 0.93).
- **The stacked classifier** is best on every summary, by small margins. It
  is slightly underconfident in the middle bins (stated 0.75, observed
  0.82). Some of its weights are negative, which reflects overlap among the
  component methods rather than anything interpretable.

## Conclusions

| Pair | Assignment favoured by all well-calibrated methods | Accuracy at that evidence strength | Does the choice matter? |
|---|---|---|---|
| 1 | line 864 = hag 5 m 06/05–06/19; line 863 = btl 50 m 06/19–07/03 | 81–87% | No: under 0.05 log units either way |
| 2 | line 942 = non 50 m 07/17–07/31; line 944 = kal 5 m 07/17–07/31 | 92–95% | Yes: 0.54 log units at non 50 m and 0.73 at kal 5 m |

The variance-modelling GAM is undecided on pair 1. The multiplicative model
favoured the opposite assignment for pair 1 but was badly calibrated.

## Limitations

- **The premise for pair 1 is doubtful.** Both rows have far fewer big
  midges than expected at hag 5 m in early June, so neither may be hag 5 m's
  sample. No method can test whether the curation note's hypothesis is
  right; every method only chooses between the two assignments it offers.
- **Test pairs are not the real situation.** They are random pairs of known
  samples, which aren't guaranteed to resemble an actual mislabelling. In
  (13), each test sample was held out alone, so its partner stayed in
  training; in the real case both slots are empty. In (14, 16), 10-fold
  cross-validation removes about 10% of samples at random, so neighbouring
  samples are usually, but not always, present.
- **Accuracy per evidence bin is an estimate.** In (16) each bin rests on
  1,500–2,400 test pairs for the bins the real pairs fall in; in (13), on
  35–57.
- **The residual correlations** for joint scoring were estimated from all
  cross-validation residuals, not re-estimated without each test pair. This
  is three parameters from about 1,000 samples, so the leakage is negligible.
- **The stacked classifier's weights** are not interpretable as the
  importance of each response or method.
- **Training data differ between scripts.** (13) trained on all Mývatn sites;
  (14, 16) on the five analysed sites. Pair results barely changed.
- **`ranger`** isn't in this project's renv library; the scripts load it from
  the global R library. Nothing was installed.

## Options

1. **Middle option:** drop both pair 1 rows, and assign pair 2 (line 942 to
   non 50 m, line 944 to kal 5 m). Dropping pair 1 costs little and avoids
   the doubtful premise; assigning pair 2 avoids losing kal 5 m's late-July
   coverage, and is supported at 92–95%.
2. **Assign both pairs** as in the conclusions table.
3. **Drop both pairs.** This makes no assumption about labels, but probably
   understates kal 5 m's 2012 rate by about half.

In every case, drop the NA rows of pairs 3 and 4.

## Decision

Both pairs were assigned as predicted (option 2). `00-clean-data-round1.R`
relabels line 863 as btl 50 m 2012/06/19–07/03 and line 944 as kal 5 m
2012/07/17–07/31, identifying rows by their values rather than position.
Lines 864 and 942 stay at hag 5 m and non 50 m. Pairs 3–4's NA copies are
dropped by the NA-count filter added for E3.

The rebuilt data match the expected values (18). The 2012 midge rates are
hag 5 m 155.5, btl 50 m 46.9, non 50 m 40.6, and kal 5 m 77.8 midges/day.
hag 5 m, btl 50 m, and non 50 m each have six samples and kal 5 m has four,
with no duplicates or overlaps.

When reporting this, say that the pair 1 assignment rests on weaker evidence
(81–87% accuracy at that evidence strength) and a doubtful premise, but
changes rates by under 0.05 log units.

## References

Breiman, L. (1996). Stacked regressions. *Machine Learning* 24:49–64.
https://doi.org/10.1007/BF00117832

Breiman, L. (2001). Random forests. *Machine Learning* 45:5–32.
https://doi.org/10.1023/A:1010933404324

Brümmer, N., & du Preez, J. (2006). Application-independent evaluation of
speaker detection. *Computer Speech & Language* 20:230–275.
https://doi.org/10.1016/j.csl.2005.08.001 (OpenAlex lists 2005, the online
date.) Logistic regression to combine log-likelihood-ratio scores; closest
match to the stacked classifier.

Platt, J. (1999). Probabilistic outputs for support vector machines and
comparisons to regularized likelihood methods. In *Advances in Large Margin
Classifiers*. MIT Press. Logistic regression to turn a score into a
calibrated probability.

Wolpert, D. H. (1992). Stacked generalization. *Neural Networks* 5:241–259.
https://doi.org/10.1016/S0893-6080(05)80023-1

Wood, S. N. (2011). Fast stable restricted maximum likelihood and marginal
likelihood estimation of semiparametric generalized linear models. *Journal
of the Royal Statistical Society: Series B* 73:3–36. (In Zotero as
`Wood2011`.)

Wright, M. N., & Ziegler, A. (2017). ranger: a fast implementation of random
forests for high dimensional data in C++ and R. *Journal of Statistical
Software* 77(1). https://doi.org/10.18637/jss.v077.i01

Yao, Y., Vehtari, A., Simpson, D., & Gelman, A. (2018). Using stacking to
average Bayesian predictive distributions (with discussion). *Bayesian
Analysis* 13:917–1007. https://doi.org/10.1214/17-BA1091

Existence and bibliographic details of these references were checked against
OpenAlex. Their contents weren't read for this note; check the relevant
sections before citing one for a specific claim.

## Scripts

| Script | Output | Content |
|---|---|---|
| `05-overlap-checks.R` | `05-overlap-checks-output.txt` | Finds the four duplicate pairs |
| `07-e5-infall-duplicates.R` | `07-e5-infall-duplicates-output.txt` | Raw rows, curation notes, missing 2012 samples |
| `12-e5-assign-pairs.R` | `12-e5-assign-pairs-output.txt` | Multiplicative model, calibration, rate sensitivity |
| `13-e5-random-forest.R` | `13-e5-random-forest-output.txt` | Random forest, calibration, seed stability, fit at hag 5 m |
| `14-e5-method-comparison.R` | `14-e5-method-comparison-output.txt` | Random forest vs. quantile forest vs. GAM, 20,000 test pairs |
| `15-e5-drop-cost.R` | `15-e5-drop-cost-output.txt` | Cost of dropping both rows of pairs 1–2 |
| `16-e5-method-comparison-2.R` | `16-e5-method-comparison-2-output.txt` | Variance-modelling GAM, joint scoring, stacked classifier |
| `18-check-e5-fix.R` | `18-check-e5-fix-output.txt` | Checks the reassignment in `00-clean-data-round1.R` |
