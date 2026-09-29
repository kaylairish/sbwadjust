# sbwadjust 0.3.0

* Removed `get_ipws_for_study()` and the `weight_type` / `ipw_use_glm`
  arguments of `boot_km_ratio()`, which now always uses SBW weights. The
  package is scoped to SBW; inverse-probability weighting was only a
  comparison method for the paper's simulations. Code that needs it can
  install v0.2.0 (`remotes::install_github("kaylairish/sbwadjust@v0.2.0")`).
* Renamed `km_ratio_loglog_greenwood()` to `km_ratio_greenwood()`. Its CI
  was always a Wald CI on the log scale; the old name described an
  intermediate log-log step that cancels out. That step is gone:
  `se(log S)` is now computed directly as `se(S) / S`, which gives identical
  results except that an arm with no events by `t0` (S = 1) now contributes
  zero to `se_log` instead of making it `NaN`.
* `boot_km_ratio()`: dropped the unused `verbose` argument; the SBW clipping
  summaries are now `NA` (not `-Inf`/`NaN` with a warning) when no bootstrap
  resample's SBW fit succeeds; the help page now documents that a failed
  resample uses the unadjusted KM ratio.
* `km_ratio_greenwood()` now accepts a factor treatment indicator, and its
  help page notes that the SE treats the weights as fixed.
* `sbw_estimate()` now stops with a clear message when the outcome has a
  different length from the data the weights were fit on (e.g. an outcome
  variable missing from `data` that R found elsewhere), instead of failing
  later with an unrelated bootstrap error.
* Fixed: `sbw_estimate()` returned an `NaN` standard error and CI when some
  bootstrap resamples gave an infinite estimate (e.g. `"RR"` with a rare
  outcome, where a resample can draw no control events). Those resamples now
  count as failed and are reported in `boot_fail_rate`.
* `sbw_estimate()` now gives a clear error when a `Surv()` outcome is used
  with an estimand other than `"survival_ratio"`.
* `sbw_weights()` now accepts the treatment column name as a string
  (`treatment = "arm"`, or a variable holding it), not only unquoted.
  Previously a quoted name failed with a misleading "exactly two levels"
  error.
* `sbw_weights()` now gives clear errors for a treatment with missing values
  or with only one arm, which previously failed inside the solver with
  unrelated messages (e.g. "system is exactly singular").
* `?sbw_weights` now states which level of a factor or character treatment
  is treated (the second level; alphabetical for character), and
  `print.sbw_fit()` names the treated and control levels in that case.
* `sbw_weights()` now gives clear errors when the balance covariates are
  collinear within an arm (including a factor level that never occurs in one
  arm) or when exact balance with nonnegative weights is infeasible, instead
  of the solver's "system is exactly singular" or quadprog's "constraints are
  inconsistent, no solution!".
* `print.sbw_fit()` now reports how many units in each arm got weight 0,
  replacing the "weight(s) clipped at 0" note, which counted rounding noise
  in the solver rather than dropped units. `summary()` does the same: its
  `n_clipped` / `max_abs_clipped` elements are replaced by `n_zero` (by
  arm), and its balance table prints to 4 significant digits.
* `summary()` no longer reports an effective sample size. The Kish ESS
  measures how concentrated the weights are, not the precision of the
  treatment-effect estimate, and was easy to misread as the latter.

# sbwadjust 0.2.0

"Usable by a stranger" release: a user-facing formula API on top of the v0.1 core.

* New `sbw_weights()`: formula/data/treatment front end to the core SBW solver,
  returning an `sbw_fit` object with `print()`, `summary()`, `plot()`, and
  `weights()` methods.
* New `sbw_estimate()`: treatment-effect estimation from an `sbw_fit`, with a
  bootstrap confidence interval, for a closed menu of estimands: average
  treatment effect (`"ATE"`), relative risk (`"RR"`), survival ratio
  (`"survival_ratio"`), Mann-Whitney win probability (`"mann_whitney"`,
  uncensored outcomes only), and quantile contrasts (`"quantile_diff"` /
  `"quantile_ratio"`).
* Quantile contrasts use a new `.weighted_quantile()` internal helper: the
  generalized-inverse (step-function / type 1) SBW-weighted empirical
  quantile, matching the empirical-process framework the paper's
  differentiability results use.
* `sbw_estimate(estimand = "survival_ratio")` now warns when the SBW point
  estimate or bootstrap SE is non-finite and the unadjusted Kaplan-Meier
  ratio is returned in its place (previously flagged only by `mc_fail`).
* `sbw_estimate()` gains a `data` argument, so weights can be fit on
  baseline data before outcomes exist and outcomes supplied at analysis
  time. Rows must line up one-to-one with the fitted data; row count and any
  shared columns are checked.
* Fixed: `sbw_estimate()`'s bootstrap failed on every resample when
  `sbw_weights()` was given `treatment` as a vector rather than a column
  name. The bootstrap now resamples the stored 0/1 treatment directly.
* `Surv()` in a `survival_ratio` outcome formula now resolves without
  attaching the survival package, and missing survival times are rejected
  instead of silently dropped.
* `?sbw_estimate` now notes that the row bootstrap assumes simple
  randomization and does not account for stratified or covariate-adaptive
  designs.
* CRAN Repository Policy compliance pass: fixed a roxygen markdown escaping
  bug (`\%in\%` rendering with stray backslashes in `?sbw_estimate`), added a
  runnable `@examples` block to `sbw_estimate()`, fixed the `DESCRIPTION`
  citation style, added `URL`/`BugReports` fields, and added `inst/WORDLIST`
  for `spelling::spell_check()`.

# sbwadjust 0.1.0

Initial public release. Core stable-balancing-weight machinery, consolidated
from the four near-duplicate copies used across the paper's simulations:

* `get_weights_for_group_neg()`, `get_weights_for_group_nonneg()`,
  `get_sbws_for_study()` — closed-form SBW solve with a nonnegative
  quadratic-programming fallback.
* `km_ratio_loglog_greenwood()`, `boot_km_ratio()` — weighted Kaplan-Meier
  survival-ratio point estimates and bootstrap confidence intervals.
* `get_ipws_for_study()` — inverse-probability weights, used as a comparison
  method.
