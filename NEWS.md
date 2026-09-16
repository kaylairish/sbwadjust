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
