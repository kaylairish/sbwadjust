# sbwadjust

Core routines behind *"Simple Covariate Adjustment for Many Estimands Using
Stable Balancing Weights"* (Irish, Zubizarreta, Luedtke;
[arXiv:2609.01638](https://arxiv.org/abs/2609.01638)): fitting stable balancing
weights (SBW) for a two-arm study, and weighted Kaplan–Meier survival-ratio
estimation with bootstrap confidence intervals.

> **Scope.** This is currently the paper's method code, not a general-purpose
> covariate-adjustment package. The functions expose the exact choices made in
> the paper — exact balance (imbalance tolerance fixed at zero), balance of the
> supplied covariate columns to the pooled-sample mean, 0/1 treatment. A
> general-purpose interface (user-set tolerance, balance functions as an
> argument, an estimand-agnostic weighting workflow) is in development.

The simulation studies and data application that build on these routines live in
[`sbw-covariate-adjustment-code`](https://github.com/kaylairish/sbw-covariate-adjustment-code).

## Installation

```r
# install.packages("remotes")
remotes::install_github("kaylairish/sbwadjust")
```

## What's in it

| File | Functions |
|---|---|
| `R/weights.R` | `get_weights_for_group_neg`, `get_weights_for_group_nonneg`, `get_sbws_for_study` — fit SBW for one arm or a full two-arm study: a closed-form solve first, falling back to a nonnegative quadratic program when the closed-form weights go negative. |
| `R/km_ratio.R` | `km_ratio_loglog_greenwood`, `boot_km_ratio` — weighted Kaplan–Meier survival-ratio point estimates and bootstrap CIs (Wald or percentile). |
| `R/ipw.R` | `get_ipws_for_study` — inverse-probability weights, used as the `weight_type = "IPW"` comparison in `boot_km_ratio`. |

## Tests

```r
# from a local clone
devtools::test()
```

Each exported function has a regression test pinned to a fixed seed / toy input
(`tests/testthat/`).
