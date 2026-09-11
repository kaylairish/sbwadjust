# sbwadjust

Covariate adjustment for randomized trials using stable balancing weights
(SBW), from *"Simple Covariate Adjustment for Many Estimands Using Stable
Balancing Weights"* (Irish, Zubizarreta, Luedtke;
[arXiv:2609.01638](https://arxiv.org/abs/2609.01638)).

```r
library(sbwadjust)

## 1. Design stage -- before unblinding. Outcome-blind, prespecifiable.
sbw <- sbw_weights(
  balance   = ~ age + sex + bmi + region,   # covariates to balance
  data      = trial_baseline,
  treatment = arm
)
summary(sbw)          # balance table + effective-sample-size diagnostics

## 2. Analysis stage -- after unblinding.
sbw_estimate(sbw, Y ~ 1, estimand = "RR")                             # relative risk
sbw_estimate(sbw, Surv(time, status) ~ 1,
             estimand = "survival_ratio", horizon = 52)

## An estimand we don't cover? Take the weights, use your own estimator.
w <- weights(sbw)
```

`sbw_weights()` fits exact-balance SBW (imbalance tolerance fixed at zero,
matching the paper) via a closed-form solve with a nonnegative
quadratic-programming fallback. `sbw_estimate()` computes one of a closed
menu of estimands from the fitted weights, with a bootstrap confidence
interval: average treatment effect (`"ATE"`), relative risk (`"RR"`),
survival ratio (`"survival_ratio"`, needs a `horizon` argument), Mann-Whitney
win probability (`"mann_whitney"`, finite/uncensored outcomes), and quantile
contrasts (`"quantile_diff"` / `"quantile_ratio"`, with a `probs` argument).
No user-supplied functional is accepted — an uncovered estimand means: take
`weights(sbw)` and run (and bootstrap) your own estimator.

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
| `R/sbw_weights.R` | `sbw_weights` — formula/data/treatment front end to `get_sbws_for_study`, returning an `sbw_fit` object with `print`, `summary`, `plot`, and `weights` methods. |
| `R/sbw_estimate.R` | `sbw_estimate` — treatment-effect estimation from an `sbw_fit`: ATE, RR, survival ratio, Mann-Whitney, and quantile contrasts, each with a bootstrap CI. |
| `R/weights.R` | `get_weights_for_group_neg`, `get_weights_for_group_nonneg`, `get_sbws_for_study` — fit SBW for one arm or a full two-arm study: a closed-form solve first, falling back to a nonnegative quadratic program when the closed-form weights go negative. |
| `R/km_ratio.R` | `km_ratio_loglog_greenwood`, `boot_km_ratio` — weighted Kaplan–Meier survival-ratio point estimates and bootstrap CIs (Wald or percentile); `sbw_estimate(..., estimand = "survival_ratio")` wraps `boot_km_ratio`. |
| `R/ipw.R` | `get_ipws_for_study` — inverse-probability weights, used as the `weight_type = "IPW"` comparison in `boot_km_ratio`. |

> **Estimand coverage note.** ATE, RR, and survival ratio have full worked
> estimators and simulations in the paper. `mann_whitney` implements only the
> finite/uncensored case given in the JASA supplement (right-censored
> outcomes are future work). Quantile contrasts use the natural SBW-weighted
> empirical-quantile plug-in (a generalized-inverse weighted quantile,
> arm-specific difference/ratio, bootstrap CI); this recipe isn't spelled out
> verbatim in the paper, unlike the other estimands. RMST is not included in
> this release.

## Tests

```r
# from a local clone
devtools::test()
```

Each exported function has a regression test pinned to a fixed seed / toy input
(`tests/testthat/`).
