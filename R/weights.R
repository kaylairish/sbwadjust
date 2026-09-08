# Shared stable-balancing-weight (SBW) solver used by the binary-ATE,
# survival-risk-ratio, coherence, and illustration-of-approach (Table 1)
# simulations. Consolidated 2026-07-20 (CODE_PACKAGE_PLAN.md Step 1) from
# four near-duplicate copies; see that file's "Step 0 findings" for the
# audit this is based on.
#
# Return contract: get_weights_for_group_nonneg / get_sbws_for_study both
# return a list (`w`, `n_clipped`, `max_abs_clipped`), not a plain vector —
# this matches what survival-risk-ratio and table_1.R already did; binary-ate
# and coherence's call sites were updated to unwrap `$w`.
#
# Tolerance: get_sbws_for_study accepts the closed-form solve as long as it's
# nonnegative within 1e-10 (zeroing out negligible negative floating-point
# dust) rather than requiring exact nonnegativity, before falling back to the
# slower QP solve. The closed-form solve occasionally
# returns weights that are negative only by numerical noise.

#' Closed-form SBW group weights (negative weights allowed)
#'
#' Solves for weights on one treatment-arm's covariate data that exactly
#' balance its mean to the target `X_n`, minimizing sum of squared weights,
#' via the closed-form solution (no nonnegativity constraint).
#'
#' @param data Covariate data frame for the group (one treatment arm).
#' @param X_n Target covariate mean to balance to.
#' @return Numeric vector of weights, one per row of `data`.
#' @export
get_weights_for_group_neg = function(data, X_n) {
  Na = nrow(data)
  X  = cbind(1, data.matrix(data))
  Xn = c(1, X_n)
  XtX = crossprod(X)
  w = rep(1, Na) - X %*% solve(XtX, t(X) %*% rep(1, Na) - Na * Xn)
  as.numeric(w)
}

#' Clip small/negative numerical noise out of a QP weight solution
#'
#' @param w Numeric vector of weights (possibly containing negative
#'   floating-point dust from the QP solve).
#' @return A list with `w` (clipped weights), `n_clipped` (how many entries
#'   were negative before clipping), and `max_abs_clipped` (largest
#'   magnitude among the clipped entries, 0 if none).
#' @keywords internal
.clip_weights = function(w) {
  w = as.numeric(w)
  idx = which(w < 0)
  n_clipped = length(idx)
  max_abs_clipped = if (n_clipped == 0L) 0 else max(abs(w[idx]))
  if (n_clipped > 0L) w[idx] = 0
  # also kill tiny negative numerical dust
  w[w < 0] = 0

  list(
    w = w,
    n_clipped = n_clipped,
    max_abs_clipped = max_abs_clipped
  )
}

#' Nonnegative-QP SBW group weights, with clipping diagnostics
#'
#' Solves the same balance objective as [get_weights_for_group_neg()], but
#' constrained to nonnegative weights via quadratic programming.
#'
#' @param data Covariate data frame for the group (one treatment arm).
#' @param X_n Target covariate mean to balance to.
#' @return A list with `w` (numeric vector of weights), `n_clipped`, and
#'   `max_abs_clipped` (see [.clip_weights()]).
#' @export
get_weights_for_group_nonneg = function(data, X_n) {
  Na = nrow(data)
  Xa = cbind(1, data.matrix(data))
  Xnb = c(1, X_n)

  Dmat = 2 * diag(Na)
  dvec = 2 * rep(1, Na)

  Aeq  = Xa
  beq  = Na * Xnb

  Aineq = diag(Na)
  bineq = rep(0, Na)

  Amat = cbind(Aeq, Aineq)
  bvec = c(beq, bineq)

  sol = quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = ncol(Aeq))
  .clip_weights(sol$solution)
}

#' Study-level SBW weights: closed-form first, QP fallback if needed
#'
#' Tries the unconstrained closed-form solve first for each arm
#' ([get_weights_for_group_neg()]); if the result is nonnegative within
#' floating-point tolerance, zeroes out negligible negative dust and returns
#' it directly. Otherwise falls back to the nonnegative-QP solve
#' ([get_weights_for_group_nonneg()]).
#'
#' @param X_subset Covariate data frame for the full study sample.
#' @param A Treatment indicator (0/1) for the full study sample.
#' @return A list with `w` (numeric vector of weights, one per row of
#'   `X_subset`), `n_clipped`, and `max_abs_clipped` (both 0 if the
#'   closed-form solve was used).
#' @export
get_sbws_for_study = function(X_subset, A) {
  A = as.integer(A)
  treated = X_subset[A == 1, , drop = FALSE]
  control = X_subset[A == 0, , drop = FALSE]

  X_n = colMeans(rbind(treated, control))

  w = numeric(length(A))

  # try closed-form first
  w[A == 1] = get_weights_for_group_neg(treated, X_n)
  w[A == 0] = get_weights_for_group_neg(control, X_n)

  # accept closed-form if nonnegative within floating-point tolerance,
  # zeroing out negligible negative dust
  if (all(w >= -1e-10)) {
    w[w < 0] = 0
    return(list(w = w, n_clipped = 0L, max_abs_clipped = 0))
  }

  # fallback to nonneg QP (with clipping diagnostics)
  tr = get_weights_for_group_nonneg(treated, X_n)
  ct = get_weights_for_group_nonneg(control, X_n)

  w[A == 1] = tr$w
  w[A == 0] = ct$w

  list(
    w = w,
    n_clipped = tr$n_clipped + ct$n_clipped,
    max_abs_clipped = max(tr$max_abs_clipped, ct$max_abs_clipped)
  )
}
