## wilcoxon_ci_2_toms516_port -- two-sample Wilcoxon CI for the shift parameter
## Delta in the model X = Y + Delta, following Ryan and McKean (1977),
## TOMS Algorithm 516.
##
## Model: X_i iid F(. - Delta), Y_j iid F.  Positive Delta means X shifted
## right of Y (stats::wilcox.test convention).
##
## Inverts the Mann-Whitney U statistic
##
##   U(Delta) = #{ (i, j) : x_i - y_j > Delta }
##
## by root-finding, using the Illinois modification of regula falsi (Dowell
## and Jarratt, 1971).  Two elements make the algorithm fast:
##
##   1. fmann evaluates U(Delta) in O(m + n) once x and y are sorted, by
##      walking a single pointer through one sorted array as it sweeps
##      through the other.
##
##   2. Root-finding targets are half-integers.  U is integer-valued with
##      jumps at the pairwise differences x_i - y_j, so a target like
##      c + 0.5 is never attained; the iteration converges into the gap
##      between two adjacent d_{ij}'s.
##
## Deviations from the original Fortran, for the record:
##
##   * TMEAN starting values and the associated BRACK expansion step are
##     omitted.  We start ill with the bracket c(min(x) - max(y),
##     max(x) - min(y)), matching wilcoxon_ci_2_uniroot.
##
##   * The critical value comes from qwilcox rather than from a normal
##     approximation with continuity correction.  This also removes TOMS's
##     post-hoc adjustment to the requested confidence percent.
##
##   * Double precision throughout; tolerances are set from
##     .Machine$double.eps rather than tuned to 1970s single-precision
##     Fortran.
##
##   * Error checking is left to R's usual mechanisms rather than a CHECK
##     routine.

wilcoxon_ci_2_toms516_port <- function(x, y, conf.int = TRUE, conf.level = 0.95) {

  m <- length(x)
  n <- length(y)
  mn <- m * n

  xs <- sort(x)
  ys <- sort(y)

  lo <- xs[1L] - ys[n]
  hi <- xs[m] - ys[1L]

  tol <- .Machine$double.eps^0.5 * max(1, hi - lo)

  ## fmann(delta, ys, xs) computes #{ (i, j) : xs_i > ys_j + delta }
  ## = #{ (i, j) : xs_i - ys_j > delta } = U(delta).
  U <- function(d) fmann(d, ys, xs)

  Flo <- mn                                 # U(lo) = mn
  Fhi <- 0                                  # U(hi) = 0

  root_at <- function(target) {
    ill(target, lo, Flo, hi, Fhi, U, tol)
  }

  ## Point estimate: median of the differences.  Odd mn -> single half-
  ## integer target at mn/2.  Even mn -> average of the roots at
  ## (mn/2) +/- 0.5.
  if (mn %% 2L == 1L) {
    estimate <- root_at(mn / 2)
  } else {
    est_lo <- root_at(mn / 2 + 0.5)
    est_hi <- root_at(mn / 2 - 0.5)
    estimate <- (est_lo + est_hi) / 2
  }

  if (!conf.int) return(list(estimate = estimate))

  alpha2 <- (1 - conf.level) / 2
  c <- qwilcox(alpha2, m, n)
  achieved <- 1 - 2 * pwilcox(c, m, n)

  ci_lo <- root_at(mn - c - 0.5)
  ci_hi <- root_at(c + 0.5)

  list(estimate   = estimate,
       conf.int   = c(ci_lo, ci_hi),
       conf.level = achieved)
}
