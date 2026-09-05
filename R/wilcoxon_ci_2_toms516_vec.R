## wilcoxon_ci_2_toms516_vec -- two-sample Wilcoxon CI for the shift
## parameter Delta in the model X = Y + Delta, following Ryan and McKean
## (1977), TOMS Algorithm 516.
##
## Identical algorithm to wilcoxon_ci_2_toms516_port; the only change is
## that the Mann-Whitney U evaluator is fmann_vec (findInterval-based,
## walks the sorted arrays with a compiled C loop) rather than fmann
## (two-pointer walk in interpreted R).  Same asymptotic O(m + n) cost
## per evaluation but a large constant-factor speedup for m, n above a
## few dozen.
##
## See wilcoxon_ci_2_toms516_port for the algorithm's derivation and for
## the deviations from the original Fortran.

wilcoxon_ci_2_toms516_vec <- function(x, y, conf.int = TRUE,
                                      conf.level = 0.95) {

  m <- length(x)
  n <- length(y)
  mn <- m * n

  xs <- sort(x)
  ys <- sort(y)

  lo <- xs[1L] - ys[n]
  hi <- xs[m] - ys[1L]

  tol <- .Machine$double.eps^0.5 * max(1, hi - lo)

  U <- function(d) fmann_vec(d, ys, xs)

  Flo <- mn
  Fhi <- 0

  root_at <- function(target) {
    ill(target, lo, Flo, hi, Fhi, U, tol)
  }

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
