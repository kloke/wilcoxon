## wilcoxon_ci_2_uniroot -- two-sample Wilcoxon CI for the shift parameter
## Delta in the model X = Y + Delta, by inverting the Wilcoxon rank-sum
## statistic with stats::uniroot.
##
## Model: X_i iid F(. - Delta), Y_j iid F.  Positive Delta means X shifted
## right of Y (stats::wilcox.test convention).
##
## W(Delta) = sum of ranks of (x_i - Delta) in the combined sample
##            { x_1 - Delta, ..., x_m - Delta, y_1, ..., y_n }.
##
## W(Delta) is a non-increasing step function of Delta with jumps at the
## mn pairwise differences x_i - y_j.  Related to Mann-Whitney U by
## U = W - m(m+1)/2, so R's qwilcox / pwilcox (on U) give the critical
## values and achieved level directly.
##
## Evaluator is O((m+n) log(m+n)); no mn matrix is materialized.
##
## Root-finding targets are half-integers.  U is integer-valued and
## constant on each plateau [D_(k), D_(k+1)) between adjacent sorted
## differences, so an integer target would place uniroot's zero on an
## entire plateau; uniroot would converge to an arbitrary point in that
## plateau rather than to the CI bound (a d_{ij}).  A half-integer target
## places the sign change of W - target exactly at a jump of W.
##
## Extra arguments in ... are forwarded to stats::uniroot; the most
## useful is tol (default .Machine$double.eps^0.25 ~ 1.2e-4).

wilcoxon_ci_2_uniroot <- function(x, y, conf.int = TRUE, conf.level = 0.95,
                                  ...) {

  m <- length(x)
  n <- length(y)
  mn <- m * n
  xidx <- seq_len(m)

  W <- function(delta) sum(rank_(c(x - delta, y))[xidx])

  interval <- c(min(x) - max(y), max(x) - min(y))

  ## W = U + m(m+1)/2.
  offset <- m * (m + 1L) / 2L
  root_at <- function(U_target) {
    W_target <- U_target + offset
    uniroot(function(d) W(d) - W_target, interval = interval, ...)$root
  }

  ## Point estimate: median of the differences.  Odd mn -> single half-
  ## integer target at mn/2.  Even mn -> average of the roots at
  ## (mn/2) +/- 0.5.
  if (mn %% 2L == 1L) {
    estimate <- root_at(mn / 2)
  } else {
    estimate <- (root_at(mn / 2 - 0.5) + root_at(mn / 2 + 0.5)) / 2
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
