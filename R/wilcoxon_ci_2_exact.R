## wilcoxon_ci_2_exact -- exact two-sample Wilcoxon CI for the shift
## parameter Delta in the model X = Y + Delta, via Bauer's direct
## construction (JASA 1972).
##
## For any rank statistic S in the two-sample location problem, S(Delta) is
## a step function whose jumps occur only at the mn pairwise differences
## d_{ij} = x_i - y_j.  For the Wilcoxon statistic every jump has amount -1,
## so the CI reduces to picking order statistics of the differences, indexed
## by the critical value of the Mann-Whitney null distribution.
##
## Model: X_i iid F(. - Delta), Y_j iid F, independent.  Positive Delta
## means X is shifted to the right of Y.  This is the convention of
## stats::wilcox.test.  Point estimate is the Hodges-Lehmann estimator,
## med { x_i - y_j }.
##
## Returned confidence level is the *achieved* level 1 - 2 * pwilcox(c, m, n),
## which for exact discrete distributions is generally not equal to conf.level.

wilcoxon_ci_2_exact <- function(x, y, conf.int = TRUE, conf.level = 0.95) {

  m <- length(x)
  n <- length(y)

  ## All mn pairwise differences x_i - y_j, sorted.
  D <- sort(as.vector(outer(x, y, "-")))

  estimate <- median(D)

  if (!conf.int) return(list(estimate = estimate))

  alpha2 <- (1 - conf.level) / 2
  c <- qwilcox(alpha2, m, n)
  achieved <- 1 - 2 * pwilcox(c, m, n)

  list(estimate   = estimate,
       conf.int   = D[c(c + 1L, m * n - c)],
       conf.level = achieved)
}
