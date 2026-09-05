## wilcoxon_ci_1_uniroot -- one-sample Wilcoxon CI and Hodges-Lehmann
## point estimate via root-finding on the signed-rank gradient
##
##   W(theta) = sum_i sign(x_i - theta) * R(|x_i - theta|).
##
## W(theta) is integer-valued and non-increasing, with jumps at the
## Walsh averages (x_i + x_j) / 2.  Root-finding targets are chosen so
## the step function is never at the target value, so uniroot converges
## to the correct jump point rather than into the interior of a plateau.
##
## Extra arguments in ... are forwarded to stats::uniroot; the most
## useful is tol (default .Machine$double.eps^0.25 ~ 1.2e-4).

wilcoxon_ci_1_uniroot <- function(x, conf.int = TRUE, conf.level = 0.95,
                                  ...) {

  int <- range(x)

  n <- length(x)
  N <- n * (n + 1) / 2

  ## Point estimate: median of the Walsh averages.  srgrad takes values
  ## with the same parity as N.  For odd N, srgrad = 0 is not attained
  ## and uniroot converges to A_((N+1)/2) = median.  For even N,
  ## srgrad = 0 is attained on the entire plateau [A_(N/2), A_(N/2+1)),
  ## so uniroot at target 0 returns an arbitrary point in that plateau.
  ## Target the boundaries with srgrad = +/-1 (both odd, so not attained
  ## when N is even) and average.  Symmetric with the two-sample uniroot.
  if (N %% 2L == 1L) {
    est <- uniroot(srgrad, x = x, interval = int, ...)$root
  } else {
    est_lo <- uniroot(function(theta) srgrad(theta, x) - 1,
                      interval = int, ...)$root
    est_hi <- uniroot(function(theta) srgrad(theta, x) + 1,
                      interval = int, ...)$root
    est <- (est_lo + est_hi) / 2
  }

  if (!conf.int) return(list(estimate = est))

  alpha_2 <- (1 - conf.level) / 2

  ## qsignrank/psignrank use an O(N^2) recursion (N = n(n+1)/2), so for
  ## n above ~50 they dominate wall time.  Switch to a normal
  ## approximation with continuity correction beyond the cutoff.
  ## Cutoff is a tuning parameter; may need adjustment.
  n_exact_cut <- 50

  if (n <= n_exact_cut) {
    c2 <- qsignrank(alpha_2, n)
    cl <- 1 - 2 * psignrank(c2, n)
  } else {
    mu  <- n * (n + 1) / 4
    sig <- sqrt(n * (n + 1) * (2 * n + 1) / 24)
    z   <- qnorm(alpha_2)                       # negative
    c2  <- ceiling(mu + z * sig - 0.5)
    cl  <- 1 - 2 * pnorm((c2 + 0.5 - mu) / sig)
  }

  target <- N - 2 * c2 - 1
  lo <- uniroot(function(theta) srgrad(theta, x) - target,
                interval = int, ...)$root
  hi <- uniroot(function(theta) srgrad(theta, x) + target,
                interval = int, ...)$root

  list(estimate   = est,
       conf.int   = c(lo, hi),
       conf.level = cl)
}
