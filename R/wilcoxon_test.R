## wilcoxon_test -- unified interface for one- and two-sample Wilcoxon
## inference.  Wraps stats::wilcox.test for the test statistic and p-value,
## and uses our own validated wilcoxon_ci_* routines for the point estimate
## and confidence interval (wilcox.test's CI indexing is off by one).
##
## Returns an object of class c("wilcoxon_test", "htest").  All htest
## downstream code (print, broom::tidy, etc.) works unchanged.

wilcoxon_test <- function(x, y = NULL,
                          alternative = c("two.sided", "less", "greater"),
                          mu = 0,
                          paired = FALSE,
                          test = TRUE,
                          conf.int = TRUE,
                          conf.level = 0.95,
                          method = c("exact", "uniroot", "toms516"),
                          ...) {

  alternative <- match.arg(alternative)
  if (!missing(method)) method <- match.arg(method)

  ## Capture data name(s) for the htest object.
  xname <- deparse(substitute(x))[1L]
  yname <- if (!is.null(y)) deparse(substitute(y))[1L] else NULL
  dname <- if (is.null(y)) xname else paste(xname, "and", yname)

  ## Drop non-finite values, matching wilcox.test.
  if (is.null(y)) {
    x <- x[is.finite(x)]
  } else if (paired) {
    if (length(x) != length(y))
      stop("'x' and 'y' must have the same length")
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]
    y <- y[ok]
  } else {
    x <- x[is.finite(x)]
    y <- y[is.finite(y)]
  }

  one_sample <- is.null(y) || paired

  if (!missing(method) && one_sample && method == "toms516")
    stop('method = "toms516" is available for the two-sample problem only')

  ## Test statistic and p-value: delegate to wilcox.test when requested.
  ## When test = FALSE, we skip that call and build a minimal htest that
  ## carries just the estimation results.
  if (test) {
    result <- stats::wilcox.test(x = x, y = y,
                                 alternative = alternative,
                                 mu = mu,
                                 paired = paired,
                                 conf.int = FALSE,
                                 ...)
  } else {
    result <- list(
      method = if (one_sample)
                 "Wilcoxon Hodges-Lehmann estimation (one sample)"
               else
                 "Wilcoxon Hodges-Lehmann estimation (two sample)"
    )
    class(result) <- "htest"
  }

  ## Compute our estimate and (optionally) CI.
  effective_x <- if (paired) x - y else x

  ci_args <- list(conf.int = conf.int, conf.level = conf.level)
  if (!missing(method)) ci_args$method <- method

  if (one_sample) {
    ci_res <- do.call(wilcoxon_ci_1,
                      c(list(x = effective_x), ci_args))
    est_name <- if (paired) "(pseudo)median of x - y" else "(pseudo)median"
  } else {
    ci_res <- do.call(wilcoxon_ci_2,
                      c(list(x = x, y = y), ci_args))
    est_name <- "difference in location"
  }

  ## Assemble htest object.  Start from what we have (from wilcox.test or
  ## the minimal skeleton) and augment with our estimate / CI.
  result$data.name <- dname

  est <- ci_res$estimate
  names(est) <- est_name
  result$estimate <- est

  if (conf.int) {
    ci <- ci_res$conf.int
    attr(ci, "conf.level") <- ci_res$conf.level
    result$conf.int <- ci
  }

  class(result) <- c("wilcoxon_test", class(result))
  result
}
