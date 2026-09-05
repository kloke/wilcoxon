## wilcoxon_ci_2 -- two-sample Wilcoxon confidence interval and point
## estimate for the shift parameter Delta in the model X = Y + Delta.
## Positive Delta means X shifted right of Y (stats::wilcox.test
## convention).  Interface to three backends:
##
##   "exact"    -- Bauer's direct construction: sort the mn differences
##                 x_i - y_j, pick order statistics indexed by qwilcox.
##                 O(mn) memory, O(mn log(mn)) time.  Reference
##                 implementation.
##
##   "uniroot"  -- Invert the Mann-Whitney U statistic by root-finding
##                 with stats::uniroot.  O(m + n) memory.
##
##   "toms516"  -- Ryan and McKean (1977), TOMS Algorithm 516.  Illinois-
##                 modified regula falsi with a two-pointer evaluator for
##                 U.  O(m + n) memory.
##
## When method is not specified, auto-selects "exact" if max(m, n) < 50
## and c(x, y) has no ties, else "uniroot".  Same rule as stats::wilcox.test
## uses to decide between exact and asymptotic p-values.  toms516 is never
## selected automatically; request it by name.
##
## Extra arguments in ... are forwarded to the uniroot backend (see
## wilcoxon_ci_2_uniroot); ignored by the exact and toms516 backends.
##
## All backends return list(estimate, conf.int, conf.level).

wilcoxon_ci_2 <- function(x, y, conf.int = TRUE, conf.level = 0.95,
                          method = c("exact", "uniroot", "toms516"), ...) {

  if (missing(method)) {
    method <- if (max(length(x), length(y)) < 50 && !anyDuplicated(c(x, y)))
                "exact"
              else
                "uniroot"
  } else {
    method <- match.arg(method)
  }

  switch(method,
    exact   = wilcoxon_ci_2_exact  (x, y, conf.int = conf.int,
                                    conf.level = conf.level),
    uniroot = wilcoxon_ci_2_uniroot(x, y, conf.int = conf.int,
                                    conf.level = conf.level, ...),
    toms516 = wilcoxon_ci_2_toms516(x, y, conf.int = conf.int,
                                    conf.level = conf.level))
}
