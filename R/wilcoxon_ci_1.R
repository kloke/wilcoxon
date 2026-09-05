## wilcoxon_ci_1 -- one-sample Wilcoxon CI and Hodges-Lehmann point
## estimate.  Interface to two backends: exact construction via Walsh
## averages, and root-finding on the signed-rank gradient via uniroot.
##
## When method is not specified, auto-selects "exact" if length(x) < 50
## and x has no ties, else "uniroot".  Same rule as stats::wilcox.test
## uses to decide between exact and asymptotic p-values.  Users who call
## wilcoxon_ci_1(x, method = "exact") get "exact" regardless of n.
##
## Extra arguments in ... are forwarded to the uniroot backend (see
## wilcoxon_ci_1_uniroot); ignored by the exact backend.

wilcoxon_ci_1 <- function(x, conf.int = TRUE, conf.level = 0.95,
                          method = c("exact", "uniroot"), ...) {

  if (missing(method)) {
    method <- if (length(x) < 50 && !anyDuplicated(x))
                "exact"
              else
                "uniroot"
  } else {
    method <- match.arg(method)
  }

  switch(method,
    exact   = wilcoxon_ci_1_exact  (x, conf.int = conf.int,
                                    conf.level = conf.level),
    uniroot = wilcoxon_ci_1_uniroot(x, conf.int = conf.int,
                                    conf.level = conf.level, ...))
}
