## wilcoxon_ci_2_toms516 -- dispatcher for the toms516 family.
##
##   "vec"   -- current default.  Uses fmann_vec (findInterval, compiled
##              inner loop) for the U evaluator.  Same algorithm as the
##              port; typically much faster for m, n above a few dozen.
##
##   "port"  -- faithful port of the original TOMS 516 Fortran, with the
##              two-pointer FMANN walk written in R.  Kept for teaching
##              and comparison; slower than "vec" in typical use.
##
## Both variants return list(estimate, conf.int, conf.level) with the
## same convention as the other backends.  See wilcoxon_ci_2 for the
## outer dispatcher.

wilcoxon_ci_2_toms516 <- function(x, y, conf.int = TRUE, conf.level = 0.95,
                                  variant = c("vec", "port")) {

  variant <- match.arg(variant)

  switch(variant,
    vec  = wilcoxon_ci_2_toms516_vec (x, y, conf.int = conf.int,
                                      conf.level = conf.level),
    port = wilcoxon_ci_2_toms516_port(x, y, conf.int = conf.int,
                                      conf.level = conf.level))
}
