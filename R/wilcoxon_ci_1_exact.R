wilcoxon_ci_1_exact <- function(x, conf.int = TRUE, conf.level = 0.95) {

  A <- outer(x, x, "+") / 2
  A <- A[lower.tri(A, diag = TRUE)]

  if (!conf.int) return(list(estimate = median(A)))

  A <- sort(A)
  n <- length(x)
  N <- n * (n + 1) / 2

  c2 <- qsignrank((1 - conf.level) / 2, n)
  cl <- 1 - 2 * psignrank(c2, n)

  list(estimate   = median(A),
       conf.int   = A[c(c2 + 1, N - c2)],
       conf.level = cl)
}
