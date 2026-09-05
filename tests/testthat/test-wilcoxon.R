## Small starter test set for wilcoxon.
##
## Covers:
##   1. Backends agree with each other (approximately).
##   2. Point estimate matches wilcox.test's HL estimate.
##   4. Backend return-structure contracts.
##   5. wilcoxon_test(test = FALSE) preserves the htest contract.
##   7. Ties in the input do not error.

set.seed(1)
x1 <- rnorm(20)
y1 <- rnorm(15, mean = 0.5)

## Tolerance for root-finder vs exact.  uniroot's default absolute tolerance
## is .Machine$double.eps^0.25 ~ 1.2e-4; toms516's is ~1.5e-8 * (bracket
## width).  1e-3 is comfortably above both.
tol <- 1e-3

## ---- 1. Backends agree with each other -----------------------------------

test_that("one-sample backends agree", {
  ex <- wilcoxon_ci_1_exact  (x1)
  un <- wilcoxon_ci_1_uniroot(x1)
  expect_equal(un$estimate,   ex$estimate,   tolerance = tol)
  expect_equal(un$conf.int,   ex$conf.int,   tolerance = tol)
  expect_equal(un$conf.level, ex$conf.level)
})

test_that("two-sample backends agree", {
  ex <- wilcoxon_ci_2_exact  (x1, y1)
  un <- wilcoxon_ci_2_uniroot(x1, y1)
  to <- wilcoxon_ci_2_toms516(x1, y1)
  expect_equal(un$estimate,   ex$estimate,   tolerance = tol)
  expect_equal(to$estimate,   ex$estimate,   tolerance = tol)
  expect_equal(un$conf.int,   ex$conf.int,   tolerance = tol)
  expect_equal(to$conf.int,   ex$conf.int,   tolerance = tol)
  expect_equal(un$conf.level, ex$conf.level)
  expect_equal(to$conf.level, ex$conf.level)
})

## ---- 2. Point estimate matches wilcox.test -------------------------------

test_that("one-sample HL matches wilcox.test", {
  wt <- suppressWarnings(stats::wilcox.test(x1, conf.int = TRUE))
  ex <- wilcoxon_ci_1_exact(x1)
  expect_equal(unname(wt$estimate), ex$estimate, tolerance = tol)
})

test_that("two-sample HL matches wilcox.test", {
  wt <- suppressWarnings(stats::wilcox.test(x1, y1, conf.int = TRUE))
  ex <- wilcoxon_ci_2_exact(x1, y1)
  expect_equal(unname(wt$estimate), ex$estimate, tolerance = tol)
})

## ---- 4. Return-structure contracts ---------------------------------------

test_that("backends return list(estimate, conf.int, conf.level) with conf.int = TRUE", {
  for (f in list(wilcoxon_ci_1_exact, wilcoxon_ci_1_uniroot)) {
    r <- f(x1)
    expect_named(r, c("estimate", "conf.int", "conf.level"))
    expect_length(r$conf.int, 2L)
  }
  for (f in list(wilcoxon_ci_2_exact, wilcoxon_ci_2_uniroot, wilcoxon_ci_2_toms516)) {
    r <- f(x1, y1)
    expect_named(r, c("estimate", "conf.int", "conf.level"))
    expect_length(r$conf.int, 2L)
  }
})

test_that("backends return list(estimate) with conf.int = FALSE", {
  for (f in list(wilcoxon_ci_1_exact, wilcoxon_ci_1_uniroot)) {
    r <- f(x1, conf.int = FALSE)
    expect_named(r, "estimate")
  }
  for (f in list(wilcoxon_ci_2_exact, wilcoxon_ci_2_uniroot, wilcoxon_ci_2_toms516)) {
    r <- f(x1, y1, conf.int = FALSE)
    expect_named(r, "estimate")
  }
})

## ---- 5. htest contract preserved when test = FALSE -----------------------

test_that("wilcoxon_test(test = FALSE) preserves the htest contract", {
  res <- wilcoxon_test(x1, y1, test = FALSE)
  expect_s3_class(res, "wilcoxon_test")
  expect_s3_class(res, "htest")
  expect_null(res$statistic)
  expect_null(res$p.value)
  expect_null(res$parameter)
  expect_null(res$null.value)
  expect_false(is.null(res$estimate))
  expect_false(is.null(res$data.name))
  expect_false(is.null(res$method))
  ## print.htest must not error on the reduced object.
  expect_output(print(res), "estimation")
})

## ---- 7. Ties in input do not error ---------------------------------------
##
## Note: exact discrete critical values assume no ties, so the reported
## achieved level under ties is an approximation.  These tests just guard
## against errors and check that the numeric output is finite.

test_that("ties in one-sample input do not error", {
  x_ties <- c(1, 2, 2, 3, 4, 4, 5, 6, 7, 8)
  res <- wilcoxon_ci_1(x_ties)
  expect_true(is.finite(res$estimate))
  expect_true(all(is.finite(res$conf.int)))
})

test_that("ties in two-sample input do not error", {
  set.seed(2)
  x_ties <- c(rnorm(10), 0, 0)
  y_ties <- c(rnorm(10), 0, 0)
  res <- wilcoxon_ci_2(x_ties, y_ties)
  expect_true(is.finite(res$estimate))
  expect_true(all(is.finite(res$conf.int)))
})
