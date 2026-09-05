## Internal helpers, not exported.

## Fast rank without the argument checking of base::rank.  Copied under
## GPL from base R, following the pattern in vig.Rmd.
rank_ <- function(x) sort.list(sort.list(x))

## Signed-rank "gradient".  W(theta) = sum_i sign(x_i - theta) * R(|x_i - theta|).
## Integer-valued, non-increasing step function of theta with jumps at the
## Walsh averages (x_i + x_j) / 2.
srgrad <- function(theta, x) {
  e <- x - theta
  drop(crossprod(rank_(abs(e)), sign(e)))
}

## fmann: Mann-Whitney U evaluator.  Two-pointer walk on sorted arrays,
## faithful to the FMANN subroutine of TOMS 516.
##
## Given sorted xs, sorted ys, and shift delta, returns
##
##   U = #{ (i, j) : y_j > x_i + delta }.
##
## The original FMANN counts y_j <= x_i + delta and reports that; we invert
## the count to match the definition of U used in wilcoxon_ci_2_toms516.
## O(m + n) once xs and ys are sorted.
fmann <- function(delta, xs, ys) {
  m <- length(xs)
  n <- length(ys)
  jle <- 0L
  ile <- 0L
  for (i in seq_len(m)) {
    threshold <- xs[i] + delta
    while (jle < n && ys[jle + 1L] <= threshold) {
      jle <- jle + 1L
    }
    ile <- ile + jle
    if (jle >= n) {
      ile <- ile + (m - i) * n
      break
    }
  }
  m * n - ile
}

## fmann_vec: same output as fmann, computed with findInterval.
##
## findInterval(v, sorted) walks a single pointer along `sorted` (which
## must be non-decreasing) and returns, for each element of v, the count
## of `sorted` values <= that element.  Because `sorted` is sorted and
## the walk is linear, complexity is O(m + n), the same as the two-pointer
## FMANN.  The gains over fmann are pure R-loop-overhead vs C-loop-
## overhead: for m = n in the hundreds the compiled findInterval loop is
## typically two to three orders of magnitude faster than the interpreted
## fmann loop.
fmann_vec <- function(delta, xs, ys) {
  length(xs) * length(ys) - sum(findInterval(xs + delta, ys))
}

## ill: Illinois-modified regula falsi.  Solves Fun(t) = target for a
## monotone Fun, given a bracket (x1, x2) with function values (f1, f2)
## straddling the target.  Faithful translation of the ILL subroutine of
## TOMS 516; see Dowell and Jarratt (1971) for the Illinois modification.
ill <- function(target, x1, f1, x2, f2, Fun, tol) {
  f1 <- f1 - target
  f2 <- f2 - target
  bisect_next <- FALSE
  repeat {
    if (abs(x2 - x1) < tol) break
    x3 <- if (bisect_next) (x1 + x2) / 2
          else             x2 - f2 * (x2 - x1) / (f2 - f1)
    bisect_next <- FALSE
    f3 <- Fun(x3) - target
    if (f3 * f2 <= 0) {
      x1 <- x2; f1 <- f2
      x2 <- x3; f2 <- f3
    } else {
      x2 <- x3; f2 <- f3
      f1 <- f1 / 2
      if (abs(f2) > abs(f1)) {
        f1 <- 2 * f1
        bisect_next <- TRUE
      }
    }
  }
  (x1 + x2) / 2
}
