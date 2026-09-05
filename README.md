# wilcoxon

Wilcoxon signed-rank and rank-sum inference for the one- and two-sample
location problems, with Hodges-Lehmann point estimates and distribution-free
confidence intervals.

The main driver, `wilcoxon_test`, returns an `htest`-compatible object.  The
test statistic and p-value come from `stats::wilcox.test`; the point estimate
and confidence interval use validated implementations that avoid the
off-by-one indexing in `wilcox.test`'s own CI.

## Installation

From GitHub:

```r
remotes::install_github("kloke/wilcoxon")
```

## Example

```r
library(wilcoxon)

set.seed(1)
x <- rnorm(20)
y <- rnorm(25, mean = 0.5)

wilcoxon_test(x, y)
```

Direct access to the CI functions is also available for tuning or programmatic
use:

```r
wilcoxon_ci_1(x)                                  # one-sample
wilcoxon_ci_2(x, y, method = "uniroot", tol = 1e-8)  # tuned root-finder
```

## Backends

Three algorithms are available for the confidence interval:

- `"exact"` — Bauer's direct construction; default for small samples with no
  ties.
- `"uniroot"` — root-finding via `stats::uniroot`; default for larger samples
  or when ties are present.  Extra arguments are forwarded to `uniroot`.
- `"toms516"` — R port of the Illinois-modified regula falsi algorithm of
  Ryan and McKean (1977, TOMS Algorithm 516).  Two-sample only; never chosen
  automatically.

The two-sample problem follows the `stats::wilcox.test` convention: the model
is X = Y + Δ, so positive Δ means X shifted to the right of Y.

## Citation

Kloke, J.D. and McKean, J.W. (2024), *Nonparametric Statistical Methods Using
R, 2nd Edition*.  Boca Raton, FL: Chapman-Hall.

## Fine Print

As with R itself, wilcoxon is free software and comes with ABSOLUTELY NO
WARRANTY.
