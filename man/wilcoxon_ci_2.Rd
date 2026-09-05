\name{wilcoxon_ci_2}
\alias{wilcoxon_ci_2}
\alias{wilcoxon_ci_2_exact}
\alias{wilcoxon_ci_2_uniroot}
\alias{wilcoxon_ci_2_toms516}
\alias{wilcoxon_ci_2_toms516_port}
\alias{wilcoxon_ci_2_toms516_vec}
\title{Two-Sample Wilcoxon Confidence Interval and Hodges-Lehmann Estimate}
\description{
  Distribution-free confidence interval for the shift parameter
  \eqn{\Delta}{Delta} in the two-sample location model
  \eqn{X_i \sim F(x - \Delta)}{X_i ~ F(x - Delta)},
  \eqn{Y_j \sim F}{Y_j ~ F}, and the associated Hodges-Lehmann point
  estimate: the median of the \eqn{mn}{mn} pairwise differences
  \eqn{d_{ij} = x_i - y_j}{d[i,j] = x[i] - y[j]}.

  Positive \eqn{\Delta}{Delta} means \eqn{X}{X} is shifted to the right
  of \eqn{Y}{Y}.  This is the convention of
  \code{\link[stats]{wilcox.test}}: the null value \code{mu} refers to
  the shift of \eqn{X}{X} relative to \eqn{Y}{Y}, and
  \code{alternative = "greater"} tests whether \eqn{X}{X} is stochastically
  larger than \eqn{Y}{Y}.

  \code{wilcoxon_ci_2} is a dispatcher over three backends; call it and
  let it pick, or call a backend by name.
}
\usage{
wilcoxon_ci_2(x, y, conf.int = TRUE, conf.level = 0.95,
              method = c("exact", "uniroot", "toms516"), \dots)

wilcoxon_ci_2_exact  (x, y, conf.int = TRUE, conf.level = 0.95)
wilcoxon_ci_2_uniroot(x, y, conf.int = TRUE, conf.level = 0.95, \dots)
wilcoxon_ci_2_toms516(x, y, conf.int = TRUE, conf.level = 0.95,
                      variant = c("vec", "port"))
wilcoxon_ci_2_toms516_vec (x, y, conf.int = TRUE, conf.level = 0.95)
wilcoxon_ci_2_toms516_port(x, y, conf.int = TRUE, conf.level = 0.95)
}
\arguments{
  \item{x, y}{numeric vectors of observations from the two samples.}
  \item{conf.int}{logical; if \code{TRUE} (the default), a confidence
    interval is returned along with the point estimate.}
  \item{conf.level}{requested confidence level.  The achieved level (a
    component of the returned list) may differ because of the
    discreteness of the exact null distribution.}
  \item{method}{character; the backend to use.  If not supplied,
    \code{wilcoxon_ci_2} chooses \code{"exact"} when
    \code{max(length(x), length(y)) < 50} and \code{c(x, y)} has no
    ties, else \code{"uniroot"}.  \code{"toms516"} is never chosen
    automatically; request it by name.}
  \item{\dots}{additional arguments forwarded to
    \code{\link[stats]{uniroot}} by the \code{"uniroot"} backend.  The
    most useful is \code{tol} (default
    \code{.Machine$double.eps^0.25} \eqn{\approx 1.2 \times 10^{-4}}).
    Ignored by the \code{"exact"} and \code{"toms516"} backends.}
  \item{variant}{character; only for \code{wilcoxon_ci_2_toms516}.
    Selects the \eqn{U}{U} evaluator: \code{"vec"} (default,
    \code{findInterval}-based) or \code{"port"} (interpreted R
    two-pointer walk).  Results are identical up to numerical noise;
    the choice is a speed/pedagogy trade-off.}
}
\details{
  \strong{Exact backend.}  Materialises all \eqn{mn}{mn} pairwise
  differences \eqn{d_{ij} = x_i - y_j}{d[i,j] = x[i] - y[j]}, sorts
  them, and picks order statistics indexed by
  \code{\link[stats]{qwilcox}}.  This is Bauer's (1972) direct
  construction and is the reference implementation.  Memory and time
  are \eqn{O(mn)}{O(mn)}.

  \strong{Uniroot backend.}  Solves for the confidence bounds by
  root-finding on the rank-sum statistic
  \eqn{W(\Delta) = \sum_i R(x_i - \Delta)}{W(Delta) = sum_i R(x_i - Delta)}
  computed on the combined sample
  \eqn{\{x_1 - \Delta, \dots, x_m - \Delta, y_1, \dots, y_n\}}, using
  \code{\link[stats]{uniroot}}.  Root-finding targets are half-integers
  so the step function is never at the target value.  Evaluator is
  \eqn{O((m+n)\log(m+n))}{O((m+n) log(m+n))}; no
  \eqn{mn}{mn}-element matrix is materialised.

  \strong{TOMS 516 backend.}  R port of the Illinois-modified regula
  falsi algorithm of Ryan and McKean (1977).  The Mann-Whitney
  \eqn{U(\Delta) = \#\{(i,j) : x_i - y_j > \Delta\}}{U(Delta) = #{(i,j) : x[i] - y[j] > Delta}}
  is evaluated in \eqn{O(m+n)}{O(m+n)} by a two-pointer walk over the
  sorted arrays, and the Illinois-modified regula falsi root-finder
  (\code{ill}) inverts \eqn{U}{U} at half-integer targets.  Two variants
  are provided: \code{"vec"} (the default) uses
  \code{\link[base]{findInterval}} for the \eqn{U}{U} evaluator, giving
  a compiled inner loop and a large constant-factor speedup over
  \code{"port"}, which is a faithful line-by-line port of the original
  Fortran with the two-pointer walk written in R.  Deviations from the
  original Fortran are documented in the source of
  \code{wilcoxon_ci_2_toms516_port}.

  \strong{Why not use \code{wilcox.test}'s confidence interval?}  See
  the \sQuote{Details} section of \code{\link{wilcoxon_test}}.
}
\value{
  A list with components
  \item{estimate}{the Hodges-Lehmann point estimate
    \eqn{\mathrm{median}\{x_i - y_j\}}{median{x_i - y_j}}.}
  \item{conf.int}{present only when \code{conf.int = TRUE}: a length-two
    numeric vector giving the confidence interval.}
  \item{conf.level}{present only when \code{conf.int = TRUE}: the
    achieved confidence level.}
}
\references{
  Bauer, D.\ F.\ (1972).  Constructing confidence sets using rank
  statistics.  \emph{Journal of the American Statistical Association}
  \bold{67}(339), 687--690.

  Hodges, J.\ L., and Lehmann, E.\ L.\ (1963).  Estimates of location
  based on rank tests.  \emph{Annals of Mathematical Statistics}
  \bold{34}(2), 598--611.

  McKean, J.\ W., and Ryan, T.\ A.\ (1977).  An algorithm for obtaining
  confidence intervals and point estimates based on ranks in the
  two-sample location problem.  \emph{ACM Transactions on Mathematical
  Software} \bold{3}(2), 183--185.

  Dowell, M., and Jarratt, P.\ (1971).  A modified regula falsi method
  for computing the root of an equation.  \emph{BIT} \bold{11}(2),
  168--174.
}
\seealso{
  \code{\link{wilcoxon_test}}, \code{\link{wilcoxon_ci_1}},
  \code{\link[stats]{qwilcox}}, \code{\link[stats]{uniroot}}
}
\examples{
set.seed(1)
x <- rnorm(20)
y <- rnorm(25, mean = 0.5)

wilcoxon_ci_2(x, y)
wilcoxon_ci_2(x, y, method = "uniroot")
wilcoxon_ci_2(x, y, method = "toms516")
wilcoxon_ci_2_exact(x, y)
wilcoxon_ci_2(x, y, conf.int = FALSE)

## Tighter uniroot tolerance
wilcoxon_ci_2(x, y, method = "uniroot", tol = 1e-8)
}
\keyword{htest}
\keyword{nonparametric}
