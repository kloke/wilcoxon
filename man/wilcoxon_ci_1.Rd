\name{wilcoxon_ci_1}
\alias{wilcoxon_ci_1}
\alias{wilcoxon_ci_1_exact}
\alias{wilcoxon_ci_1_uniroot}
\title{One-Sample Wilcoxon Confidence Interval and Hodges-Lehmann Estimate}
\description{
  Distribution-free confidence interval for the location parameter of a
  symmetric continuous distribution, and the associated Hodges-Lehmann
  point estimate: the median of the Walsh averages
  \eqn{A_{ij} = (x_i + x_j)/2}{A[i,j] = (x[i] + x[j])/2},
  \eqn{i \le j}{i <= j}.

  \code{wilcoxon_ci_1} is a dispatcher over two backends; call it and let
  it pick, or call a backend by name.
}
\usage{
wilcoxon_ci_1(x, conf.int = TRUE, conf.level = 0.95,
              method = c("exact", "uniroot"), \dots)

wilcoxon_ci_1_exact  (x, conf.int = TRUE, conf.level = 0.95)
wilcoxon_ci_1_uniroot(x, conf.int = TRUE, conf.level = 0.95, \dots)
}
\arguments{
  \item{x}{numeric vector of observations.}
  \item{conf.int}{logical; if \code{TRUE} (the default), a confidence
    interval is returned along with the point estimate.}
  \item{conf.level}{requested confidence level.  The achieved level (a
    component of the returned list) may differ because of the
    discreteness of the exact null distribution.}
  \item{method}{character; the backend to use.  If not supplied,
    \code{wilcoxon_ci_1} chooses \code{"exact"} when \code{length(x) < 50}
    and \code{x} has no ties, else \code{"uniroot"}.  Same rule
    \code{\link[stats]{wilcox.test}} uses to decide between exact and
    asymptotic p-values.}
  \item{\dots}{additional arguments forwarded to
    \code{\link[stats]{uniroot}} by the \code{"uniroot"} backend.  The
    most useful is \code{tol} (default
    \code{.Machine$double.eps^0.25} \eqn{\approx 1.2 \times 10^{-4}}).
    Ignored by the \code{"exact"} backend.}
}
\details{
  \strong{Exact backend.}  Materialises the \eqn{n(n+1)/2}{n(n+1)/2}
  Walsh averages, sorts them, and picks order statistics indexed by
  \code{\link[stats]{qsignrank}}.  This is Bauer's (1972) direct
  construction and is the reference implementation.  Memory and time are
  \eqn{O(n^2)}{O(n^2)}.

  \strong{Uniroot backend.}  Solves for the confidence bounds by
  root-finding on the signed-rank gradient
  \eqn{W(\theta) = \sum_i \mathrm{sign}(x_i - \theta) \cdot R(|x_i - \theta|)}{W(theta) = sum_i sign(x_i - theta) * R(|x_i - theta|)}
  with \code{\link[stats]{uniroot}}.  Root-finding targets are chosen so
  the step function \eqn{W} is never at the target value, so
  \code{uniroot} converges to the correct jump point rather than into
  the interior of a plateau.  For \eqn{n \le 50}{n <= 50} the critical
  value comes from \code{\link[stats]{qsignrank}}; beyond that a normal
  approximation with continuity correction is used because
  \code{qsignrank}'s \eqn{O(N^2)}{O(N^2)} recursion (\eqn{N = n(n+1)/2}{N = n(n+1)/2})
  otherwise dominates wall time.

  \strong{Why not use \code{wilcox.test}'s confidence interval?}  See
  the \sQuote{Details} section of \code{\link{wilcoxon_test}}.
}
\value{
  A list with components
  \item{estimate}{the Hodges-Lehmann point estimate.}
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
}
\seealso{
  \code{\link{wilcoxon_test}}, \code{\link{wilcoxon_ci_2}},
  \code{\link[stats]{qsignrank}}, \code{\link[stats]{uniroot}}
}
\examples{
set.seed(1)
x <- rnorm(20)

wilcoxon_ci_1(x)
wilcoxon_ci_1(x, method = "uniroot")
wilcoxon_ci_1_exact(x)
wilcoxon_ci_1(x, conf.int = FALSE)

## Tighter uniroot tolerance
wilcoxon_ci_1(x, method = "uniroot", tol = 1e-8)
}
\keyword{htest}
\keyword{nonparametric}
