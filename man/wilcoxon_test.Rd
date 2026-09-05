\name{wilcoxon_test}
\alias{wilcoxon_test}
\title{Wilcoxon Tests with Hodges-Lehmann Estimation}
\description{
  Performs one- and two-sample Wilcoxon tests with Hodges-Lehmann point
  estimation and distribution-free confidence intervals for the location
  parameter (one-sample or paired) or the shift parameter (two-sample).

  The test statistic and p-value are delegated to
  \code{\link[stats]{wilcox.test}}.  The point estimate and confidence
  interval are computed by \code{\link{wilcoxon_ci_1}} (one-sample or
  paired) or \code{\link{wilcoxon_ci_2}} (two-sample), which use
  correctly-indexed exact constructions rather than the off-by-one
  scheme in \code{wilcox.test} (see \sQuote{Details}).
}
\usage{
wilcoxon_test(x, y = NULL,
              alternative = c("two.sided", "less", "greater"),
              mu = 0,
              paired = FALSE,
              test = TRUE,
              conf.int = TRUE,
              conf.level = 0.95,
              method = c("exact", "uniroot", "toms516"),
              \dots)
}
\arguments{
  \item{x}{numeric vector of data values.  Non-finite values are
    dropped.}
  \item{y}{optional numeric vector of data values.  Non-finite values
    are dropped.}
  \item{alternative}{character; one of \code{"two.sided"} (default),
    \code{"less"}, or \code{"greater"}.  Passed to
    \code{\link[stats]{wilcox.test}}.}
  \item{mu}{null value of the location parameter (one-sample or paired)
    or shift parameter (two-sample).  Passed to
    \code{\link[stats]{wilcox.test}}.  Does not enter the point estimate
    or confidence interval, which are for the parameter itself.}
  \item{paired}{logical; if \code{TRUE}, a paired test is performed and
    the problem reduces to a one-sample problem on \code{x - y}.}
  \item{test}{logical; if \code{TRUE} (the default), the test statistic
    and p-value are computed via \code{\link[stats]{wilcox.test}} and
    included in the returned object.  If \code{FALSE}, the test is
    skipped and the returned object carries only the point estimate and
    (if requested) the confidence interval.  Any arguments in
    \code{\dots} are ignored when \code{test = FALSE}.}
  \item{conf.int}{logical; if \code{TRUE} (the default), a confidence
    interval for the parameter is included in the returned object.}
  \item{conf.level}{confidence level of the interval.  The interval's
    actual level, reported via \code{attr(conf.int, "conf.level")}, may
    differ from this because of the discreteness of the exact
    distribution.}
  \item{method}{character; algorithm used to compute the point estimate
    and confidence interval.  See \sQuote{Details}.  If missing, the
    default is chosen automatically.}
  \item{\dots}{further arguments passed to
    \code{\link[stats]{wilcox.test}} (e.g.\ \code{exact}, \code{correct},
    \code{tol.root}, \code{digits.rank}).  Note that the \code{exact}
    argument here controls the p-value computation only; the
    \code{method} argument controls the CI.}
}
\details{
  \strong{Sign convention.}  For the two-sample problem the model is
  \eqn{X_i \sim F(x - \Delta)}{X_i ~ F(x - Delta)}, \eqn{Y_j \sim F}{Y_j ~ F};
  positive \eqn{\Delta}{Delta} means \eqn{X}{X} is shifted to the right
  of \eqn{Y}{Y}, so the point estimate is \eqn{\mathrm{median}\{x_i -
  y_j\}}{median{x_i - y_j}} and \code{alternative = "greater"} tests
  whether \eqn{X}{X} is stochastically larger than \eqn{Y}{Y}.  This
  matches \code{\link[stats]{wilcox.test}}.

  \strong{Why not use \code{wilcox.test}'s confidence interval?}
  \code{stats::wilcox.test} returns a confidence interval with off-by-one
  indexing into the sorted Walsh averages (one-sample) or pairwise
  differences (two-sample).  For example, with \eqn{n} observations and
  critical value \eqn{c_2} from \code{qsignrank}, the correct interval
  is \eqn{(A_{(c_2 + 1)}, A_{(N - c_2)})} where \eqn{N = n(n+1)/2}, but
  \code{wilcox.test} returns \eqn{(A_{(c_2)}, A_{(N - c_2 - 1)})}.  The
  difference is usually small in absolute terms but produces incorrect
  achieved coverage.  See the references for the direct construction.

  \strong{Method selection.}  The \code{method} argument is forwarded to
  \code{\link{wilcoxon_ci_1}} or \code{\link{wilcoxon_ci_2}}; when not
  supplied, those dispatchers auto-select between \code{"exact"} (direct
  construction) and \code{"uniroot"} (root-finding on the signed-rank or
  Mann-Whitney statistic) by the same rule
  \code{\link[stats]{wilcox.test}} uses to switch between exact and
  asymptotic p-values: exact when all samples have fewer than 50
  observations and no ties, uniroot otherwise.  \code{"toms516"} (an R
  port of the Illinois-modified regula falsi algorithm of Ryan and
  McKean, TOMS Alg.\ 516) is available for two-sample problems only and
  is never chosen automatically; request it by name.

  \strong{Tuning the uniroot backend.}  \code{wilcoxon_test} reserves
  \code{\dots} for pass-through to \code{\link[stats]{wilcox.test}}, so
  \code{uniroot}'s tolerance and other tuning knobs are not reachable
  from here.  To tune them, call \code{\link{wilcoxon_ci_1}} or
  \code{\link{wilcoxon_ci_2}} directly with \code{method = "uniroot"} and
  any \code{uniroot} arguments (e.g.\ \code{tol = 1e-8}); those
  dispatchers forward \code{\dots} to \code{uniroot}.  No computation is
  duplicated: if the test statistic and p-value are also needed, get
  them from \code{wilcoxon_test(\dots, conf.int = FALSE)} or
  \code{\link[stats]{wilcox.test}}.

  \strong{Paired data.}  If \code{paired = TRUE}, the problem is
  reduced to a one-sample problem on the differences \code{x - y} and
  \code{wilcoxon_ci_1} is used.  \code{method = "toms516"} is not
  available in this case.
}
\value{
  An object of class \code{c("wilcoxon_test", "htest")} with components
  \item{statistic}{the test statistic, as returned by
    \code{\link[stats]{wilcox.test}} (\code{V} for the one-sample or
    paired case, \code{W} for the two-sample case).}
  \item{p.value}{the p-value from \code{\link[stats]{wilcox.test}}.}
  \item{null.value}{the null value of the parameter, from
    \code{\link[stats]{wilcox.test}}.}
  \item{alternative}{the alternative hypothesis.}
  \item{method}{a character string describing the test, from
    \code{\link[stats]{wilcox.test}}.}
  \item{data.name}{a character string identifying the data.}
  \item{estimate}{the Hodges-Lehmann point estimate: the pseudo-median
    of \code{x} (one-sample), of \code{x - y} (paired), or the median of
    the mn pairwise differences \eqn{x_i - y_j} (two-sample).}
  \item{conf.int}{present only when \code{conf.int = TRUE}: a length-two
    numeric vector giving the confidence interval, with the achieved
    confidence level in an attribute named \code{"conf.level"}.}
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
}
\seealso{
  \code{\link{wilcoxon_ci_1}}, \code{\link{wilcoxon_ci_2}},
  \code{\link[stats]{wilcox.test}}
}
\examples{
## One-sample
set.seed(1)
x <- rnorm(20)
wilcoxon_test(x)

## Two-sample
y <- rnorm(15, mean = 0.5)
wilcoxon_test(x, y)

## Paired
y2 <- x + rnorm(20, mean = 0.3)
wilcoxon_test(x, y2, paired = TRUE)

## Requesting a specific CI backend
wilcoxon_test(x, y, method = "toms516")

## Test without CI
wilcoxon_test(x, conf.int = FALSE)

## Estimation only, no test
wilcoxon_test(x, y, test = FALSE)
}
\keyword{htest}
\keyword{nonparametric}
