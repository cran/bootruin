\name{pvaldens}
\alias{pvaldens}
\title{Density Estimation of Data in the Unit Interval}

\description{
This function computes density estimators for densities with the unit interval as support.  One
example of data with such a density are p-values.  Currently, two methods are implemented that
differ in the kernel function used for estimation.
}

\usage{
pvaldens(x, bw, rho, method = c("jh", "chen"))
}

\arguments{
    \item{x}{a numeric vector of data points between 0 and 1.}
    \item{bw}{a number indicating the bandwidth used for the density estimation.}
    \item{rho}{a number determining the correlation coefficient, only used if \code{method = "jh"}}
    \item{method}{a character string determining the kernel function that is used, see Details.}
}

\details{
Depending on which \code{method} is selected, a different kernel function is used for the
estimation.  Since the support of the estimated function is bounded, those kernel functions are
location-dependent.

If \code{method = "jh"}, a Gaussian copula-based kernel function according to Jones and Henderson
(2007) is used. In this case the bandwidth can either be specified directly or as correlation
coefficient: if \eqn{\rho>0}{rho > 0} denotes the correlation coefficient and \eqn{h>0}{h > 0} the
bandwidth, then \eqn{h^2=1-\rho}{h^2 = 1 - rho}. Note that \code{rho} and \code{bw} are mutually
exclusive.

For \code{method = "chen"}, the kernel function is based on a beta density, according to Chen
(1999).

See the cited articles for more details.
}

\value{
A function with a single vector-valued argument that returns the estimated density at any given
point(s).
}

\references{
Jones, M. C. and Henderson, D. A. (2007) \emph{Kernel-Type Density Estimation on the Unit
Interval}. Biometrika, \bold{94}(4), pp. 977--984.

Chen, S. X. (1999) \emph{A Beta Kernel Estimation for Density Functions}. Computational Statistics
and Data Analysis, \bold{31}(2), pp. 131--145.
}

\seealso{
    \code{\link[stats]{density}}
}

\examples{
require(graphics)

x <- rbeta(100, 2, 5)
fhat <- pvaldens(x, rho = 0.9, method = "jh")

hist(x, freq = FALSE, xlim = c(0, 1))
curve(fhat(x), from = 0, to = 1, add = TRUE, col = 2)
box()
}

\keyword{dgraph}
\keyword{distribution}
\keyword{univar}
