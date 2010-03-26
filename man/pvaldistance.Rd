\encoding{utf8}
\name{pvaldistance}
\alias{pvaldistance}
\title{Distance measures of empirical probability functions}

\description{
This function provides a framework to evaluate various measures of distance between an empirical
distribution (induced by the dataset provided) and a theoretical probability distribution.
}

\usage{
pvaldistance(x, method = c("ks", "cvm"), dist.to = c("uniform"))
}

\arguments{
    \item{x}{a numeric vector containing a data sample.}
    \item{method}{a character string indicating which measure of distance is computed.}
    \item{dist.to}{a character string determining the (theoretical) probability distribution that is
    used as a reference.}
}

\details{
\code{method = "ks"} gives the Kolmogorov-Smirnov distance.

\code{method = "cvm"} yields the Cramér-von-Mises criterion (scaled with the sample size).
}

\value{
A positive real number giving the distance measure.
}

\note{
At the moment, \code{dist.to = "uniform"} (the uniform distribution on the unit interval) is the
only valid option for the theoretical distribution, and hence the members of \code{x} have to lie in
the unit interval.
}

\seealso{
See \code{\link[stats]{ks.test}} for the Kolmogorov-Smirnov test.
}

\examples{
# A sample from the standard uniform distribution
x <- runif(100, 0, 1)

# Distance to uniformity should be small
pvaldistance(x, "ks")
pvaldistance(x, "cvm")

# A sample from the Beta(2, 7) distribution
y <- rbeta(100, 2, 7)

# Distance to uniformity should be much larger here
pvaldistance(y, "ks")
pvaldistance(y, "cvm")
}

\keyword{distribution}
\keyword{univar}
