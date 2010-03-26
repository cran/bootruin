\name{rppvalue}
\alias{rppvalue}
\title{P-values for the test of the probability of ruin}

\description{
This function provides p-values for the test of the probability of ruin using one of two different
methods.
}

\usage{
rppvalue(x, method = c("bootstrap", "normal"), x.boot)
}

\arguments{
    \item{x}{The observed values of the test statistic as numeric vector or matrix, see Details.}
    \item{method}{A character string determining the method used.}
    \item{x.boot}{The bootstrap replications that \code{x} should be compared with, see Details.}
}

\details{
This function is not intended to be used by itself, but rather in combination with
\code{\link{rpteststat}}. Hence, ideally, both \code{x} and \code{x.boot} stem from a call of the
latter.

If \code{method = "bootstrap"}, then bootstrap p-values are computed. The values of \code{x} are
compared to those of \code{x.boot}. The number of rows of \code{x.boot} has to match the number of
columns of \code{x} (or its length if it is a vector). For most applications, however, \code{x} will
be a single number and \code{x.boot} will be the bootstrap replications of \code{x}

For \code{method = "normal"} the p-values are computed using the asymptotic normal approximation of
the test statistic. \code{x} can be a vector, a matrix or an array of numerics.

The elements of \code{x} are interpreted as statistics of separate, independent tests, and adjusting
the p-values for multiple comparison may be necessary.
}

\value{
A numeric (vector, matrix or array) with the same dimension as \code{x}.
}

\references{
Baumgartner, B. and Gatto, R. (2010) \emph{A Bootstrap Test for the Probability of Ruin in the
Compound Poisson Risk Process}. To appear.
}

\note{
If \code{method = "normal"}, the argument \code{x.boot} is not used and a warning is issued if it is
still provided.
}

\seealso{
See \code{\link{rpteststat}} for the computation of the test statistics, \code{\link{ruinprob}} for
computating the probability of ruin, and \code{\link{rpjack}} for the computation of the standard
errors.

\code{\link[stats]{p.adjust}} from the package \pkg{stats} provides methods to adjust p-values for
multiple testing.
}

\keyword{htest}
\keyword{nonparametric}
\keyword{internal}
