\name{rpdataconv}
\alias{rpdataconv}
\title{Convert a List of Numerics to a Matrix}

\description{
This function converts a list whose members are numeric vectors (possibly of different length) to a
matrix.
}

\usage{
rpdataconv(x)
}

\arguments{
    \item{x}{A list of numerics, a numeric vector or a matrix}
}

\details{
If the list entries do not have the same length, the shorter elements are \emph{not} recycled, but
filled with \code{NA} instead.
}

\value{
A matrix whose columns contain the members of \code{x}.
}

\note{
For consistency, \code{x} can also be a numeric vector or matrix. The return value is then the same
as the one of \code{\link[base]{as.matrix}}.
}

\seealso{
See \code{\link{rpdatasim}} for creating a similar data matrix with simulated data, and
\code{\link{rpdataboot}} for creating boostrap replications of such a data matrix.
}

\examples{
# Gemerate samples of different size from an
# exponential distribution with different parameters
x <- list(rexp(10, 0.2), rexp(7, 0.1), rexp(12, 0.5))

# Write x into a matrix that can be further processed
x.rp <- bootruin:::rpdataconv(x)
}

\keyword{array}
\keyword{manip}
\keyword{internal}
