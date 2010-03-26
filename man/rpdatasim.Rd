\name{rpdatasim}
\alias{rpdatasim}
\title{Simulating data, shaped into a matrix}

\description{
This is a wrapper function for simulating data that have the correct structure for further
processing with \code{\link{rpdataboot}}.
}

\usage{
rpdatasim(n, replications, rdist, ...)
}

\arguments{
    \item{n}{The number of observations per sample}
    \item{replications}{the number of samples/replications}
    \item{rdist}{A function for simulating random variables}
    \item{\dots}{Further arguments for \code{rdist}}
}

\value{
A numeric matrix with \code{n} rows and \code{replications} columns.
}

\note{
Typical choices for \code{rdist} are \code{\link[stats]{rexp}} or \code{\link[stats]{rlnorm}}.
}

\seealso{
See \code{\link{rpdataconv}} for converting an existing dataset to a matrix, and
\code{\link{rpdataboot}} for creating boostrap replications of such a data matrix.
}

\examples{
# Generate 5 independent samples of size 10 from
# an exponential distribution with mean 10
x <- bootruin:::rpdatasim(n = 10, replications = 5, rexp, rate = 1/10)
}

\keyword{array}
\keyword{datagen}
\keyword{distribution}
\keyword{internal}
