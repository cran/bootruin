\name{rpdataboot}
\alias{rpdataboot}
\title{Creating Bootstrap Replications from an Matrix of Observations}

\description{
This function provides a simple way to create bootstrap replications of a dataset. The replication
is either non-parametrical or parametrical (for exponential or logarithmic normal data).
}

\usage{
rpdataboot(x, b, method = c("nonp", "exp", "lnorm"))
}

\arguments{
    \item{x}{A matrix containing numeric values}
    \item{b}{The number of bootstrap replications}
    \item{method}{The replication method applied to the data}
}

\details{
The input matrix \code{x} is supposed to contain (independent) observations in each column. The
bootstrap replication take this into account is done column-wise.

Depending on how the boostrap replications are further processed, the boostrap resampling should be
done either non-parametrically (\code{method = "nonp"}) or parametrically.

In the non-parametrical case, the bootstrap replications are samples drawn from the empirical
distribution of the original observation, this is equivalent to drawing with replacement.

For the parametrical bootstrap replications there are currently two options: With \code{method =
"exp"} each bootstrap replication is a vector simulated from an exponential distribution function
whose parameter is estimated by the original observation. For \code{method = "lnorm"} the resampling
is done by simulating from a logarithmic normal distribution whose log-mean and log-variance are
estimated from the original observation.
}

\value{
An array of dimension \code{c(dim(x), b)} containing column-wise bootstrap replications of \code{x}
}

\note{
\code{NA}'s are propagated consistently. More precisely, only the non-\code{NA} values undergo the
resampling and thus, missing values remain unchanged in the bootstrap replications.
}

\seealso{
\code{\link{rpdataconv}} for creating a suitable data matrix from a list of observation vectors, and
\code{\link{rpdatasim}} for creating such a matrix by simulation.
}

\examples{
# Generate a data matrix of 5 samples with 10 observations each.
x <- matrix(rexp(50), nrow = 10, ncol = 5)

# Create (parametric) bootstrap replications
x.boot.par <- bootruin:::rpdataboot(x, b = 50, method = "exp")

# Create (non-parametric) bootstrap replications
x.boot.nonp <- bootruin:::rpdataboot(x, b = 50, method = "nonp")
}

\keyword{array}
\keyword{datagen}
\keyword{distribution}
\keyword{internal}
