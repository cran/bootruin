\name{rpjack}
\alias{rpjack}
\title{A Jackknife Estimator of the Standard Error of the Probability of Ruin}

\description{
This function computes the jackknife estimator of the standard error of the estimator of the
probability of ruin given observed claims.
}

\usage{
rpjack(x, ...)
}

\arguments{
    \item{x}{A vector or a matrix of numeric data.}
    \item{\dots}{Further options that are passed on to \code{\link{ruinprob}}.}
}

\details{
If \code{x} is a vector of observed claims, for each element of \code{x} the probability of ruin is
estimated with said element of \code{x} left out. The resulting vector of ruin probabilities has the
same length as \code{x}, and its standard error, properly rescaled, is used to approximate the
standard error of the estimator of the probability of ruin of \code{x}.

This procedure is applied column-wise if \code{x} is a matrix.
}

\value{
A numeric vector of length \code{ncol(as.matrix(x))}.
}

\references{
Efron, B. and Tibshirani, R. (1993) \emph{An Introduction to the Bootstrap}. Chapman and Hall.

Tukey, J. W. (1958) \emph{Bias and Confidence in Not Quite Large Samples}. The Annals of
Mathematical Statistics, \bold{29}(2), p. 614.

Quenouille, M. H. (1956) \emph{Note on Bias in Estimation}. Biometrika, \bold{43}, pp. 353--360.
}

\note{
The calculation of the jackknife standard error can be computationally intensive. In most cases the
computation time can be drastically reduced at the price of a slightly lower accuracy, viz. a higher
value for the \code{interval} argument of \code{ruinprob}.
}

\seealso{
\code{\link{ruinprob}} for valid options that can be used for \code{\dots}.
}

\keyword{nonparametric}
\keyword{distribution}
\keyword{internal}
