\name{rpteststat}
\alias{rpteststat}
\title{Computation of the Test Statistic for the Probability of Ruin}

\description{
This function computes the test statistic for the test of the probability of ruin
}

\usage{
rpteststat(x, x.null, se)
}

\arguments{
    \item{x}{The observed/estimated probability of ruin}
    \item{x.null}{The null probability of ruin}
    \item{se}{The standard error of \code{x}}
}

\details{
This function studentizes (i.e., centers and rescales) the observed probabilities of ruin in a way
that they can be further processed by \code{\link{rppvalue}}. This latter computes p-values for the
test with the null hypothesis that the probability of ruin is equal to \code{x.null} versus the
one-sided alternative that probability of ruin is smaller than \code{x.null}

There are two possible scenarios how to use this function:

The first one is to compute the test statistic for vector of estimated probabilities of ruin (where
each element corresponds to an independent sample of observed claims) using a pre-determined value
of \code{x.null}. In this case the output is a matrix with \code{length(x)} columns and
\code{length(x.null)} rows.

The second use is to compute the bootstrap replications of the test statistic. In that case the
estimated probabilities of ruin are used for \code{x.null}, and the numbers of rows of \code{x} has
to match the length of \code{x.null}, i.e. there is one row of bootstrap replications of the
probability of ruin for each sample of claim sizes.
}

\value{
A numeric matrix that contains the test statistic. Its dimension depends on the input, see Details.
}

\note{
The dimensions (or lengths) of \code{x} and \code{se} have to be the same.
}

\seealso{
See \code{\link{ruinprob}} on how to obtain the estimators of the probabilities of ruin,
\code{\link{rpjack}} on how to get approximate standard errors of the former, and
\code{\link{rppvalue}} on how to compute the p-values of the test described here.
}

\keyword{htest}
\keyword{nonparametric}
\keyword{internal}
