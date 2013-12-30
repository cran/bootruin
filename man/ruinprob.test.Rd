\name{ruinprob.test}
\alias{ruinprob.test}
\title{A Bootstrap Test for the Probability of Ruin in the Classical Risk Process}

\description{
This function provides a testing framework for the probability of ruin in the classical, compound
Poisson risk process. The test can be performed using the bootstrap method or using normal
approximation.
}

\usage{
ruinprob.test(x, prob.null, type = c("bootstrap", "normal"),
    nboot, bootmethod = c("nonp", "exp", "lnorm"), \dots)
}

\arguments{
    \item{x}{a numeric vector of data values (claims)}
    \item{prob.null}{a number indicating the hypothesized true probability of ruin.}
    \item{type}{a character string determining the type of test that is performed.}
    \item{nboot}{a number indicating the number of bootstrap replications.}
    \item{bootmethod}{a character string determining how the bootstrap replications are created.}
    \item{\dots}{further arguments to be passed to \code{\link{ruinprob}}.}
}

\details{
The null hypothesis is that the probability of ruin is equal to \code{prob.null} versus the
one-sided alternative that probability of ruin is smaller than \code{prob.null}.

If \code{type = "bootstrap"}, a bootstrap test is performed. The arguments \code{nboot} and
\code{bootmethod} have to be specified. \code{bootmethod} determines the kind of bootstrap:
\code{"nonp"} creates the usual nonparametric bootstrap replications, while \code{"exp"} and
\code{"lnorm"} create parametric bootstrap replications, the former assuming exponentially
distributed claims, the latter log-normally distributed ones.

\code{type = "normal"} makes use of an asymptotic normal approximation. The computations are a lot
faster, but from a theoretical point of view the bootstrap method is more accurate, see References.

For details about the necessary and valid arguments that might have to be supplied for \code{\dots},
see \code{\link{ruinprob}}.
}

\value{
A list with class \code{"htest"} containing the following components:
    \item{statistic}{the value of the studentized probability of ruin, i.e. the test statistic.}
    \item{parameter}{additional parameters.}
    \item{p.value}{the p-value for the test.}
    \item{estimate}{the estimated probability of ruin.}
    \item{null.value}{the specified hypothesized value of the probability of ruin.}
    \item{alternative}{a character string describing the alternative hypothesis.}
    \item{method}{a character string indicating what type of test was performed.}
    \item{data.name}{a character string giving the name of the data.}
}

\references{
Baumgartner, B. and Gatto, R. (2010) \emph{A Bootstrap Test for the Probability of Ruin in the
Compound Poisson Risk Process}. ASTIN Bulletin, \bold{40}(1), pp. 241--255.
}

\note{
Using the bootstrap method is computationally intensive. Values for \code{nboot} should not be too
large, usually numbers between 50 and 200 are reasonable choices.
}

\seealso{
\code{\link{ruinprob}}
}

\examples{
# Generating a sample of 50 exponentially distributed claims with mean 10
x <- rexp(50, 0.1)

\dontrun{
# Given this sample, test whether the probability of ruin is smaller than
# 0.1 using a bootstrap test with 100 bootstrap replications.
ruinprob.test(
    x = x, prob.null = 0.10, type = "bootstrap",
    loading = 0.2, reserve = 100, interval = 1,
    bootmethod = "nonp", nboot = 100
)
}

# The same test using normal approximation. This is a lot faster.
ruinprob.test(
    x = x, prob.null = 0.15, type = "normal",
    loading = 0.2, reserve = 100, interval = 1
)
}

\keyword{htest}
