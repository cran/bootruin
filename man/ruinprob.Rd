\encoding{utf8}
\name{ruinprob}
\alias{ruinprob}
\title{The Probability of Ruin in the Classical Risk Process}

\description{
This function calculates or estimates the probability of ruin in the classical (compund Poisson)
risk process using several different methods.
}

\usage{
ruinprob(x, param.list, compmethod = c("dg", "exp"),
    flmethod = c("nonp", "exp", "lnorm", "custom"),
    reserve, loading, fl = NA, interval = 0.5, ...)
}

\arguments{
    \item{x}{a numeric vector, matrix or array of individual claims.}
    \item{param.list}{a named list of parameters. It might contain any of the arguments except
    \code{x} and \code{\dots}}
    \item{compmethod}{a character string determining the algorithm for the computation.}
    \item{flmethod}{a character string indicating what cumulative probability distribution function
    is used for the increments of the running maximum of the aggregate loss process if
    \code{compmethod = "dg"}, see also Details and References.}
    \item{reserve}{a number indicating the initial surplus.}
    \item{loading}{a number determining the relative security loading.}
    \item{fl}{a function that is used as custom cumulative probability distribution to be used for
    the discretization if \code{flmethod = "custom"}.}
    \item{interval}{a number determining the approximation precision, viz. the mesh width of the
    discretization if \code{compmethod = "dg"}, see Details and References.}
    \item{\dots}{further arguments are passed to \code{fl}.}
}

\details{
The classical risk process, also called Cramér-Lundberg risk process, is a stochastic model for an
insurer's surplus over time and, for any \eqn{t\geq0}{t >= 0}, it is given by
\deqn{Y_{t} = r_{0} + ct - Z_{t},}{Y_t = r_0 + ct - Z_t,}
where \eqn{Z_{t}}{Z_t} is a compund Poisson process, \eqn{r_{0} \geq 0}{r_0 >= 0} is the initial
surplus and \eqn{c > 0} is the constant premium rate.

This function calculates, approximates or estimates (depending on what options are given) the
probability of ruin in the infinite time horizon, i.e. the probability that \eqn{Y_{t}}{Y_t} ever
falls below 0.

Currently there are two options for the \code{compmethod} argument. If \code{compmethod = "exp"},
the claims are assumed to be from an exponential distribution. In that case, the probability of ruin
is given by
\deqn{\frac{1}{1 + \beta} \exp\left\{-\frac{\beta}{1+\beta} \frac{r_{0}}{\mu}\right\},}{%
      exp\{-beta * r_0 / (mu (1 + beta))\} / (1 + beta),}
where \eqn{\mu}{mu} is the mean claim size (estimated from \code{x}) and \eqn{\beta}{beta} is the
relative security loading.

For \code{compmethod = "dg"}, the recursive algorithm due to Dufresne and Gerber (1989) is used.  In
this case, the parameter \code{flmethod} determines what cumulative distribution function is used
for the discretization. The possible choices are either a non-parametric estimator, parametric
estimators for exponential or log-normal claims, or a user-supplied function (in which the argument
\code{fl} must be specified).  See the reference for more details on how this algorithm works.
}

\value{
The estimated or calculated probability of ruin. The shape and dimension of the output depends on
the specifics of the claim data \code{x}. If \code{x} is a vector, the output is a single numeric
value. In general, the dimension of the output is one less than that of \code{x}. More precisely, if
\code{x} is an array, then the output value is an array of dimension \code{dim(x)[-1]}, see the note
below.
}

\references{
Dufresne, F. and Gerber, H.-U. (1989) \emph{Three Methods to Calculate the Probability of Ruin}.
ASTIN Bulletin, \bold{19}(1), pp. 71--90.
}

\note{
If \code{x} is an array rather than a vector, the function acts as if it was called through
\code{\link[base]{apply}} with \code{MARGIN = 2:length(dim(x))}

If an option is given both explicitly and as part of the \code{param.list} argument, then the value
given explicitly takes precedence. This way the parameter list, saved as a variable, can be reused,
but modifications of one or more parameter values are still possible.
}

\seealso{
\code{\link{ruinprob.test}}
}

\examples{
# Claims have an exponential distribution with mean 10
x <- rexp(10, 0.1)
print(x)

# The estimated probability of ruin
ruinprob(x, reserve = 100, loading = 0.2, interval = 0.25)

# The true probability of ruin of the risk process
ruinprob(
    10, reserve = 100, loading = 0.2,
    flmethod = "exp", compmethod = "exp"
)
}

\keyword{nonparametric}
