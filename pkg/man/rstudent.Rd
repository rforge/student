\name{rstudent}
\alias{rstudent}
\title{Random Number Generator for the Multivariate Student t Distribution}
\description{
  Generate vectors of random variates from the multivariate Student \emph{t}
  distribution.
}
\usage{
rstudent(n, loc = rep(0, nrow(factor)), sigma, df = 3.5,
         factor = factorize(sigma))
}
\arguments{
  \item{n}{sample size \eqn{n} (positive integer).}
  \item{loc}{location vector of dimension \eqn{d}.}
  \item{sigma}{covariance matrix of dimension \eqn{(d,d)}.}
  \item{df}{degrees of freedom (positive real number or \code{Inf} in
    which case samples from the multivariate normal distribution with
    mean vector \code{loc} and covariance matrix \code{sigma} are
    drawn).}
  \item{factor}{factorization matrix of the covariance matrix
    \code{sigma}; a matrix \eqn{R} with \eqn{d} rows such that \eqn{R^T
    R} equals \code{sigma}. This matrix \eqn{R} is internally
    multiplied from the right to an \eqn{(n, d)}-matrix of independent
    standard normals to obtain a sample from a multivariate normal
    with zero mean vector and covariance matrix \code{sigma}.}
}
\value{
  \code{rstudent()} returns an \eqn{(n, d)}-matrix containing \eqn{n}
  samples of the \eqn{d}-dimensional Student \emph{t} distribution with
  \code{df} degrees of freedom, location vector \code{loc} and
  covariance matrix \code{sigma}.
}
\details{
  Internally used is \code{factor}, so \code{sigma} is not required
  to be provided if \code{factor} is given.

  The default factorization used to obtain \code{factor} is the Cholesky
  decomposition. To this end, \code{sigma} needs to have full rank.
  See \code{\link{factorize}()} for other options.

  The number of rows of \code{factor} are the dimension \eqn{d} of
  the sample. Typically (but not necessarily), \code{factor} is square.
%%    The default for the degrees of freedom were chosen for two
%%    reasons: first, because of financial return data, which often
%%    tends to have a degrees of freedom around 3.5; second, for
%%    pedagogical reasons, to show that the degrees of freedom
%%    parameter does not have to be an integer.
}
\author{Marius Hofert}
\references{
  McNeil, A. J., Frey, R., and Embrechts, P. (2015).
  \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
  Princeton University Press.
}
\seealso{
  \code{\link{factorize}()}
}
\examples{
## Generate a random correlation matrix in dimension 5
d <- 5
set.seed(271)
A <- matrix(runif(d * d), ncol = d)
P <- cov2cor(A \%*\% t(A))
## Draw random variates and compare
n <- 1000
set.seed(157)
X <- rstudent(n, sigma = P) # providing sigma
set.seed(157)
X. <- rstudent(n, factor = chol(P)) # directly the factor
stopifnot(all.equal(X, X.))
}
\keyword{distribution}