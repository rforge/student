\name{rstudent}
\alias{rstudent}
\title{Random Number Generator for the Multivariate Student t Distribution}
\description{
  Generate vectors of random variates from the multivariate Student \eqn{t}
  distribution.
}
\usage{
rstudent(n, loc = rep(0, ncol(sigma)), sigma = diag(2), df = 3.5,
         factor = chol(sigma, pivot = TRUE))
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
    \code{sigma}; a matrix \eqn{R} with \eqn{d} rows which is internally
    multiplied from the right to an \eqn{(n, d)}-matrix of independent
    standard normals to obtain a sample to obtain a sample from
    a multivariate normal with zero mean vector and covariance matrix
    \code{sigma}.}
}
\value{
  \code{rstudent()} returns an \eqn{(n, d)}-matrix containing \eqn{n}
  samples of the \eqn{d}-dimensional Student $t$ distribution with
  \code{df} degrees of freedom, location vector \code{loc} and
  covariance matrix \code{sigma}.
}
\details{
  TODO: fix 't' problem above; mention why df 3.5; make nice example below
}
\author{Marius Hofert}
%% \references{
%%   TODO: QRM book (everywhere; make other help files, too)
%% }
\examples{
## Generate a random correlation matrix in dimension 5
d <- 5
set.seed(271)
A <- matrix(runif(d * d), ncol = d)
P <- cov2cor(A \%*\% t(A))
## Draw random variates
n <- 1000
X <- rstudent(n, sigma = P)
## Compare with providing A directly (different sample than default) and with providing R directly (same sample)
}
\keyword{distribution}