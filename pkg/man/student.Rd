\name{student}
\alias{pstudent}
\title{Multivariate Student t Distribution}
\description{
  Evaluation of the multivariate Student t distribution function
  (including non-integer degrees of freedom).
}
\usage{
pstudent(a, b, R, nu, gam = 3.3, eps = 0.001, Nmax = 1e8, N = 10,
         n_init = 2^5, precond = TRUE)
}
\arguments{
  \item{a}{vector of length \eqn{q}.}
  \item{b}{vector of length \eqn{q}.}
  \item{R}{positive definite \eqn{(q,q)}-covariance matrix.}
  \item{nu}{degress of freedom (any positive value).}
  \item{eps}{error tolerance.}
  \item{gam}{determines the stopping criterion of the algorithm; it will
    run until \eqn{err < gam * eps}.}
  \item{Nmax}{maximum number of function evaluations, can be used to
    control run time.}
  \item{N}{Number of repetitions to get an error estimate in the
    randomized quasi-Monte Carlo approach.}
  \item{n_init}{size of the first point set being used to estimate
    the probability.}
  \item{precond}{\code{\link{logical}} indicating if preconditioning
    is applied, that is, reordering the variables.}
}
\value{
  \code{pstudent()} returns a list of length four, containing the
  the estimated probabilities, the number of iterations, the total
  number of function evaluations and an error estimate.
}
\details{
  Note that this procedure calls underlying C code. Currently, the
  dimension \eqn{q} cannot exceed 16510. If \eqn{q=1}, the function
  calls the univarite \code{pt()}.
}
\author{Marius Hofert, Erik Hintz and Christiane Lemieux}
\references{
  TODO
}
\examples{
## Generate a random correlation matrix in dimension 5
dim <- 5
R <- matrix(runif(dim*dim), ncol = dim)
R <- cov2cor(R\%*\%t(R))

## Generate random limits
a <- runif(dim) * sqrt(dim) * (-3)
b <- runif(dim) * sqrt(dim) * 3

## Specify degrees of freedom
nu <- 1.1

## Evaluate the cdf
pstudent(a = a, b = b, R = R, nu = nu)
}
\keyword{distribution}