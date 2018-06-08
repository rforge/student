\name{student}
\alias{pstudent}
\title{Multivariate t distribution }
\description{
  Functionality for the multivariate t distribution (including evaluation for non-integer degrees of freedom).
}
\usage{
pstudent (a, b, R, nu, gam = 3.3, eps = 0.001, Nmax = 1e8, N = 10, n_init = 2^5, precond = TRUE)
}
\arguments{
  \item{a}{Vector of length \eqn{q}.}
  \item{b}{Vector of length \eqn{q}.}
  \item{R}{\eqn{(q,q)} covariance matrix. Must be positive definit.}
  \item{nu}{Degress of freedom. Any value \eqn{>0}.}
  \item{eps}{Error tolerance.}
  \item{gam}{Algorithm will run until \eqn{err < gam * eps}.}
  \item{Nmax}{Maximum number of function evaluations, can be used to control run time.}
  \item{N}{Number of repetitions to get error estimate in RQMC approach.}
  \item{n_init}{Size of the first point set being used to estimate the probability.}
  \item{precond}{A \code{\link{logical}} indicating if preconditioning (i.e. reordering the variables) is to be performed.}
  
}
\value{
  \code{pstudent()} returns a list of 4, including the estimated probability, the number of iterations, the total number of function evaluations and an error estimate. 
}
\details{
  Note that this procedure calls fast C code. Currently, the dimension \eqn{q} cannot exceed 16510. If \eqn{q=1}, the function calls the univarite \code{pt()}. 
}
\author{Marius Hofert, Erik Hintz and Christiane Lemieux}
\references{
  To come ...
}
%% \seealso{
%% }
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