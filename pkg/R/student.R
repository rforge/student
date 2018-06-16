### dStudent() #################################################################

##' @title Density of a Multivariate Student t Distribution
##' @param x (n, d)-matrix of evaluation points
##' @param df degrees of freedom (positive real or Inf in which case the density
##'        of a N(loc, scale) is evaluated)
##' @param loc location vector of dimension d
##' @param scale covariance matrix of dimension (d, d)
##' @param factor factorization matrix of the covariance matrix scale;
##'        caution: this has to be an *upper triangular* matrix R
##'        such that R^T R = scale here (otherwise det(scale) not computed correctly)
##' @return n-vector with t_nu(loc, scale) density values
##' @author Marius Hofert
dStudent <- function(x, df, loc = rep(0, d), scale,
                     factor = tryCatch(factorize(scale), error = function(e) e), # needs to be triangular!
                     log = FALSE)
{
    if(!is.matrix(x)) x <- rbind(x)
    n <- nrow(x)
    d <- ncol(x)
    stopifnot(df > 0, length(loc) == d)
    notNA <- apply(!is.na(x), 1, all)
    lres <- rep(-Inf, n)
    lres[!notNA] <- NA
    x <- x[notNA,] # available points
    tx <- t(x) # (d, n)-matrix
    if(inherits(factor, "error") || is.null(factor)) {
        lres[notNA & (colSums(tx == loc) == d)] <- Inf
    } else {
        ## Solve R^T * z = x - mu for z, so z = (R^T)^{-1} * (x - mu) (a (d, d)-matrix)
        ## => z^2 (=> componentwise) = z^T z = (x - mu)^T * ((R^T)^{-1})^T (R^T)^{-1} (x - mu)
        ##                           = z^T z = (x - mu)^T * R^{-1} (R^T)^{-1} (x - mu)
        ##                           = (x - mu)^T * (R^T R)^{-1} * (x - mu)
        ##                           = (x - mu)^T * scale^{-1} * (x - mu) = quadratic form
        z <- backsolve(factor, tx - loc, transpose = TRUE)
        maha <- colSums(z^2) # = sum(z^T z); Mahalanobis distance from x to mu w.r.t. scale
        ## log(sqrt(det(scale))) = log(det(scale))/2 = log(det(R^T R))/2 = log(det(R)^2)/2
        ## = log(prod(diag(R))) = sum(log(diag(R)))
        lrdet <- sum(log(diag(factor)))
        lres[notNA] <- if(is.finite(df)) {
                           df.d.2 <- (df + d) / 2
                           lgamma(df.d.2) - lgamma(df/2) - (d/2) * log(df * pi) - lrdet - df.d.2 * log1p(maha / df)
                       } else {
                           - (d/2) * log(2 * pi) - lrdet - maha/2
                       }
    }
    if(log) lres else exp(lres) # also works with NA, -Inf, Inf
}


### rStudent() #################################################################

##' @title Random Number Generator for a Multivariate Student t Distribution
##' @param n sample size
##' @param df degrees of freedom (positive real or Inf in which case samples
##'        from N(loc, scale) are drawn).
##' @param loc location vector of dimension d
##' @param scale covariance matrix of dimension (d, d)
##' @param factor factorization matrix of the covariance matrix scale; a matrix
##'        R such that R^T R = scale. R is multiplied to the (n, d)-matrix of
##'        independent N(0,1) random variates in the construction from the *right*
##'        (hence the notation).
##' @return (n, d)-matrix with t_nu(loc, scale) samples
##' @author Marius Hofert
rStudent <- function(n, df, loc = rep(0, d), scale, factor = factorize(scale))
{
    ## Checks
    d <- nrow(as.matrix(factor))
    stopifnot(n >= 1, df > 0)
    ## Generate Z ~ N(0, I)
    Z <- matrix(rnorm(n * d), ncol = d) # (n, d)-matrix of N(0, 1)
    ## Generate Y ~ N(0, scale)
    Y <- Z %*% factor # (n, d) %*% (d, k) = (n, k)-matrix of N(0, scale); allows for different k
    ## Generate Y ~ t_nu(0, scale)
    ## Note: sqrt(W) for W ~ df/rchisq(n, df = df) but rchisq() calls rgamma(); see ./src/nmath/rchisq.c
    ##       => W ~ 1/rgamma(n, shape = df/2, rate = df/2)
    if(is.finite(df)) {
        df2 <- df/2
        Y <- Y / rgamma(n, shape = df2, rate = df2) # also fine for different k
    }
    ## Generate X ~ t_nu(loc, scale)
    sweep(Y, 2, loc, "+") # X
}


### pStudent() #################################################################

##' @title Re-order Variables According to their Expected Integration Limits
##'        (Precondition). See [genzbretz2002, p. 957]
##' @param C Cholesky (lower triangular) factor of the covariance matrix R
##' @param q dimension of the problem
##' @param ... all the other parameters same as in pStudent
##' @return list a, b, R, C of integration limits/covariance matrix/...
##'         after reordering has been performed
##' @return d vector indicating the new ordering
##' @author Erik Hintz
##' @note the outermost integral has the smallest limits, etc.
precond_t <- function(a, b, R, C, q, nu)
{
    y <- rep(0, q-1)
    sumysq <- 0
    d <- 1:q

    ## Find i = argmin_j { <expected length of interval> }
    i <- which.min( apply( pt( cbind(a, b) / sqrt(diag(R) ), nu), 1, diff) )

    if(i != 1) {
        ## Swap 1 and i
        tmp <- swap(a,b,R,1,i)
        a <- tmp$a
        b <- tmp$b
        R <- tmp$R
        d[c(1,i)] <- d[c(i,1)]
    }

    ## Store y1
    y[1] <- gamma( (nu+1)/2 )/(gamma( nu/2 ) * sqrt( nu*pi ) ) * nu / (nu+3) *
        ( (1+ a[1]^2/2)^(-(nu+3)/2 ) - (1+b[1]^2/2)^(-(nu+3)/2)) / ( pt(b[1],nu) - pt(a[1],nu) )
    ## Initialize sum of y^2
    sumysq <- y[1]^2

    ## Update Cholesky
    C[1,1] <- sqrt(R[1,1])
    C[2:q,1] <- as.matrix(R[2:q,1]/C[1,1])

    for(j in 2:(q-1)){
        ## c0 <- rowSums(as.matrix(C[j:q,1:(j-1)])^2)
        c01<- sqrt( (nu+j-1)/(nu+sumysq))
        c1 <- c01/sqrt(diag(R)[j:q]-rowSums(as.matrix(C[j:q,1:(j-1)])^2))
        c2 <- as.matrix(C[j:q,1:j-1]) %*%y [1:(j-1)]

        ## Find i = argmin { <expected length of interval j> }
        i <- which.min(pt(c1 * (b[j:q] - c2), nu+j-1) - pt(c1 * (a[j:q] - c2), nu+j-1))+j-1

        if(i!=j){
            ## Swap i and j
            tmp <- swap(a,b,R,i,j)
            a <- tmp$a
            b <- tmp$b
            R <- tmp$R

            C[c(i,j),]    <- as.matrix(C[c(j,i),])
            C[j, (j+1):i] <- as.matrix(0, ncol = i-j, nrow = 1)

            d[c(i,j)] <- d[c(j,i)]
        }
        ## Update Cholesky
        C[j,j] <- sqrt(R[j,j]-sum(C[j,1:(j-1)]^2))
        if(j< (q-1)) C[(j+1):q,j] <- (R[(j+1):q,j]-as.matrix(C[(j+1):q,1:(j-1)])%*%C[j,1:(j-1)])/C[j,j]
        else C[(j+1):q,j] <- (R[(j+1):q,j]-C[(j+1):q,1:(j-1)]%*%C[j,1:(j-1)])/C[j,j]

        ## Get yj
        ajbj <- c01 * c( (a[j]-C[j,1:(j-1)]%*%y[1:(j-1)]),(b[j]-C[j,1:(j-1)]%*%y[1:(j-1)]))/C[j,j]
        y[j] <- exp(lgamma((nu+j)/2)-lgamma((nu+j-1)/2)) / (sqrt((nu+j-1)*pi))*(nu+j-1)/(nu+j+2)*
            ( (1+ajbj[1]^2/2)^(-(nu+j+2)/2) - (1+ajbj[2]^2/2)^(-(nu+j+2)/2))/(pt(ajbj[2],nu+j-1)-pt(ajbj[1],nu+j-1))*
            sqrt( (nu+sumysq)/(nu+j-1))

        sumysq <- sumysq + y[j]^2
    }
    C[q,q] <- sqrt(R[q,q] - sum(C[q,1:(q-1)]^2))
    list(a = a, b = b, R = R, C = C, d = d)
}

##' @title Evaluating the Integrand
##' @param n sample size
##' @param skip skip for sobol
##' @param C Cholesky factor
##' @param q dimension of the problem
##' @param ONE largest number x < 1 such that x != 1
##' @param ZERO smallest number x>0 such that x != 0
##' @param ... all the other parameters same as pStudent
##' @return mean((f(U)+f(1-U))/2) (scalar)
##' @author Erik Hintz
##' @note Switches between f_ and g_ depending on formt
func <- function(n, skip, a, b, C, q, nu, ONE, ZERO, allinfina)
{

    U <- sobol(n = n, d = (q - 1), randomize = TRUE, skip = skip)
    if(allinfina){
        m <- .Call("evalfbonly_",
                   n    = as.integer(n),
                   q    = as.integer(q),
                   U    = as.double(U),
                   b    = as.double(b),
                   C    = as.double(C),
                   nu   = as.double(nu),
                   ONE  = as.double(ONE),
                   ZERO = as.double(ZERO))
    } else {
        m <- .Call("evalf_",
                   n    = as.integer(n),
                   q    = as.integer(q),
                   U    = as.double(U),
                   a    = as.double(a),
                   b    = as.double(b),
                   C    = as.double(C),
                   nu   = as.double(nu),
                   ONE  = as.double(ONE),
                   ZERO = as.double(ZERO))
    }
    return(m)
}

##' @title Distribution Function of the Multivariate t Distribution
##' @param a vector of lower limits
##' @param b vector of upper limits
##' @param R covariance matrix
##' @param nu degrees of freedom
##' @param gam,eps,Nmax,N,n_init tuning parameters. N=number of rep'ns to get sigmahat, algorithm runs unitl gam*sighat < eps or total number evaluation >= Nmax. First loop uses n_init samples in each of the N runs.
##' @param precond logical. If TRUE, preconditioning as described in [genzbretz2002] pp. 955-956 is performed.
##' @author Erik Hintz
pStudent <- function(a, b, R, nu, gam = 3.3, eps = 0.001, Nmax = 1e8, N = 10, n_init = 2^5, precond = TRUE)
{
    if(!is.matrix(R)) R <- as.matrix(R)
    q <- dim(R)[1] # dimension of the problem
    ## Checks
    if( length(a) != length(b) ) stop("Lenghts of a and b differ")
    if( any(a > b) ) stop("a needs to be smaller than b")
    if( q != length(a) ) stop("Dimension of R does not match dimension of a")

    ## Find infinite limits
    infina  <-  (a == -Inf)
    infinb  <-  (b == Inf)
    infinab <-  infina * infinb

    ## Check if all a_i are -Inf
    allinfina <- (sum(infina)==q)

    ## Remove double infinities
    if( sum(infinab) >0 )
    {
        whichdoubleinf <- which( infinab == TRUE)
        a <- a[ -whichdoubleinf ]
        b <- b[ -whichdoubleinf ]
        R <- R[ -whichdoubleinf, -whichdoubleinf ]
        ## Update dimension
        q <- dim(R)[1]
    }

    ## Deal with the univariate case
    if(q == 1) return( pt(b, df = nu) - pt(a, df = nu) )

    ## Get Cholesky factor (lower triangular)
    C <- t(chol(R))

    ## precondtioning (resorting the limits (cf precond_t); only for q > 2)
    if(precond && q>2) {
        temp <- precond_t(a, b, R, C, q, nu)
        a <- temp$a
        b <- temp$b
        C <- temp$C
    }

    gam <- gam / sqrt(N) # instead of dividing sigma by sqrt(N) each time
    n. <- n_init # initial n
    T. <- rep(NA, N) # vector to safe RQMC estimates

    ONE <- 1-.Machine$double.neg.eps
    ZERO <- .Machine$double.eps
    seed <- .Random.seed # need to reset to the seed later when sobol is being used.

    for(l in 1:N)
        T.[l] <- func(n = n., skip = 0, a = a, b = b, C = C, q = q, nu = nu, ONE = ONE, ZERO = ZERO, allinfina = allinfina)
        ## func returns the average of f(u) and f(1-u)
    N. <- 2 * N * n. # N. will count the total number of f. evaluations

    sig <- sd(T.)
    err <- gam * sig # note that this gam is actually gamma/sqrt(N)

    i. <- 0
    while(err > eps && N. < Nmax)
    {
        .Random.seed <- seed # reset seed to have the same shifts in sobol( ... )
        for(l in 1:N)
            T.[l] <- ( T.[l] + func(n = n., skip = n.,  a = a, b = b, C = C, q = q, nu = nu, ONE = ONE, ZERO = ZERO, allinfina = allinfina) )/2
        ## Note that T.[l] and func(...) both depend on n. evaluations; hence we can just average them
        N. <- N. + 2 * N * n.
        n. <- 2 * n.
        sig <- sd(T.)
        err <- gam * sig # note that this gam is actually gamma/sqrt(N)
        i. <- i. + 1
    }
    T <- mean(T.)
    list(Prob = T, N = N., i = i., ErrEst = err)
}

## ##' @title Evaluate t Copulas
## ##' @param u (n, d) matrix of evaluation points
## ##' @param R covariance matrix
## ##' @param nu degrees of freedom
## ##' @return list of four: estimated prob, total number of function evaluations needed, number of iterations
## ##'         in the while loop, estimated error
## ##' @author Erik Hintz and Marius Hofert
## t_cop_prob <- function(u, R, nu)
## {
##     if(!is.matrix(u)) u <- rbind(u)
##     t_prob(rep(-Inf, ncol(u)), b = qt(u, nu), R = R, nu = nu)
## }
