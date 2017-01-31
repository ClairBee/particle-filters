
# SETUP FUNCTIONS                                                                               ####

#' Create parameter list
#' 
#' @export
#' 
create.params <- function(d = 25, seed = 42, 
                          mu.rng = c(1:6), sig.0.par = list("diag" = 2, "od" = 0.4),
                          phi = 0.97, sig.v.par = list("diag" = sample(3:12, d, rep = T)/10, "od" = 0.3),
                          rho = 1, sig.w.par = list("diag" = sample(8:12, d, rep = T)/10, "od" = 0)) {
    
    set.seed(seed)
    
    mu.0 <- sample(mu.rng, d, replace = T)
    phi <- phi
    rho <- rho
    
    sig.0 <- matrix(sig.0.par$od, ncol = d, nrow = d)
    diag(sig.0) <- sig.0.par$diag
    
    sig.v <- matrix(sig.v.par$od, ncol = d, nrow = d)
    diag(sig.v) <- sig.v.par$diag
    
    sig.w <- matrix(sig.w.par$od, ncol = d, nrow = d)
    diag(sig.w) <- sig.w.par$diag
    jitter(rep(0.3, 100))
    
    list(mu.0 = mu.0, sig.0 = sig.0,
         phi = phi, sig.v = sig.v,
         rho = rho, sig.w = sig.w)
}




#' Synthetic state space model
#' 
#' Given a set of parameters, simulate hidden states & observed values of a Gaussian HMM.
#' @details NB in 1 dimension, Sigma is given as a variance, not a standard deviation
#' @export
#' 
synthesise.data <- function(TT = 100, params, seed = 42) {
    
    set.seed(seed)
    
    mu.0 <- params$mu.0; sig.0 <- params$sig.0; phi <- params$phi
    sig.v <- params$sig.v; sig.w <- params$sig.w; rho <- params$rho
    
    d <- length(mu.0)
    x <- y <- array(dim = c(TT, d))
    
    if (d == 1) {
        sig.0 <- matrix(sig.0, ncol = 1)
        sig.v <- matrix(sig.w, ncol = 1)
        sig.w <- matrix(sig.w, ncol = 1)
    }
    
    x[1,] <- rmvnorm(1, mu.0, sig.0, method = "chol")      # stationary distribution
    
    # sample state X_t sequentially
    for(i in 2:TT) {
        x[i,] <- x[i-1,] * phi + rmvnorm(1, rep(0, d), sig.v, method = "chol")
    }
    
    # calculate observations Y_t as batch (dependent only on X_t)
    y <- x * rho + rmvnorm(TT,  rep(0, d), sig.w, method = "chol")
    
    list(x = x, y = y)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# FILTERS                                                                                       ####

#' SIS/SIR with bootstrap proposal
#' 
#' Run Sequential Importance Sampler/Resampler for Gaussian HMM, with known parameters.
#' @details NB in 1 dimension, Sigma is given as a variance, not a standard deviation
#' @export
#' 
bootstrap.SIR <- function(y, param, N = 100, seed = 42, resample = T) {
    
    oldw <- getOption("warn")
    options(warn = -1)
    
    set.seed(seed)
    
    if ((class(y) == "numeric") || (class(y) == "matrix" && (ncol(y) == 1))) {
        
        # UNIVARIATE
        
        y <- c(y)
        TT <- length(y) + 1
        d <- 1
        
        # initialise particles
        particles <- rnorm(N, param$mu.0, sqrt(param$sig.0))
        w <- rep(1/N, N)
        
        ll <- log(mean(w))
        ESS <- 1/sum(w^2)
        
        filt.M <- sum(w * particles)
        filt.V <- sum(w * particles^2) - filt.M^2
        
        for (t in 2:TT) {
            
            if (resample) {
                resamp <- sample(1:N, prob = w, size = N, rep = T)
                particles <- particles[resamp]
            }
            
            particles <- rnorm(N, param$phi * particles, sqrt(param$sig.v))
            
            w.u <- dnorm(y[t-1], param$rho * particles, sqrt(param$sig.w))
            w <- w.u / sum(w.u)
            
            ESS[t] <- 1/sum(w^2)
            ll[t] <- ll[t-1] + log(mean(w.u))
            
            filt.M[t] <- sum(w * particles)
            filt.V[t] <- sum(w * particles^2) - filt.M[t]^2
        }
        
    } else {
        
        # MULTIVARIATE
        
        TT <- nrow(y) + 1
        d <- ncol(y)
        
        filt.M <- array(dim = c(TT, d))
        filt.V <- array(dim = c(TT, d, d))
        
        particles <- rmvnorm(N, param$mu.0, param$sig.0, method = "chol")
        w <- rep(1/N, N)
        
        ll <- log(mean(w))
        ESS <- 1/sum(w^2)
        
        filt.M[1,] <- apply(sweep(particles, 1, w, "*"), 2, sum)
        
        filt.V[1,,] <- apply(sweep(aaply(particles, 1, function(x) x %*% t(x), .drop = F), 1, w, "*"), 2:3, sum) - 
            (filt.M[1,] %*% t(filt.M[1,]))
        
        for (t in 2:TT) {
            
            if (resample) {
                resamp <- sample(1:N, prob = w, size = N, rep = T)
                particles <- particles[resamp,, drop = F]
            }
            
            particles <- aaply(particles, 1, function(x) rmvnorm(1, param$phi * x, param$sig.v, method = "chol"))
            
            w.u <- apply(particles, 1, function(x) dmvnorm(y[t-1,], param$rho * x, param$sig.w))
            w <- w.u / sum(w.u)
            
            ESS[t] <- 1/sum(w^2)
            ll[t] <- ll[t-1] + log(mean(w.u))
            
            filt.M[t,] <- apply(sweep(particles, 1, w, "*"), 2, sum)
            
            filt.V[t,,] <- apply(sweep(aaply(particles, 1, function(x) x %*% t(x)), 1, w, "*"), 2:3, sum) - 
                (filt.M[t,] %*% t(filt.M[t,]))
        }
    }
    
    filt.V[filt.V < 0] <- 0     # regularize negative variances
    options(warn = oldw)        # return warnings to original setting
    
    return(list("mean" = filt.M, "var" = filt.V, "ESS" = ESS, "ll" = ll))
}



#' SIS/SIR with optimal proposal
#' 
#' @export
#'
optimal.SIR <- function(y, param, N = 100, seed = 42, resample = T) {
    
    oldw <- getOption("warn")
    options(warn = -1)
    
    set.seed(seed)
    
    if ((class(y) == "numeric") || (class(y) == "matrix" && (ncol(y) == 1))) {
        
        # UNIVARIATE
        
        y <- c(y)
        TT <- length(y) + 1
        d <- 1
        
        # initialise particles
        particles <- rnorm(N, param$mu.0, sqrt(param$sig.0))
        w <- rep(1/N, N)
        
        ll <- log(mean(w))
        ESS <- 1/sum(w^2)
        filt.M <- sum(w * particles)
        filt.V <- max(sum(w * particles^2) - filt.M^2, 0)
        
        # time-invariant variances, so can calculate offline
        S <- param$rho^2 * param$sig.v + param$sig.w
        Sig <- param$sig.v - param$rho^2 * param$sig.v %*% solve(S) %*% param$sig.v
        
        # run sequential importance sampler
        for (t in 2:TT) {
            
            # optional resampling step
            if (resample) {
                resamp <- sample(1:N, prob = w, size = N, rep = T)
                particles <- particles[resamp]
                w <- rep(1/N, N)
            }
            
            # state prediction
            x.n <- param$phi * particles
            b.n <- param$rho * x.n
            
            # update weights
            lh <- dnorm(y[t-1], b.n, sqrt(S))
            w.u <- w * lh
            w <- w.u / sum (w.u)
            
            # importance sampling
            a.n <- x.n + Sig * t(param$rho) * solve(param$phi) * (y[t-1] - b.n)
            particles <- rnorm(N, a.n, sqrt(Sig))
            
            # calculate output statistics
            ESS[t] <- 1/sum(w^2)
            ll[t] <- ll[t-1] + log(mean(lh))
            filt.M[t] <- sum(w * particles)
            filt.V[t] <- sum(w * particles^2) - filt.M[t]^2
        }
    } else {
        
        # MULTIVARIATE
        
        TT <- nrow(y) + 1
        d <- ncol(y)
        
        filt.M <- array(dim = c(TT, d))
        filt.V <- array(dim = c(TT, d, d))
        
        # reshape rho to diagonal matrix (reflects covariance structure of observation error)
        rho <- diag(param$rho, ncol = d, nrow = d)
        
        # reshape phi to vector (if scalar)
        phi <- rep(param$phi, d)[1:d]
        
        # time-invariant variances, so can calculate offline
        S <- (rho %*% param$sig.v) %*% t(rho) + param$sig.w
        Sig <- param$sig.v - param$sig.v %*% t(rho) %*% solve(S) %*% rho %*% param$sig.v
        
        # initialise filter
        particles <- rmvnorm(N, param$mu.0, param$sig.0, method = "chol")
        w <- rep(1/N, N)
        
        ll <- log(mean(w))
        ESS <- 1/sum(w^2)
        
        filt.M[1,] <- apply(sweep(particles, 1, w , "*"), 2, sum)
        filt.V[1,,] <- apply(sweep(aaply(particles, 1, function(x) x %*% t(x)), 1, w, "*"), 2:3, sum) - 
            filt.M[1,] %*% t(filt.M[1,])   
        
        # run sequential importance sampler
        for (t in 2:TT) {
            
            # optional resample step
            if (resample) {
                resamp <- sample(1:N, prob = w, size = N, rep = T)
                particles <- particles[resamp,]
                w <- rep(1/N, N)
            }
            
            # state prediction
            x.n <- sweep(particles, 2, phi, "*")
            b.n <- aaply(x.n, 1, function(x) rho %*% x)
            
            # update weights
            lh <- aaply(b.n, 1, function(pp) dmvnorm(y[t-1,], pp, S))
            w.u <- w * lh
            w <- w.u / sum (w.u)
            
            # importance sampling
            a.n <- aaply(abind("x" = x.n, "b" = b.n, along = 0), 2, function(pp) {
                pp["x",] + Sig %*% t(rho) %*% solve(param$sig.w) %*% (y[t-1,] - pp["b",])
            })   
            particles <- aaply(a.n, 1, function(a) rmvnorm(1, a, Sig, method = "chol"))
            
            # calculate output statistics
            ESS[t] <- 1/sum(w^2)
            ll[t] <- ll[t-1] + log(mean(w.u))
            
            filt.M[t,] <- apply(sweep(particles, 1, w , "*"), 2, sum)
            filt.V[t,,] <- apply(sweep(aaply(particles, 1, function(x) x %*% t(x)), 1, w, "*"), 2:3, sum) - 
                filt.M[t,] %*% t(filt.M[t,])   
        }
    }
    
    filt.V[filt.V < 0] <- 0     # regularize negative variances
    options(warn = oldw)        # return warnings to original setting
    
    return(list("mean" = filt.M, "var" = filt.V, "ESS" = ESS, "ll" = ll))    
}



#' Auxiliary particle filter with non-optimal proposal
#' 
#' Log-likelihood is not calculated correctly. Need to look at this
#' @details NB in 1 dimension, Sigma is given as a variance, not a standard deviation. 
#' @export
aux.pf <- function(y, param, N = 100, seed = 42) {
    
    oldw <- getOption("warn")
    options(warn = -1)
    
    set.seed(seed)
    
    if ((class(y) == "numeric") || (class(y) == "matrix" && (ncol(y) == 1))) {
        
        # UNIVARIATE
        
        y <- c(y)
        TT <- length(y) + 1
        d <- 1
        
        filt.M <- NA
        filt.V <- NA
        
        ll <- log(mean(rep(1/N, N)))
        ESS <- TT-1
        
        # initialise particles
        particles <- rnorm(N, param$mu.0, sqrt(param$sig.0))
        
        # run importance sampler
        for (t in 2:TT) {
            
            w <- dnorm(y[t-1], param$rho * particles, sqrt(param$sig.w))
            k <- sample(1:N, size = N, replace = TRUE, prob = w)
            
            aux.particles <- rnorm(N, param$phi * particles, sqrt(param$sig.v))
            lh <- dnorm(y[t-1], param$rho * aux.particles, sqrt(param$sig.w))
            w.aux <- lh / w[k]
            
            resamp <- sample(1:N, size = N, replace = T, prob = w.aux)
            
            particles <- aux.particles[resamp]
            
            filt.M[t] <- mean(particles)
            filt.V[t] <- var(particles)
            
            w.norm <- w.aux / sum(w.aux)
            ESS[t] <- 1/sum(w.norm^2)
            ll[t] <- ll[t-1] + log(mean(lh))
        }
    } else {
        
        # MULTIVARIATE
        
        TT <- nrow(y) + 1
        d <- ncol(y)
        
        # reshape rho to diagonal matrix (reflects covariance structure of observation error)
        rho <- diag(param$rho, ncol = d, nrow = d)
        
        # create arrays for summary output statistics
        filt.M <- array(NA, dim = c(TT, d))
        filt.V <- array(NA, dim = c(TT, d, d))
        
        ll <- log(mean(rep(1/N, N)))
        ESS <- 100
        
        # initialise particles
        particles <- rmvnorm(N, param$mu.0, param$sig.0, method = "chol")
        
        # run importance sampler
        for (t in 2:TT) {
            
            w <- apply(particles, 1, function(x) dmvnorm(y[t-1,], param$rho %*% x, param$sig.w))
            k <- sample(1:N, size = N, replace = TRUE, prob = w)
            
            aux.particles <- aaply(particles[k,], 1, function(x) rmvnorm(1, param$phi * x, param$sig.v))
            lh <- apply(aux.particles, 1, function(xx) dmvnorm(y[t-1,], rho %*% xx, param$sig.w))
            w.aux <- lh / w[k]
            
            resamp <- sample(1:N, size = N, replace = T, prob = w.aux)
            
            particles <- aux.particles[resamp,]
            
            filt.M[t,] <- apply(particles, 2, mean)
            filt.V[t,,] <- cov(particles)
            
            w.norm <- w.aux / sum(w.aux)
            ESS[t] <- 1/sum(w.norm^2)
            ll[t] <- ll[t-1] + log(mean(lh))
        }
        
        filt.V[filt.V < 0] <- 0     # regularize negative variances
        options(warn = oldw)        # return warnings to original setting
    }
    return(list("mean" = filt.M, "var" = filt.V, "ESS" = ESS, "ll" = ll)) 
}