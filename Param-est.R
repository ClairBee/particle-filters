
#' Skeleton model for MARSS fitting
#' 
#' @export
#' 
skeleton <- function() {
    list(x0 = "unconstrained",               # mu.0
         V0 = "unconstrained",               # sig.0
         B = "diagonal and equal",           # phi
         U = "zero",                         # state intercept
         Q = "unconstrained",                # sig.v
         Z = "identity",                     # rho
         A = "zero",                         # obs. intercepts
         R = "diagonal and unequal")         # sig.w
}


#' Wrapper function to fit MARSS model and return useful statistics
#' 
#' @export
#' 
fit.marss <- function(y, par, control = list(conv.test.slope.tol = 0.5, abstol = .001, maxit = 500), 
                      method = "kem", silent = T, ...) {
    
    # default control options are based on package defaults
    d <- ncol(y)
    
    if (method == "BFGS") {
        control <- control[!names(control) %in% c("conv.test.slope.tol", "abstol")]
    }
    
    t <- system.time(mm <- MARSS(t(y), model = skeleton(), control = control, 
                                 silent = silent, ...))["elapsed"]
    
    par.est <- list(mu.0 = c(coef(mm)$x0),
                    sig.0 = matrix(coef(mm)$V0[c(1,2,2,3)], d, d),
                    phi = c(coef(mm)$B),
                    sig.v = matrix(coef(mm)$Q[c(1,2,2,3)], d, d),
                    sig.w = diag(c(coef(mm)$R)))
    
    if (method == "BFGS") {
        perf.stats <- data.frame("d" = d, "TT" = nrow(y), "method" = method, 
                                 "runtime" = t, "n.iter" = mm$numIter, "conv" = mm$convergence)
    } else {
        perf.stats <- data.frame("d" = d, "TT" = nrow(y), "method" = method, 
                                            "runtime" = t, "n.iter" = mm$numIter, "conv" = mm$convergence,
                                            "sl.tol" = control$conv.test.slope.tol,
                                            "abs.tol" = control$abstol)
    }
    
    ctrl <- control[!names(control) %in% c("conv.test.slope.tol", "abstol")]
    
    # also mm$iter.record for details
    return(list("par.est" = par.est, "org.par" = par, "perf" = perf.stats,
                "control" = ctrl))
}


#' Top-level wrapper function to generate data, fit MARSS model, and save output
#' 
#' @export
#' 
run.marss <- function(d, hv = T, TT = 100, keep = F, method = "kem", silent = T, 
                      control = list(conv.test.slope.tol = 0.5, abstol = .001, maxit = 500), ...) {
    
    if (hv) {
        # high variances
        param <- create.params(d = d, sig.v.par = list("diag" = sample(3:12, d, rep = T), "od" = 3),
                               sig.w.par = list("diag" = sample(3:12, d, rep = T), "od" = 0))
        v <- "hv"
    } else {
        # mid-low variances (standard)
        param <- create.params(d)
        v <- "lv"
    }
    
    dat <- synthesise.data(100, param)
    
    if (method == "BFGS") {
        control <- control[!names(control) %in% c("conv.test.slope.tol", "abstol")]
    }
    
    t <- system.time(mm <- MARSS(t(dat$y), method = method, model = skeleton(), control = control, 
                                 silent = silent, ...))["elapsed"]
    
    par.est <- list(mu.0 = c(coef(mm)$x0),
                    sig.0 = matrix(coef(mm)$V0[c(1,2,2,3)], d, d),
                    phi = c(coef(mm)$B),
                    sig.v = matrix(coef(mm)$Q[c(1,2,2,3)], d, d),
                    sig.w = diag(c(coef(mm)$R)))
    
    if (method == "BFGS") {
        perf.stats <- data.frame("d" = d, "TT" = TT, "method" = method, "lvl" = v,
                                 "runtime" = t, "n.iter" = mm$numIter, "conv" = mm$convergence)
    } else {
        perf.stats <- data.frame("d" = d, "TT" = TT, "method" = method, "lvl" = v,
                                 "runtime" = t, "n.iter" = mm$numIter, "conv" = mm$convergence,
                                 "sl.tol" = control$conv.test.slope.tol,
                                 "abs.tol" = control$abstol)
    }
    
    ctrl <- control[!names(control) %in% c("conv.test.slope.tol", "abstol")]
    
    # also mm$iter.record for details
    op <- list("par.est" = par.est, "org.par" = par, "perf" = perf.stats,
                "control" = ctrl)

    saveRDS(op, paste0("../obj/p-est/", v, "-", method, "-", d,"-", TT, ".rds"))
    
    if (keep) {
        assign(paste0("mm.", v, ".", method, ".", d,".", TT), op, pos = .GlobalEnv)
    }
}