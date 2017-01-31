
# code used to produce output for coursework - parameter estimation for hidden Markov models

library("particle.filters")
setwd("~/Documents/PhD/Courses/Computational-Statistics/Cwk/doc/")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# estimate parameters                                                                           ####

# options to investigate:
#   method: EM or BFGS (numerical optimisation)
#   magnitude of sig matrices

# useful outputs:
#   time to convergence
#   log-likelihood at convergence?
#   iterations to convergence

# tried to run BFGS from defaults with D = 10, TT = 100, hh. Stopped after 48 minutes, no convergence.

# 2-dimensional, high variances (EM)
D <- 1; x.var <- 10; y.var <- 10; TT <- 100; v <- "hh"; set.seed(101)

run.kem.bfgs <- function(D, max.it = 5000, v = "hh", TT = 100, seed = 101) {
    
    x.var <- switch(substring(v, 1, 1),
                    "h" = 10,
                    "l" = 1)
    y.var <- switch(substring(v, 2,2),
                  "h" = 10,
                  "l" = 1)
    
    # create parameters & synthetic data
    {
        param <- create.params(d = D, 
                               sig.v.par = list("diag" = sample(3:12,D, rep = T) / 10 * x.var,
                                                "od" = 3 / 10 * x.var),
                               sig.w.par = list("diag" = sample(3:12, D, rep = T) / 10 * y.var,
                                                "od" = 0))
        dat <- synthesise.data(TT, param)
    }
    
    # run KEM
    {
        t <- system.time(mm <- MARSS(t(dat$y), model = skeleton(), silent = T,
                                     control = list(maxit = max.it)))["elapsed"]
        
        perf.stats <- data.frame("d" = D, "TT" = nrow(dat$y), "method" = mm$method, lvl = v, 
                                 "runtime" = t, "n.iter" = mm$numIter, ll = mm$logLik,
                                 "conv" = mm$convergence)
        
        saveRDS(list("par.est" = coef(mm), "org.par" = param, "perf" = perf.stats, "ctrl" = mm$control,
                     "trace" = mm$iter.record),
                paste0("../obj/p-est/", v, "-", mm$method, "-", D,"-", TT, ".rds"))
    }
    
    if (mm$conv > 20) {
        cat("KEM numerical error. BFGS not run \n")
        return(NULL)
    }
    # run BFGS using EM as starting values
    {
        t2 <- system.time(mm2 <- MARSS(t(dat$y), model = skeleton(), silent = T, 
                                       method = "BFGS", inits = coef(mm, form = "marss")))["elapsed"]
        
        perf.stats <- data.frame("d" = D, "TT" = nrow(dat$y), "method" = paste0("kem/", mm2$method), lvl = v,
                                 "runtime" = t2, "n.iter" = mm2$numIter, "conv" = mm2$convergence, "ll" = mm2$logLik)
        
        saveRDS(list("par.est" = coef(mm2), "org.par" = param, "perf" = perf.stats,
                     "trace" = mm2$iter.record),
                paste0("../obj/p-est/", v, "-kem", mm2$method, "-", ncol(dat$y),"-", nrow(dat$y), ".rds"))
    }
}

run.bfgs <- function(D, max.it = 5000, v = "hh", TT = 100, seed = 101) {
    
    x.var <- switch(substring(v, 1, 1),
                    "h" = 10,
                    "l" = 1)
    y.var <- switch(substring(v, 2,2),
                    "h" = 10,
                    "l" = 1)
    
    # create parameters & synthetic data
    {
        param <- create.params(d = D, 
                               sig.v.par = list("diag" = sample(3:12,D, rep = T) / 10 * x.var,
                                                "od" = 3 / 10 * x.var),
                               sig.w.par = list("diag" = sample(3:12, D, rep = T) / 10 * y.var,
                                                "od" = 0))
        dat <- synthesise.data(TT, param)
    }
    
    # run BFGS using defaults as starting values (mu0 = mean of y[1])
    {
        t3 <- system.time(mm3 <- MARSS(t(dat$y), model = skeleton(), silent = T, method = "BFGS", 
                                       inits = list(Z=1, B=1, U=0, Q=0.05,
                                                    A=0, R=0.05, x0=mean(dat$y[1,]), V0=1)))["elapsed"]
        
        perf.stats <- data.frame("d" = D, "TT" = nrow(dat$y), "method" = mm3$method, lvl = v,
                                 "runtime" = t3, "n.iter" = mm3$numIter, "conv" = mm3$convergence, ll = mm3$logLik)
        
        saveRDS(list("par.est" = coef(mm3), "org.par" = param, "perf" = perf.stats,
                     "trace" = mm3$iter.record),
                paste0("../obj/p-est/", v, "-", mm3$method, "-", ncol(dat$y),"-", nrow(dat$y), ".rds"))
    }
}

run.kem.bfgs.trunc <- function(D, max.it = 100, v = "hh", TT = 100, seed = 101) {
    
    x.var <- switch(substring(v, 1, 1),
                    "h" = 10,
                    "l" = 1)
    y.var <- switch(substring(v, 2,2),
                    "h" = 10,
                    "l" = 1)
    
    # create parameters & synthetic data
    {
        param <- create.params(d = D, 
                               sig.v.par = list("diag" = sample(3:12,D, rep = T) / 10 * x.var,
                                                "od" = 3 / 10 * x.var),
                               sig.w.par = list("diag" = sample(3:12, D, rep = T) / 10 * y.var,
                                                "od" = 0))
        dat <- synthesise.data(TT, param)
    }
    
    # try BFGS after partial KEM convergence as well
    {
        t.src <- system.time(mm.src <- MARSS(t(dat$y), model = skeleton(), silent = T,
                                             control = list(maxit = max.it)))["elapsed"]
        
        perf.stats <- data.frame("d" = D, "TT" = nrow(dat$y), "method" = paste0(mm.src$method, "-", mm.src$control$maxit), lvl = v,
                                 "runtime" = t.src, "n.iter" = mm.src$numIter, "conv" = mm.src$convergence, "ll" = mm.src$logLik)
        
        saveRDS(list("par.est" = coef(mm.src), "org.par" = param, "perf" = perf.stats,
                     "trace" = mm.src$iter.record),
                paste0("../obj/p-est/", v, "-", mm.src$method, mm.src$control$maxit, "-", ncol(dat$y),"-", nrow(dat$y), ".rds"))
        
        t.bf <- system.time(mm.bf <- MARSS(t(dat$y), model = skeleton(), silent = T, 
                                           method = "BFGS", inits = coef(mm.src, form = "marss")))["elapsed"]
        
        perf.stats <- data.frame("d" = D, "TT" = nrow(dat$y), "method" = paste0("kem-",  mm.src$control$maxit, "/", mm.bf$method), lvl = v,
                                 "runtime" = t.bf, "n.iter" = mm.bf$numIter, "conv" = mm.bf$convergence, "ll" = mm.bf$logLik)
        
        saveRDS(list("par.est" = coef(mm.bf), "org.par" = param, "perf" = perf.stats,
                     "trace" = mm.bf$iter.record),
                paste0("../obj/p-est/", v, "-kem", mm.src$control$maxit,"-", mm.bf$method, "-", ncol(dat$y),"-", nrow(dat$y), ".rds"))
    }
}

mapply(run.kem.bfgs, D = c(1:3), max.it = 5000)
mapply(run.bfgs, D = c(1:3), max.it = 5000)
mapply(run.kem.bfgs.trunc, D = c(1:7), max.it = 20)

run.kem.bfgs(5, max.it = 5000)      # truncated before crash
run.kem.bfgs(6, max.it = 5000)      # truncated before crash
run.kem.bfgs(7, max.it = 5000)      # truncated before crash
run.bfgs(4)

run.kem.bfgs.trunc(10, max.it = 100)


run.kem.bfgs(dd, max.it = 5000)
# run.bfgs(dd)
run.kem.bfgs.trunc(dd, max.it = 20)

run.kem.bfgs(c(1:4), max.it = 5000, TT = 500)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# runtime                                                                                       ####

p.est.runtimes <- function() {
    
    mf <- rbind.fill(lapply(lapply(list.files("../obj/p-est", full.names = T), readRDS), "[[", "perf"))
    mf <- mf[order(mf$d, mf$TT, mf$method),]
    mf$cv <- sapply(mf$conv, function(cc)
        switch(toString(cc), "0" = 1, "1" = 2, "10" = 2, "11" = 2, "52" = 3))
    mf$final <- NA
    mf$final[grep("/BFGS", mf$method)] <- mf$cv[grep("/BFGS", mf$method)]
    mf
}
    
m.cols <- c("black", "blue", "cyan3", "red3", "orange", "forestgreen", "green2")
c.pch <- c(20, 1, 4) 

mf <- p.est.runtimes()
trunc.20 <- ddply(mf[grep("kem-20", mf$method),], .(d), summarise, rt = sum(runtime),
                  cv = min(final, na.rm = T))
trunc.100 <- ddply(mf[grep("kem-100", mf$method),], .(d), summarise, rt = sum(runtime),
                  cv = min(final, na.rm = T))
cbn <- ddply(mf[mf$method %in% c("kem", "kem/BFGS"),], .(d), summarise, runtime = sum(runtime),
                 conv = min(final, na.rm = T))

fpdf("p-est/runtime"); {
    plot(mf$d[mf$method == "BFGS"], mf$runtime[mf$method == "BFGS"] / 60, xlim = range(mf$d),
         col = m.cols[1], pch = c.pch[mf$cv[mf$method == "BFGS"]], ylim = c(0, max(mf$runtime / 60) * 1.1),  
         xlab = "Dimension (d)", ylab = "Time elapsed (m)")
    
    points(mf$d[mf$method == "kem"], mf$runtime[mf$method == "kem"] / 60, col = m.cols[2], 
           pch = c.pch[mf$cv[mf$method == "kem"]])
    
    points(mf$d[mf$method == "kem-20"], mf$runtime[mf$method == "kem-20"] / 60, col = m.cols[4], 
           pch = c.pch[mf$cv[mf$method == "kem-20"]])
    
    points(mf$d[mf$method == "kem-100"], mf$runtime[mf$method == "kem-100"] / 60, col = m.cols[6], 
           pch = c.pch[mf$cv[mf$method == "kem-100"]])
    
    points(cbn$d, cbn$runtime / 60, col = m.cols[3], pch = c.pch[cbn$conv])
    points(trunc.20$d, trunc.20$rt / 60, col = m.cols[5], pch = c.pch[trunc.20$cv])
    points(trunc.100$d, trunc.100$rt / 60, col = m.cols[7], pch = c.pch[trunc.100$cv])
    
    legend("topleft", pch = 20, col = m.cols, cex = 0.7, bty = "n",
           legend = c("BFGS only", "EM", "EM + BFGS", "EM-20", "EM-20 + BFGS", "EM-100", "EM-100 + BFGS"))
}; dev.off()

fpdf("p-est/log-runtime", width = 4); {
    plot(trunc.20$d, log(trunc.20$rt), col = m.cols[5], pch = c.pch[trunc.20$cv],
         xlab = "Dimension (d)", ylab = "Log(time elapsed (s))", ylim = range(-1, log(trunc.100$rt)))
    points(trunc.100$d, log(trunc.100$rt), col = m.cols[7], pch = c.pch[trunc.100$cv])
    
    points(mf$d[mf$method == "kem-20"], log(mf$runtime[mf$method == "kem-20"]), col = m.cols[4], 
           pch = c.pch[mf$cv[mf$method == "kem-20"]])
    
    points(mf$d[mf$method == "kem-100"], log(mf$runtime[mf$method == "kem-100"]), col = m.cols[6], 
           pch = c.pch[mf$cv[mf$method == "kem-100"]])
    
    legend("topleft", pch = 20, col = m.cols[c(5, 7)], cex = 0.7, bty = "n",
           legend = c("EM-20 + BFGS","EM-100 + BFGS"))
}; dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# rate of convergence                                                                           ####
# rate of convergence in all unmodified KEM fittings
cr <- rbind.fill(lapply(list.files("../obj/p-est", full.names = T, pattern = ".+kem-.+"), function(fnm) {
    o <- readRDS(fnm)
    data.frame(t(unlist(c("d" = o$perf$d, "TT" = o$perf$TT, o$trace$logLik))))
}))
matplot(t(cr[,3:12]), type = "o", pch = 20, xlab = "", ylab = "")
# add points to show final log-likelihood.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# accuracy of predictions                                                                       ####

compare.phi <- function() {
    hh <- rbind.fill(lapply(list.files("../obj/p-est", full.names = T), function(fnm) {
        o <- readRDS(fnm)
        if(is.null(o$par.est$phi)) {
            if(is.null(o$par.est$B)) {
                pe <- NA
            } else {
                pe <- o$par.est$B
            }
        } else {
            pe <- o$par.est$phi
        }
        data.frame("d" = o$perf$d, "TT" = o$perf$TT, "method" = o$perf$method, "conv" = o$perf$conv, 
                   "phi.est" = pe, "phi.org" = o$org.par$phi)
    }))
    hh$phi.err <- (hh$phi.est - hh$phi.org)
    hh
}

compare.sig.v <- function() {
    
    ff <- list.files("../obj/p-est", full.names = T)
    ff <- ff[c(grep("kem20-BFGS", ff), grep("kem100-BFGS", ff))]
    hh <- rbind.fill(lapply(ff, function(fnm) {
        o <- readRDS(fnm)
        if(is.null(o$par.est$sig.v)) {
            if(is.null(o$par.est$Q)) {
                sv <- NA
                type <- "NA"
            } else {
                sv <- o$par.est$Q
                if (length(sv) > 1) {
                    v.diag <- sv[gsub("\\(", "", lapply(strsplit(rownames(sv), ","), "[[", 1)) == gsub("\\)", "", lapply(strsplit(rownames(sv), ","), "[[", 2))]
                    v.od <- sv[gsub("\\(", "", lapply(strsplit(rownames(sv), ","), "[[", 1)) != gsub("\\)", "", lapply(strsplit(rownames(sv), ","), "[[", 2))]
                } else {
                    v.diag <- sv
                    v.od <- NA
                }
                type <- "Q"
            }
        } else {
            sv <- o$par.est$phi
            cat(fnm, "\n")
            type <- "V"
        }
        
        diag.err <- v.diag - diag(o$org.par$sig.v)
        od.err <- v.od - o$org.par$sig.v[upper.tri(o$org.par$sig.v)]
        
        suppressWarnings({data.frame("d" = o$perf$d, "TT" = o$perf$TT, "method" = o$perf$method, "conv" = o$perf$conv, 
                   "diag.rmse" = sqrt(mean(diag.err^2)), "od.rmse" = sqrt(mean(od.err^2)),
                   "diag.min.err" = min(diag.err, na.rm = T), "diag.max.err" = max(diag.err, na.rm = T),
                   "od.min.err" = min(od.err, na.rm = T), "od.max.err" = max(od.err, na.rm = T),
                   type)})
    }))

        hh[order(hh$method, as.integer(hh$d), as.integer(hh$TT)),]
}

compare.sig.w <- function() {
    
    ff <- list.files("../obj/p-est", full.names = T)
    ff <- ff[c(grep("kem20-BFGS", ff), grep("kem100-BFGS", ff))]
    hh <- rbind.fill(lapply(ff, function(fnm) {
        o <- readRDS(fnm)
            
        sw <- c(o$par.est$R)
        
        if (length(o$par.est$R) > 1) {
            osw <- diag(o$org.par$sig.w)
        } else {
            osw <- o$org.par$sig.w
        }
        err <- sw - osw
        
        data.frame("d" = o$perf$d, "TT" = o$perf$TT, "method" = o$perf$method, "conv" = o$perf$conv, 
                                     "w.rmse" = sqrt(mean(err^2)))
    }))
    
    hh[order(hh$method, as.integer(hh$d), as.integer(hh$TT)),]
}

compare.V0 <- function() {
    hh <- rbind.fill(lapply(list.files("../obj/p-est", full.names = T), function(fnm) {
        o <- readRDS(fnm)
        
        data.frame("d" = o$perf$d, "TT" = o$perf$TT, "method" = o$perf$method, "conv" = o$perf$conv, 
                   "max.V0.est" = max(o$par.est$V0))
    }))
    hh
}

c.phi <- compare.phi()
plot(c.phi$d, c.phi$phi.err, pch = 20, 
     col = c("black", "forestgreen", "green2", "blue", "red", "orange", "cyan3")[c.phi$method])
abline(h = 0, col = "grey")

ddply(c.phi, .(method), summarise, sqrt(mean(phi.err)^2))   # RMSE
ddply(c.phi, .(method), summarise, mean(phi.err))           # mean error (not absolute)

v.diag <- sv[gsub("\\(", "", lapply(strsplit(rownames(sv), ","), "[[", 1)) == gsub("\\)", "", lapply(strsplit(rownames(sv), ","), "[[", 2))]
v.od <- sv[gsub("\\(", "", lapply(strsplit(rownames(sv), ","), "[[", 1)) != gsub("\\)", "", lapply(strsplit(rownames(sv), ","), "[[", 2))]

c.sig.v <- compare.sig.v()

c.V0 <- compare.V0()
max(c.V0$max.V0.est)

fpdf("p-est/phi-density", width = 4); {
    plot(density(c.phi$phi.est), xlab = "", ylab = "", main = "")
    abline(v = 0.97, col = "red")
}; dev.off()

fpdf("p-est/sig-v-rmse"); {
    plot(c.sig.v$d[2:10], c.sig.v$diag.rmse[2:10], pch = 20, col = "steelblue", 
         ylim = range(0, c.sig.v$diag.rmse[2:10]), xlab = "Dimension (d)", ylab = "RMSE")
    points(c.sig.v$d[2:10], c.sig.v$od.rmse[2:10], pch = 1, col = "steelblue")
    
    points(c.sig.w$d[2:10], c.sig.w$w.rmse[2:10], pch = 20, col = "green3")
    
    legend("bottomleft", bty = "n", cex = 0.7, pch = c(20, 1, 20),
           col = c("steelblue","steelblue", "green3"),
           legend = c(expression(paste(sigma[v], " diagonal")),
                      expression(paste(sigma[v], " off-diagonal")),
                      expression(paste(sigma[w], " diagonal"))))
}; dev.off()

c.sig.w <- compare.sig.w()

plot(c.sig.w$d[2:10], c.sig.w$w.rmse[2:10], pch = 20, col = "steelblue", 
     ylim = range(0, c.sig.w$w.rmse[2:10]), xlab = "Dimension (d)", ylab = "RMSE")

mm <- readRDS("../obj/p-est/hh-kem100-BFGS-4-100.rds")
mm$par
mm$org.par

compare.mu0 <- function() {
    ff <- list.files("../obj/p-est", full.names = T)
    ff <- ff[c(grep("kem20-BFGS", ff), grep("kem100-BFGS", ff))]
    
    rbind.fill(lapply(ff, function(fnm) {
        o <- readRDS(fnm)

        mu.est <- c(o$par.est$x0)
        mu.org <- c(o$org.par$mu.0)
        y1 <- run.dat(length(mu.est))
        
        data.frame("d" = o$perf$d, "TT" = o$perf$TT, "method" = o$perf$method, "conv" = o$perf$conv, 
                   "v.mu0" = sqrt(mean((mu.est - mu.org)^2)),
                   "v.y1" = sqrt(mean((mu.est - y1)^2)))
    }))
}

c.mu <- compare.mu0()

mean(c.mu$v.mu0); mean(c.mu$v.y1)
o$org.par$mu.0
c(o$par.est$x0)

run.dat <- function(D, TT = 100) {
    param <- create.params(d = D, 
                           sig.v.par = list("diag" = sample(3:12,D, rep = T),
                                            "od" = 3),
                           sig.w.par = list("diag" = sample(3:12, D, rep = T),
                                            "od" = 0))
    dat <- synthesise.data(TT, param)
    dat$y[1,]
}