
# code used to produce output for coursework - particle filters

library("particle.filters")
setwd("~/Documents/PhD/Courses/Computational-Statistics/Cwk/doc/")

mapply(make.save, 5, c(50, 75, 100, 200, 500))
mapply(make.save, c(25, 75, 100, 200), 500)
mapply(make.save, c(5, 25, 37), 1000)

make.save(20,100, s = 24747)

make.save(38, 1000, s = 12)


# redo runtimes after adding extra data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Quickly batch-run filters                                                                     ####

make.save(1,2, keep = T)    # store first version to reuse data & parameters
invisible(mapply(make.save, d = 1, N = c(10, 50, 75, 100, 200)))

make.save(5, 100, keep = T)
mapply(make.save, 5, c(50, 75, 100, 200, 500, 1000))

# remainder not yet run - need to confirm happy with output so far, then run remainder
make.save(2,5, keep = T)
invisible(mapply(make.save, 2, c(10, 50, 75, 100, 200)))

make.save(10,1000, keep = T, s = 24747)
invisible(mapply(make.save, 10, c(50, 75, 100, 200, 500, 1000)))

make.save(50, 5000)
invisible(mapply(make.save, 50, c(250, 1000)))

# often fails at dimension 100
make.save(100, 500, keep = T)
invisible(mapply(make.save, 100, c(1000)))
    # at 100d, often fails on particle filter - 
    # can't generate auxiliary particles based on initial population

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Runtimes                                                                                      ####

rt <- pf.runtimes()
rt[is.na(pf.rmse())] <- NA

# remove univariate case - different function used
rt <- rt[rt$d > 1,]
#rt <- rt[rt$d <= 50,]

D <- 50
matplot(rt$N[rt$d == D], rt[rt$d == D, 3:7], type = "o", col = pf.col(), pch = pf.pch(),
        xlab = "# particles (N)", ylab = "", lwd = 2, lty = pf.lty(),
        xlim = c(0, 1000), ylim = range(0, rt[rt$d == D,3:7]))

fpdf("runtime/d2"); {
    
    matplot(rt$N[rt$d == 2], rt[rt$d == 2, 3:7], type = "o", col = pf.col(), pch = pf.pch(),
            xlab = "# particles (N)", ylab = "", lwd = 2, lty = pf.lty(),
            xlim = c(0, 1000), ylim = range(0, rt[rt$d == 2,3:7]))
    
    legend("topleft", bty = "n", cex = 0.7, lty = pf.lty(), col = pf.col(), legend = pf.list(), lwd = 2)
    
}; dev.off()

fpdf("runtime/d50"); {
    
    matplot(rt$N[rt$d == 50], rt[rt$d == 50, 3:7], type = "o", col = pf.col(), pch = pf.pch(),
            xlab = "# particles (N)", ylab = "", lwd = 2, lty = pf.lty(), 
            xlim = c(0,1000), ylim = range(0, rt[rt$d == 50, 3:7]))
    
    legend("topleft", bty = "n", cex = 0.7, lty = pf.lty(), col = pf.col(), legend = pf.list())
    
}; dev.off()

fpdf("runtime/N100"); {
    
    matplot(rt$d[rt$N == 100], rt[rt$N == 100, 3:7], type = "o", col = pf.col(), pch = pf.pch(),
            xlab = "Dimension (d)", ylab = "", lwd = 2, lty = pf.lty(),
            xlim = c(0,50), ylim = range(0, rt[rt$N == 100, 3:7], na.rm = T))
    legend("topleft", bty = "n", cex = 0.7, lty = pf.lty(), col = pf.col(), legend = pf.list(), lwd = 2)
    
}; dev.off()

fpdf("runtime/N100-log"); {
    
    matplot(rt$d[rt$N == 100], log(rt[rt$N == 100, 3:7]), type = "o", col = pf.col(), pch = pf.pch(),
            xlab = "Dimension (d)", ylab = "Log(t)", lwd = 2, lty = pf.lty(),
            xlim = c(0,100))
    legend("topleft", bty = "n", cex = 0.7, lty = pf.lty(), col = pf.col(), legend = pf.list(), lwd = 2)
    
}; dev.off()


fpdf("runtime/N1000"); {
    
    matplot(rt$d[rt$N == 1000], rt[rt$N == 1000, 3:7], type = "o", col = pf.col(), pch = pf.pch(),
            xlab = "Dimension (d)", ylab = "", lwd = 2, lty = pf.lty(),
            xlim = c(0,50), ylim = range(0, rt[rt$N == 1000, 3:7], na.rm = T))
    legend("topleft", bty = "n", cex = 0.7, lty = pf.lty(), col = pf.col(), legend = pf.list(), lwd = 2)
    
}; dev.off()

fpdf("runtime/N1000-log"); {
    
    matplot(rt$d[rt$N == 1000], log(rt[rt$N == 1000, 3:7]), type = "o", col = pf.col(), pch = pf.pch(),
            xlab = "Dimension (d)", ylab = "Log(t)", lwd = 2, lty = pf.lty(),
            xlim = c(0,50))
    legend("topleft", bty = "n", cex = 0.7, lty = pf.lty(), col = pf.col(), legend = pf.list(), lwd = 2)
    
}; dev.off()


# single plot showing effect of d & N across a single filter 
# may be better displayed in a table?
# still need to sort out colours
{
    matplot(as.integer(dimnames(runtimes)[[2]])[], t(apply(runtimes, 1:2, mean)),
            type = "o", pch = 15, lty = 1, cex = 0.6, col = mapply(adjustcolor, "steelblue", c(1:9)/10),
            xlab = "N. particles",  ylab = "Time elapsed (1/100s)", cex.axis = 0.7, cex.lab = 0.8)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


# Filter means                                                                                  ####

# show closeness of filter means to KF optimum for varying dimensions

read.FO(1, 10)
read.FO(1, 100)
read.FO(10, 100)
read.FO(10, 1000)

fpdf("filt-M/d1-N10"); {
    plot(1:100, fo.1.10$kf$mean[1,], type = "l", ylab = "", xlab = "Time t", ylim = c(-5,12))
    lines(0:100, fo.1.10$filters$boot.sis$mean, col = pf.col(1), lty = pf.lty(1))
    lines(0:100, fo.1.10$filters$boot.sir$mean, col = pf.col(2), lty = pf.lty(2))
    
    lines(0:100, fo.1.10$filters$opt.sis$mean, col = pf.col(3), lty = pf.lty(3))
    lines(0:100, fo.1.10$filters$opt.sir$mean, col = pf.col(4), lty = pf.lty(4))
    
    lines(0:100, fo.1.10$filters$aux$mean, col = pf.col(5), lty = pf.lty(5))
    
    legend("topright", bty = "n", cex = 0.7, lty = c(1, pf.lty()), 
           col = c("black", pf.col()), legend = c("Kalman filter", pf.list()))
}; dev.off()

fpdf("filt-M/d1-N100"); {
    plot(1:100, fo.1.100$kf$mean[1,], type = "l", ylab = "", xlab = "Time t", ylim = c(-5,12))
    lines(0:100, fo.1.100$filters$boot.sis$mean, col = pf.col(1), lty = pf.lty(1))
    lines(0:100, fo.1.100$filters$boot.sir$mean, col = pf.col(2), lty = pf.lty(2))
    
    lines(0:100, fo.1.100$filters$opt.sis$mean, col = pf.col(3), lty = pf.lty(3))
    lines(0:100, fo.1.100$filters$opt.sir$mean, col = pf.col(4), lty = pf.lty(4))
    
    lines(0:100, fo.1.100$filters$aux$mean, col = pf.col(5), lty = pf.lty(5))
    
    legend("topright", bty = "n", cex = 0.7, lty = c(1, pf.lty()), 
           col = c("black", pf.col()), legend = c("Kalman filter", pf.list()))
}; dev.off()

fpdf("filt-M/d50-N100"); {
    plot(1:100, fo.50.100$kf$mean[1,], type = "l", lwd = 2, ylab = "", xlab = "Time t")
    lines(0:100, fo.50.100$filters$boot.sis$mean[,1], col = pf.col(1), lty = pf.lty(1), type = "l")
    lines(0:100, fo.50.100$filters$boot.sir$mean[,1], col = pf.col(2), lty = pf.lty(2))
    
    lines(0:100, fo.50.100$filters$opt.sis$mean[,1], col = pf.col(3), lty = pf.lty(3))
    lines(0:100, fo.50.100$filters$opt.sir$mean[,1], col = pf.col(4), lty = pf.lty(4))
    
    lines(0:100, fo.50.100$filters$aux$mean[,1], col = pf.col(5), lty = pf.lty(5))
    
    legend("topright", bty = "n", cex = 0.7, lty = c(1, pf.lty()), 
           col = c("black", pf.col()), legend = c("Kalman filter", pf.list()))
}; dev.off()

fpdf("filt-M/d50-N1000"); {
    plot(1:100, fo.50.1000$kf$mean[1,], type = "l", lwd = 2, ylab = "", xlab = "Time t")
    lines(0:100, fo.50.1000$filters$boot.sis$mean[,1], col = pf.col(1), lty = pf.lty(1), type = "l")
    lines(0:100, fo.50.1000$filters$boot.sir$mean[,1], col = pf.col(2), lty = pf.lty(2))
    
    lines(0:100, fo.50.1000$filters$opt.sis$mean[,1], col = pf.col(3), lty = pf.lty(3))
    lines(0:100, fo.50.1000$filters$opt.sir$mean[,1], col = pf.col(4), lty = pf.lty(4))
    
    lines(0:100, fo.50.1000$filters$aux$mean[,1], col = pf.col(5), lty = pf.lty(5))
    
    legend("topright", bty = "n", cex = 0.7, lty = c(1, pf.lty()),
           col = c("black", pf.col()), legend = c("Kalman filter", pf.list()))
}; dev.off()

fpdf("filt-M/d10-N100"); {
    plot(1:100, fo.10.100$kf$mean[1,], type = "l", lwd = 2, ylab = "", xlab = "Time t", ylim = c(-6,10))
    lines(0:100, fo.10.100$filters$boot.sis$mean[,1], col = pf.col(1), lty = pf.lty(1), type = "l")
    lines(0:100, fo.10.100$filters$boot.sir$mean[,1], col = pf.col(2), lty = pf.lty(2))
    
    lines(0:100, fo.10.100$filters$opt.sis$mean[,1], col = pf.col(3), lty = pf.lty(3))
    lines(0:100, fo.10.100$filters$opt.sir$mean[,1], col = pf.col(4), lty = pf.lty(4))
    
    lines(0:100, fo.10.100$filters$aux$mean[,1], col = pf.col(5), lty = pf.lty(5))
    
    legend("topright", bty = "n", cex = 0.7, lty = c(1, pf.lty()), lwd = 2,
           col = c("black", pf.col()), legend = c("Kalman filter", pf.list()))
}; dev.off()

fpdf("filt-M/d10-N1000"); {
    plot(1:100, fo.10.1000$kf$mean[1,], type = "l", lwd = 2, ylab = "", xlab = "Time t", ylim = c(-6,10))
    lines(0:100, fo.10.100$filters$boot.sis$mean[,1], col = pf.col(1), lty = pf.lty(1), type = "l")
    lines(0:100, fo.10.1000$filters$boot.sir$mean[,1], col = pf.col(2), lty = pf.lty(2))
    
    lines(0:100, fo.10.1000$filters$opt.sis$mean[,1], col = pf.col(3), lty = pf.lty(3))
    lines(0:100, fo.10.1000$filters$opt.sir$mean[,1], col = pf.col(4), lty = pf.lty(4))
    
    lines(0:100, fo.10.1000$filters$aux$mean[,1], col = pf.col(5), lty = pf.lty(5))
    
    legend("topright", bty = "n", cex = 0.7, lty = c(1, pf.lty()), lwd = 2,
           col = c("black", pf.col()), legend = c("Kalman filter", pf.list()))
}; dev.off()

# Table of RMSE for all dimensions & population sizes tested                                    ####

df.plot.by.d(pf.rmse(), 50)

df.plot.by.d(rmse, 50)
df.plot.by.N(rmse, 100)
df.plot.N.d(rmse)

rmse <- pf.rmse()
rmse <- rmse[rmse$d > 1,]

fpdf("rmse/d1"); {
    matplot(rmse$N[rmse$d == 1], rmse[rmse$d == 1, 3:7], type = "o", col = pf.col(), pch = pf.pch(),
            xlab = "# particles (N)", ylab = "", lty = 1, lwd = 2, ylim = c(0,3))
    
    legend("topright", bty = "n", cex = 0.7, lty = pf.lty(), col = pf.col(), legend = pf.list())
}; dev.off()

fpdf("rmse/d10"); {
    matplot(rmse$N[rmse$d == 10], rmse[rmse$d == 10, 3:7], type = "o", col = pf.col(), pch = pf.pch(), 
            lwd = 2, lty = 1, ylim = c(0,3), xlab = "# particles (N)", ylab = "")
    
    legend("topright", bty = "n", cex = 0.7, lty = pf.lty(), col = pf.col(), legend = pf.list())
}; dev.off()

fpdf("rmse/d50"); {
    matplot(rmse$N[rmse$d == 50], rmse[rmse$d == 50, 3:7], type = "o", col = pf.col(), pch = pf.pch(), 
        lwd = 2, lty = 1, ylim = c(0,3), xlab = "# particles (N)", ylab = "", xlim = c(0, max(rmse$N)))
    
    legend("topright", bty = "n", cex = 0.7, lty = pf.lty(), col = pf.col(), legend = pf.list())
}; dev.off()

rmse <- rmse[rmse$d <= 50,]
fpdf("rmse/N100"); {
    matplot(rmse$d[rmse$N == 100], rmse[rmse$N == 100, 3:7], type = "o", col = pf.col(), pch = pf.pch(), 
            lwd = 2, lty = 1, ylim = c(0,3), xlab = "# dimensions (d)", ylab = "RMSE", xlim = c(0, 50))
    legend("topleft", bty = "n", cex = 0.7, lty = pf.lty(), col = pf.col(), legend = pf.list())
}; dev.off()  

fpdf("rmse/N1000"); {
    matplot(rmse$d[rmse$N == 1000], rmse[rmse$N == 1000, 3:7], type = "o", col = pf.col(), pch = pf.pch(), 
        lwd = 2, lty = 1, ylim = c(0,3), xlab = "# dimensions (d)", ylab = "RMSE", xlim = c(0, 50))
    legend("topleft", bty = "n", cex = 0.7, lty = pf.lty(), col = pf.col(), legend = pf.list())
}; dev.off()

read.FO(50,100)

plot(fo.50.100$kf$mean[1,], type = "l", lwd = 2)
lines(fo.50.100$fiters$opt.sir$mean[,1], lwd = 2, col = "magenta3")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Filter variances                                                                              ####

read.FO(1, 10)
read.FO(1, 100)
read.FO(1, 1000)
read.FO(1, 50)
read.FO(1, 200)
read.FO(1, 500)

read.FO(10, 100)
read.FO(10, 1000)

# show closeness of filter means to KF optimum for varying dimensions

fpdf("filt-V/d1-N10"); {
    plot(fo.1.10$kf$var[1,1,], type = "l", lwd = 2, ylab = "", xlab = "Time t", ylim = c(0,2.5))
    lines(0:100, fo.1.10$filters$boot.sis$var, col = pf.col(1), lty = pf.lty(1))
    lines(0:100, fo.1.10$filters$boot.sir$var, col = pf.col(2), lty = pf.lty(2))
    
    lines(0:100, fo.1.10$filters$opt.sis$var, col = pf.col(3), lty = pf.lty(3))
    lines(0:100, fo.1.10$filters$opt.sir$var, col = pf.col(4), lty = pf.lty(4))
    
    lines(0:100, fo.1.10$filters$aux$var, col = pf.col(5), lty = pf.lty(5))
    
    legend("topleft", bty = "n", cex = 0.7, lty = c(1, pf.lty()), 
           col = c("black", pf.col()), legend = c("Kalman filter", pf.list()))
}; dev.off()

fpdf("filt-V/d1-N100"); {
    plot(fo.1.100$kf$var[1,1,], type = "l", lwd = 2, ylab = "", xlab = "Time t", ylim = c(0,2.5))
    lines(0:100, fo.1.100$filters$boot.sis$var, col = pf.col(1), lty = pf.lty(1))
    lines(0:100, fo.1.100$filters$boot.sir$var, col = pf.col(2), lty = pf.lty(2))
    
    lines(0:100, fo.1.100$filters$opt.sis$var, col = pf.col(3), lty = pf.lty(3))
    lines(0:100, fo.1.100$filters$opt.sir$var, col = pf.col(4), lty = pf.lty(4))
    
    lines(0:100, fo.1.100$filters$aux$var, col = pf.col(5), lty = pf.lty(5))
    
    legend("topright", bty = "n", cex = 0.7, lty = c(1, pf.lty()),
           col = c("black", pf.col()), legend = c("Kalman filter", pf.list()))
}; dev.off()

fpdf("filt-V/d1-N1000"); {
    plot(fo.1.1000$kf$var[1,1,], type = "l", lwd = 2, ylab = "", xlab = "Time t", ylim = c(0,2.5))
    lines(0:100, fo.1.1000$filters$boot.sis$var, col = pf.col(1), lty = pf.lty(1))
    lines(0:100, fo.1.1000$filters$boot.sir$var, col = pf.col(2), lty = pf.lty(2))
    
    lines(0:100, fo.1.1000$filters$opt.sis$var, col = pf.col(3), lty = pf.lty(3))
    lines(0:100, fo.1.1000$filters$opt.sir$var, col = pf.col(4), lty = pf.lty(4))
    
    lines(0:100, fo.1.1000$filters$aux$var, col = pf.col(5), lty = pf.lty(5))
    
    legend("topright", bty = "n", cex = 0.7, lty = c(1, pf.lty()),
           col = c("black", pf.col()), legend = c("Kalman filter", pf.list()))
}; dev.off()

fpdf("filt-V/d10-N100"); {
    plot(fo.10.100$kf$var[1,1,], type = "l", lwd = 2, ylab = "", xlab = "Time t", ylim = c(0, 4))
    lines(0:100, fo.10.100$filters$boot.sis$var[,1,1], col = pf.col(1), lty = pf.lty(1))
    lines(0:100, fo.10.100$filters$boot.sir$var[,1,1], col = pf.col(2), lty = pf.lty(2))
    
    lines(0:100, fo.10.100$filters$opt.sis$var[,1,1], col = pf.col(3), lty = pf.lty(3))
    lines(0:100, fo.10.100$filters$opt.sir$var[,1,1], col = pf.col(4), lty = pf.lty(4))
    
    lines(0:100, fo.10.100$filters$aux$var[,1,1], col = pf.col(5), lty = pf.lty(5))
    
    legend("topright", bg = "white", box.lwd = 1, cex = 0.7, lty = c(1, pf.lty()), 
           col = c("black", pf.col()), legend = c("Kalman filter", pf.list()))
}; dev.off()

fpdf("filt-V/d10-N1000"); {
    plot(fo.10.1000$kf$var[4,4,], type = "l", lwd = 2, ylab = "", xlab = "Time t", ylim = c(0, 2.5))
    lines(0:100, fo.10.1000$filters$boot.sis$var[,4,4], col = pf.col(1), lty = pf.lty(1))
    lines(0:100, fo.10.1000$filters$boot.sir$var[,4,4], col = pf.col(2), lty = pf.lty(2))
    
    lines(0:100, fo.10.1000$filters$opt.sis$var[,4,4], col = pf.col(3), lty = pf.lty(3))
    lines(0:100, fo.10.1000$filters$opt.sir$var[,4,4], col = pf.col(4), lty = pf.lty(4))
    
    lines(0:100, fo.10.1000$filters$aux$var[,4,4], col = pf.col(5), lty = pf.lty(5))
    
    legend("topright", bg = "white", box.lwd = 1, cex = 0.7, lty = c(1, pf.lty()), 
           col = c("black", pf.col()), legend = c("Kalman filter", pf.list()))
}; dev.off()

fpdf("filt-V/d50-N1000"); {
    plot(fo.50.1000$kf$var[4,4,], type = "l", lwd = 2, ylab = "", xlab = "Time t", ylim = c(0, 2.5))
    lines(0:100, fo.50.1000$filters$boot.sis$var[,4,4], col = pf.col(1), lty = pf.lty(1))
    lines(0:100, fo.50.1000$filters$boot.sir$var[,4,4], col = pf.col(2), lty = pf.lty(2))
    
    lines(0:100, fo.50.1000$filters$opt.sis$var[,4,4], col = pf.col(3), lty = pf.lty(3))
    lines(0:100, fo.50.1000$filters$opt.sir$var[,4,4], col = pf.col(4), lty = pf.lty(4))
    
    lines(0:100, fo.50.1000$filters$aux$var[,4,4], col = pf.col(5), lty = pf.lty(5))
    
    legend("topright", bg = "white", box.lwd = 1, cex = 0.7, lty = c(1, pf.lty()), 
           col = c("black", pf.col()), legend = c("Kalman filter", pf.list()))
}; dev.off()

# rmse of variance                                                                              ####        

var.rmse <- function() {
    
    v.rmse.h <- rbind.fill(lapply(list.files("../obj", "filter-obj-", full.names = T), function(fo) {
        o <- readRDS(fo)
        data.frame(t(unlist(c(d = o$param.in$d, N = o$param.in$N, 
                              sapply(o$filters, function(ff) {
                                  if (length(ff) == 1) {
                                      NA
                                  } else {
                                      if (class(ff$var) == "numeric") {
                                          ff$var <- array(ff$var, dim = c(length(ff$var), 
                                                                          1,1))
                                      }
                                      err <- ff$var[-1,,] - aperm(o$kf$var, 3:1)
                                      sqrt(mean((err[err > 0])^2, 
                                                na.rm = T))
                                  }
                              })))))
    }))
    v.rmse.h <- v.rmse.h[order(v.rmse.h$d, v.rmse.h$N), ]
    
    v.rmse.l <- rbind.fill(lapply(list.files("../obj", "filter-obj-", full.names = T), function(fo) {
        o <- readRDS(fo)
        data.frame(t(unlist(c(d = o$param.in$d, N = o$param.in$N, 
                              sapply(o$filters, function(ff) {
                                  if (length(ff) == 1) {
                                      NA
                                  } else {
                                      if (class(ff$var) == "numeric") {
                                          ff$var <- array(ff$var, dim = c(length(ff$var), 
                                                                          1,1))
                                      }
                                      err <- ff$var[-1,,] - aperm(o$kf$var, 3:1)
                                      sqrt(mean((err[err < 0])^2, 
                                                na.rm = T))
                                  }
                              })))))
    }))
    v.rmse.l <- v.rmse.l[order(v.rmse.l$d, v.rmse.l$N), ]
    list("h" = v.rmse.h, "l" = v.rmse.l)
}

vrms <- var.rmse()
vrms <- vrms[vrms$d > 1,]

fpdf("rmse/d10"); {
    matplot(vrms$h$N[vrms$h$d == 50], vrms$h[vrms$h$d == 50, 3:7], type = "o", col = pf.col(), pch = pf.pch(), 
            lwd = 2, lty = 1, ylim = c(0,3), xlab = "# particles (N)", ylab = "")
    
    legend("topright", bty = "n", cex = 0.7, lty = pf.lty(), col = pf.col(), legend = pf.list())
}; dev.off() 

fpdf("rmse/N1000"); {
    matplot(vrms$d[vrms$N == 1000], vrms[vrms$N == 1000, 3:7], type = "o", col = pf.col(), pch = pf.pch(), 
            lwd = 2, lty = 1, ylim = c(0,5), xlab = "# dimensions (d)", ylab = "RMSE", xlim = c(0, 50))
    legend("topleft", bty = "n", cex = 0.7, lty = pf.lty(), col = pf.col(), legend = pf.list())
}; dev.off()

fpdf("rmse/N1000"); {
    matplot(vrms$h$d[vrms$h$N == 1000], vrms$h[vrms$h$N == 1000, 3:7], type = "o", col = pf.col(), pch = pf.pch(), 
            lwd = 2, lty = 1, ylim = c(0,5), xlab = "# dimensions (d)", ylab = "RMSE", xlim = c(0, 50))
    legend("topleft", bty = "n", cex = 0.7, lty = pf.lty(), col = pf.col(), legend = pf.list())
}; dev.off()

matplot(vrms$h$d, vrms$h[,3:7], pch = pf.pch(), col = pf.col(),
        ylab = "", xlab = "Dimension d")

fpdf("filt-V/err"); {
    matplot(vrms$h$d, vrms$h[,3:7], pch = pf.pch(), col = pf.col(),
            ylab = "", xlab = "Dimension d", 
            ylim = c(-2, 8))
    matplot(vrms$l$d, -vrms$l[,3:7], pch = pf.pch(), col = pf.col(), add = T)
    abline(h = 0, lwd = 1, lty = 2)
    legend("topright", bty = "n", cex = 0.7, pch = pf.pch(), col = pf.col(), legend = pf.list())
}; dev.off()
    
fpdf("filt-V/high-err"); {
    matplot(vrms$h$d, vrms$h[,3:7], pch = pf.pch(), col = pf.col(),
            ylab = "", xlab = "Dimension d")
    legend("topright", bty = "n", cex = 0.7, pch = pf.pch(), col = pf.col(), legend = pf.list())
}; dev.off()

fpdf("filt-V/low-err"); {
    matplot(vrms$l$d, vrms$l[,3:7], pch = pf.pch(), col = pf.col(),
            ylab = "", xlab = "Dimension d")
    legend("topright", bty = "n", cex = 0.7, pch = pf.pch(), col = pf.col(), legend = pf.list())
}; dev.off()



matplot(vrms$l$d, vrms$l[,3:7], pch = pf.pch(), col = pf.col(),
        ylab = "", xlab = "Dimension d")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Effective sample size                                                                         ####

read.FO(1, 10)
read.FO(1, 100)
read.FO(50, 100)
read.FO(10, 1000)

fpdf("ESS/d1-N100"); {

    plot(fo.1.100$filters$boot.sis$ESS, col = pf.col(1), lty = pf.lty(1), type = "l",
         xlab = "t", ylab = "")
    lines(fo.1.100$filters$boot.sir$ESS, col = pf.col(2), lty = pf.lty(2))
    
    lines(fo.1.100$filters$opt.sis$ESS, col = pf.col(3), lty = pf.lty(3))
    lines(fo.1.100$filters$opt.sir$ESS, col = pf.col(4), lty = pf.lty(4))
    
    legend("bottomright", cex = 0.7, lty = pf.lty(), bg = "white", 
           col = pf.col(), legend = pf.list()[1:4])
}; dev.off()

fpdf("ESS/d1-N1000"); {
    
    plot(fo.1.1000$filters$boot.sis$ESS, col = pf.col(1), lty = pf.lty(1), type = "l",
         xlab = "t", ylab = "")
    lines(fo.1.1000$filters$boot.sir$ESS, col = pf.col(2), lty = pf.lty(2))
    
    lines(fo.1.1000$filters$opt.sis$ESS, col = pf.col(3), lty = pf.lty(3))
    lines(fo.1.1000$filters$opt.sir$ESS, col = pf.col(4), lty = pf.lty(4))
    
    legend("bottomright", cex = 0.7, lty = pf.lty(), bg = "white", 
           col = pf.col(), legend = pf.list()[1:4])
}; dev.off()

fpdf("ESS/d10-N1000"); {
    
    plot(fo.10.1000$filters$boot.sis$ESS, col = pf.col(1), lty = pf.lty(1), type = "l",
         xlab = "t", ylab = "")
    lines(fo.10.1000$filters$boot.sir$ESS, col = pf.col(2), lty = pf.lty(2))
    
    lines(fo.10.1000$filters$opt.sis$ESS, col = pf.col(3), lty = pf.lty(3))
    lines(fo.10.1000$filters$opt.sir$ESS, col = pf.col(4), lty = pf.lty(4))
    
    legend("topright", bty = "n", cex = 0.7, lty = pf.lty(), bg = "white", 
           col = pf.col(), legend = pf.list()[1:4])
}; dev.off()

# effective sample size less meaningful for auxiliary particle filter, since
# whole particle population is drawn from required distribution - so is essentially 100

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

