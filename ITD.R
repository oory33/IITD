setwd("/Users/ryo/Documents/R/R_study/IITD")
source("gammatone.R")
library(tuneR)


# len <- 1
# fc <- 400
# phi <- 0
# n <- 4
# t <- seq(0, len, length = srate * len)

# filtt <- gammatone(len, srate, fc, phi, n)

# filtt <- sigl_f
# filt <- abs(fft(filtt))
# filt <- filt[1:(length(filtt) / 2) + 1]
# f <- seq(1, srate / 2, length = length(filtt) / 2)
# plot(f, filt, type = "l", log = "x", xlim = c(200, 1000))
# png("gammatone.png", width = 600, height = 400)
# par(mgp = c(2.3, 0.7, 0))
# plot(t, filtt, type = "l", xlim = c(0, 0.05), yaxt = "n", ylab = "Amplitude", xlab = "time(sec)", cex.lab = 1.5, cex.axis = 1.3)
# dev.off()

filfcs <- c(450, 600, 750) # gammatoneの中心周波数

setwd("input")
pwd <- sprintf("%s", getwd())
fnames <- list.files(path = pwd, pattern = "*.wav")

for (fname in fnames) {
  dat <- readWave(fname)

  sigl <- dat@left
  sigr <- dat@right
  srate <- dat@samp.rate
  for (filfc in filfcs) {
    winln <- 2 # msec
    interval <- 4 # msec
    limit <- 2.22 # (1 / filfc) / 2 * 1000 # msec
    ylimit <- c(-1000 * (limit), 1000 * limit)

    winspl <- winln / 1000 * srate
    intspl <- interval / 1000 * srate
    limitspl <- round(limit / 1000 * srate)

    window <- function(L) {
      0.54 - 0.46 * cos(2 * pi * c(0:(L - 1)) / (L - 1))
    } # Hamming
    win <- window(winspl)

    filt <- fft(gammatone(length(sigl) / srate, srate, filfc, 0, 4))
    sigl_f <- Re(fft(fft(sigl) * filt, inverse = TRUE))
    sigr_f <- Re(fft(fft(sigr) * filt, inverse = TRUE))


    IITD <- numeric(length(sigl) / intspl)

    i <- limitspl + winspl / 2
    m <- 1
    while (i <= length(sigl_f) - winspl - limitspl + 1) {
      sigl_w <- sigl_f[i:(i + winspl - 1)] * win
      # msl <- mean((sigl_w^2))
      cent <- i + round((length(winspl) / 2))
      list <- seq(cent - limitspl, cent + limitspl, length = limitspl * 2 + 1)
      diflr <- numeric(length(list))
      j <- 1
      while (j <= length(list)) {
        s <- list[j] - round((length(winspl) / 2))
        sigr_w <- sigr_f[s:(s + winspl - 1)] * win
        # msr <- mean(sigr_w^2)
        diflr[j] <- mean((sigr_w - sigl_w)^2) # abs(msl - msr)
        j <- j + 1
      }
      # plot(diflr)
      min <- match(min(diflr), diflr)
      IITD[m] <- (round(length(diflr) / 2) + 1 - min) / srate * 1000 * 1000
      m <- m + 1
      i <- i + intspl
    }
    assign(sprintf("IITD_%s", filfc), IITD)
  }
  t <- seq(0, 8, length = length(IITD))
  setwd("/Users/ryo/Documents/R/R_study/IITD/output")
  png(sprintf("%s.png", substring(fname, 1, (nchar(fname) - 4))), width = 1080, height = 1080)
  par(mgp = c(2.3, 0.7, 0))
  plot(t, IITD_450, type = "p", xlim = c(0, 2), ylab = "ITD(µsec)", xlab = "time(sec)", cex.lab = 1.7, cex.axis = 1.5, ylim = ylimit)
  par(new = T)
  plot(t, IITD_600, type = "p", xlim = c(0, 2), ylab = "", xlab = "", cex.lab = 1.7, cex.axis = 1.5, ylim = ylimit, pch = 17)
  par(new = T)
  plot(t, IITD_750, type = "p", xlim = c(0, 2), ylab = "", xlab = "", cex.lab = 1.7, cex.axis = 1.5, ylim = ylimit, pch = 4)
  dev.off()
  setwd("/Users/ryo/Documents/R/R_study/IITD/input")
}
