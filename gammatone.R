
ERB <- function(fc) {
  f <- fc / 1000
  24.7 * (4.37 * f + 1)
  # 6.23 * f^2 + 93.39 * f + 28.52
}


# len in sec
gammatone <- function(len, srate, fc, phi, n) {
  a <- 1
  b <- 1.019
  tt <- seq(0, len, length = srate * len)
  a * tt^(n - 1) * exp(-2 * pi * b * ERB(fc) * tt) * cos(2 * pi * fc * tt + phi)
}
