# Trimmed Harrell-Davis median estimator based on the highest density interval of width D
thdme <- function(x, D) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n == 0) return(NA)
  if (n == 1) return(x)
  x <- sort(x)
  a <- (n + 1) / 2; b <- (n + 1) / 2
  hdi <- c(0.5 - D / 2, 0.5 + D / 2)
  hdiCdf <- pbeta(hdi, a, b)
  cdf <- function(xs) {
    xs[xs <= hdi[1]] <- hdi[1]
    xs[xs >= hdi[2]] <- hdi[2]
    (pbeta(xs, a, b) - hdiCdf[1]) / (hdiCdf[2] - hdiCdf[1])
  }
  iL <- floor(hdi[1] * n); iR <- ceiling(hdi[2] * n)
  cdfs <- cdf(iL:iR/n)
  W <- tail(cdfs, -1) - head(cdfs, -1)
  sum(x[(iL + 1):iR] * W)
}

# Standard trimmed Harrell-Davis median estimator
sthdme <- function(x) thdme(x, pnorm(1) - pnorm(-1))

# Optimal trimmed Harrell-Davis median estimator
othdme <- function(x) thdme(x, 0.861678977787423)