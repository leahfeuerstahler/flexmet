#' FMP Item Information Function
#'
#' Find FMP item information
#' for user-supplied item and person parameters.
#'
#' @param theta Vector of latent trait parameters.
#' @param bmat Items x parameters matrix of FMP item parameters (or a vector of
#' FMP item parameters for a single item).
#' @param maxncat Maximum number of response categories (the first maxncat - 1
#' columns of bmat are intercepts).
#' @param cvec Optional vector of lower asymptote parameters. If cvec = NULL,
#' then all lower asymptotes set to 0.
#' @param dvec Optional vector of upper asymptote parameters. If dvec = NULL,
#' then all upper asymptotes set to 1.
#'
#' @return Matrix of item information.
#'
#' @examples
#'
#' # plot the IIF for a dichotomous item with k = 2
#'
#' set.seed(2342)
#' bmat <- sim_bmat(n_items = 1, k = 2)$bmat
#'
#' theta <- seq(-3, 3, by = .01)
#'
#' information <- iif_fmp(theta = theta, bmat = bmat)
#'
#' plot(theta, information, type = 'l')
#'
#' @export


iif_fmp <- function(theta, bmat, maxncat = 2, cvec = NULL, dvec = NULL) {

  if (!is.matrix(bmat))
    bmat <- as.matrix(bmat)

  if (ncol(bmat) == 1)
    bmat <- t(bmat)

  if (maxncat == 2) {

    p <- irf_fmp(theta = theta, bmat = bmat, maxncat = maxncat,
                 returncat = 0:(maxncat - 1), cvec = cvec, dvec = dvec)


    pstar <- irf_fmp(theta = theta, bmat = bmat, maxncat = maxncat,
                     returncat = 0:(maxncat - 1))

    # compute first derivatives
    theta <- sapply(theta, function(x) x ^ (0:(ncol(bmat) - maxncat)))
    if (is.matrix(theta))
      theta <- t(theta)

    amat <- bmat[, -1]

    if (!is.matrix(amat))
      amat <- matrix(amat)
    if (nrow(amat) != nrow(bmat))
      amat <- t(amat)
    for (i in seq_len(ncol(amat))) amat[, i] <- i * amat[, i]

    if (maxncat == 2) {

      pprime <- theta %*% t(amat) * p[, , 1] * pstar[, , 2]

      out <- pprime ^ 2 / (p[, , 2] * p[, , 1])
    } else{

    }
  } else{

    if (!is.null(cvec) | !is.null(dvec))
      message("Beware! Asymptote parameters only available with maxncat = 2!")

    p <- irf_fmp(theta = theta, bmat = bmat, maxncat = maxncat,
                 returncat = 0:(maxncat - 1))

    amat <- bmat[, - (1:(maxncat - 1)), drop = FALSE]

    if (ncol(amat) > 1) { # if maxk > 0
      amat2 <- amat[, -1, drop = FALSE]
      for (i in seq_len(ncol(amat2))) amat2[, i] <- i * (i + 1) * amat[, i + 1]
      theta3 <- t(sapply(theta, function(x) x ^ (0:(ncol(bmat) - maxncat - 1))))
      partial_m2 <- theta3 %*% t(amat2)
    } else{ # if maxk = 0
      partial_m2 <- matrix(0, nrow = length(theta), ncol = nrow(bmat))
    }

    for (i in seq_len(ncol(amat))) amat[, i] <- i * amat[, i]
    ntheta <- length(theta)

    theta2 <- t(sapply(theta, function(x) x ^ (0:(ncol(bmat) - maxncat))))
    if (nrow(theta2) != ntheta) theta2 <- t(theta2)
    theta <- t(sapply(theta, function(x) x ^ (1:(ncol(bmat) - maxncat + 1))))
    if (nrow(theta) != ntheta) theta <- t(theta)

    partial_m <-  theta2 %*% t(amat)

    b0 <- bmat[, 1:(maxncat - 1), drop = FALSE]
    bm <- bmat[, maxncat:(ncol(bmat)), drop = FALSE]

    xis <- as.matrix(apply(b0, 1, cumsum))

    h <- h1 <- h2 <- array(NA, dim = c(ntheta, nrow(bmat), maxncat))
    h[, , 1] <- 1
    for (i in 2:maxncat) {
      h[, , i] <- exp(theta %*% t(bm) *
                        (i - 1) + rep(1, ntheta) %*% t(xis[i - 1, ]))
    }
    h[is.infinite(h)] <- 1e+308

    for (i in 1:maxncat) {
      h1[, , i] <- (i - 1) * partial_m * h[, , i]
      h2[, , i] <- (i - 1)^2 * partial_m ^ 2 * h[, , i] +
        (i - 1) * partial_m2 * h[, , i]
    }

    g <- apply(h, c(1, 2), sum, na.rm = TRUE)
    g1 <- apply(h1, c(1, 2), sum, na.rm = TRUE)
    g2 <- apply(h2, c(1, 2), sum, na.rm = TRUE)

    p1 <- p2 <- array(NA, dim = dim(p))

    for (returncat in 0:(maxncat - 1)) {
      f <- h[, , returncat + 1]
      f1 <- h1[, , returncat + 1]
      f2 <- h2[, , returncat + 1]

      p1[, , returncat + 1] <- (f1 * g - f * g1) / g^2
      p2[, , returncat + 1] <- f2 / g - (2 * f1 * g1 + f * g2) / g^2 +
        2 * f * g1^2 / g^3
    }

    out <- apply(p1^2 / p - p2, c(1, 2), sum, na.rm = TRUE)
  }

  out
}
