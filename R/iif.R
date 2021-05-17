#' FMP Item Information Function
#'
#' Find FMP item information
#' for user-supplied item and person parameters.
#'
#' @param theta Vector of latent trait parameters.
#' @param bmat Items x parameters matrix of FMP item parameters
#' (or a vector of FMP item parameters for a single item).
#' @param cvec Optional vector of lower asymptote parameters.
#' If cvec = NULL, then all lower asymptotes set to 0.
#' @param dvec Optional vector of upper asymptote parameters.
#' If dvec = NULL, then all upper asymptotes set to 1.
#'
#' @return Matrix of item information.
#'
#' @examples
#'
#' # plot the IIF for an item with k = 2
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
  
  if(maxncat == 2){
    
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
    for (i in 1:ncol(amat)) amat[, i] <- i * amat[, i]
    
    if(maxncat == 2){
      
      pprime <- theta %*% t(amat) * p[, , 1] * pstar[, , 2]
      
      out <- pprime ^ 2 / (p[, , 2] * p[, , 1])
    } else{
      
    }
  } else{
    
    if (!is.null(cvec) | !is.null(dvec))
      message("Beware! Parameters with asymptote parameters only available with maxncat = 2!")
    
    p <- irf_fmp(theta = theta, bmat = bmat, maxncat = maxncat, 
                 returncat = 0:(maxncat - 1))
    
    amat <- bmat[, -(1:(maxncat - 1)), drop = FALSE]
    
    if(ncol(amat) > 1){ # if maxk > 0
      amat2 <- amat[, -1, drop = FALSE]
      for (i in 1:ncol(amat2)) amat2[, i] <- i * (i + 1) * amat[, i + 1]
      theta3 <- t(sapply(theta, function(x) x ^ (0:(ncol(bmat) - maxncat - 1))))
      partial_m2 <- theta3 %*% t(amat2)
    } else{ # if maxk = 0
      partial_m2 <- matrix(0, nrow = length(theta), ncol = nrow(bmat))
    }
    
    for (i in 1:ncol(amat)) amat[, i] <- i * amat[, i]
    ntheta <- length(theta)
    
    theta2 <- t(sapply(theta, function(x) x ^ (0:(ncol(bmat) - maxncat))))
    if(nrow(theta2) != ntheta) theta2 <- t(theta2)
    theta <- t(sapply(theta, function(x) x ^ (1:(ncol(bmat) - maxncat + 1))))
    if(nrow(theta) != ntheta) theta <- t(theta)
    
    partial_m <-  theta2 %*% t(amat)
    
    
    b0 <- bmat[, 1:(maxncat - 1), drop = FALSE]
    whichbinary <- which(apply(b0, 1, function(x) sum(!is.na(x))) == 1)
    bm <- bmat[, maxncat:(ncol(bmat)), drop = FALSE]
    
    xis <- as.matrix(apply(b0, 1, cumsum))
    
    G <- G1 <- G2 <- array(NA, dim = c(ntheta, nrow(bmat), maxncat))
    G[, , 1] <- 1
    for(i in 2:maxncat){
      G[, , i] <- exp(theta %*% t(bm) * (i - 1) + rep(1, ntheta) %*% t(xis[i - 1, ]))
    }
    G[is.infinite(G)] <- 1e+308
    
    for(i in 1:maxncat){
      G1[, , i] <- (i - 1) * partial_m * G[, , i]
      G2[, , i] <- (i - 1)^2 * partial_m ^ 2 * G[, , i] + 
        (i - 1) * partial_m2 * G[, , i]
    }
    
    g <- apply(G, c(1, 2), sum, na.rm = TRUE)
    g1 <- apply(G1, c(1, 2), sum, na.rm = TRUE)
    g2 <- apply(G2, c(1, 2), sum, na.rm = TRUE)
    
    p1 <- p2 <- array(NA, dim = dim(p))
    
    for(returncat in 0:(maxncat - 1)){
      f <- G[, , returncat + 1]
      f1 <- G1[, , returncat + 1]
      f2 <- G2[, , returncat + 1]

      p1[, , returncat + 1] <- (f1 * g - f * g1) / g^2
      p2[, , returncat + 1] <- f2 / g - (2 *f1 * g1 + f * g2) / g^2 + 2 * f * g1^2 / g^3
    }
    
    out <- apply(p1^2 / p - p2, c(1, 2), sum, na.rm = TRUE)
  }
  
  out
}