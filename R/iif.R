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


iif_fmp <- function(theta, bmat, cvec = NULL, dvec = NULL) {

  if (!is.matrix(bmat))
    bmat <- as.matrix(bmat)

  if (ncol(bmat) == 1)
    bmat <- t(bmat)

  if (is.null(cvec))
    cvec <- rep(0, nrow(bmat))

  if (is.null(dvec))
    dvec <- rep(1, nrow(bmat))

  p <- irf_fmp(theta = theta, bmat = bmat, cvec = cvec)
  q <- 1 - p
  pstar <- irf_fmp(theta = theta, bmat = bmat)

  # compute first derivatives
  theta <- sapply(theta, function(x) x ^ (0:(ncol(bmat) - 2)))
  if (is.matrix(theta))
      theta <- t(theta)

  amat <- bmat[, -1]

  if (!is.matrix(amat))
      amat <- matrix(amat)
  if (nrow(amat) != nrow(bmat))
      amat <- t(amat)
  for (i in 1:ncol(amat)) amat[, i] <- i * amat[, i]

  pprime <- theta %*% t(amat) * q * pstar

  pprime ^ 2 / (p * q)
}
