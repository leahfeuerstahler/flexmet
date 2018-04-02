#' FMP Item Response Function
#'
#' Find FMP item response probabilities
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
#' @return Matrix of item response probabilities.
#'
#' @examples
#'
#' # plot the IRF for an item with k = 2
#'
#' set.seed(2342)
#' bmat <- sim_bmat(n_items = 1, k = 2)$bmat
#'
#' theta <- seq(-3, 3, by = .01)
#'
#' probability <- irf_fmp(theta = theta, bmat = bmat)
#'
#' plot(theta, probability, type = 'l')
#'
#' @export


irf_fmp <- function(theta, bmat, cvec = NULL, dvec = NULL) {

    if (!is.matrix(bmat))
        bmat <- as.matrix(bmat)

    if (ncol(bmat) == 1)
        bmat <- t(bmat)

    if (is.null(cvec))
        cvec <- rep(0, nrow(bmat))

    if (is.null(dvec))
        dvec <- rep(1, nrow(bmat))

    cvec <- matrix(cvec, nrow = length(theta), ncol = nrow(bmat),
                   byrow = TRUE)

    dvec <- matrix(dvec, nrow = length(theta), ncol = nrow(bmat),
                   byrow = TRUE)

    theta <- t(sapply(theta, function(x) x ^ (0:(ncol(bmat) - 1))))

    out <- cvec + (dvec - cvec) / (1 + exp(-theta %*% t(bmat)))

    out
}
