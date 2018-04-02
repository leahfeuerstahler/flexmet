#' Simulate FMP Data
#'
#' Simulate data according to user-specified FMP item parameters and
#' latent trait parameters.
#'
#' @param bmat Matrix of FMP item parameters.
#' @param theta Vector of latent trait values.
#' @param cvec Optional vector of lower asymptote parameters.
#' If cvec = NULL, then all lower asymptotes set to 0.
#' @param dvec Optional vector of upper asymptote parameters.
#' If dvec = NULL, then all upper asymptotes set to 1.
#'
#' @return Matrix of randomly generated binary item responses.
#'
#' @examples
#'
#' ## generate binary item responses for normally distributed theta
#' ##   and 5 items with k = 2
#'
#' set.seed(2342)
#' bmat <- sim_bmat(n_items = 5, k = 2)$bmat
#'
#' theta <- rnorm(50)
#' dat <- sim_data(bmat = bmat, theta = theta)
#'
#' @importFrom stats runif
#' @export


sim_data <- function(bmat, theta, cvec = NULL, dvec = NULL) {

  if (is.null(cvec))
    cvec <- rep(0, nrow(bmat))

  if (is.null(dvec))
    dvec <- rep(1, nrow(bmat))

  bmat <- as.matrix(bmat)
  n_subj <- length(theta)
  probs <- irf_fmp(theta = theta, bmat = bmat,
                   cvec = cvec, dvec = dvec)

  dat <- apply(probs, 2, function(x) as.numeric(runif(n_subj) < x))

  dat
}
