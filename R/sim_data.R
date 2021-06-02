#' Simulate FMP Data
#'
#' Simulate data according to user-specified FMP item parameters and
#' latent trait parameters.
#'
#' @param bmat Matrix of FMP item parameters.
#' @param theta Vector of latent trait values.
#' @param maxncat Maximum number of response categories  (the first maxncat - 1
#' columns of bmat are intercepts)
#' @param cvec Optional vector of lower asymptote parameters. If cvec = NULL,
#' then all lower asymptotes set to 0.
#' @param dvec Optional vector of upper asymptote parameters. If dvec = NULL,
#' then all upper asymptotes set to 1.
#'
#' @return Matrix of randomly generated binary item responses.
#'
#' @examples
#'
#' ## generate 5-category item responses for normally distributed theta
#' ##   and 5 items with k = 2
#'
#' set.seed(2342)
#' bmat <- sim_bmat(n_items = 5, k = 2, ncat = 5)$bmat
#'
#' theta <- rnorm(50)
#' dat <- sim_data(bmat = bmat, theta = theta, maxncat = 5)
#'
#' @importFrom stats runif
#' @export


sim_data <- function(bmat, theta, maxncat = 2, cvec = NULL, dvec = NULL) {

  bmat <- as.matrix(bmat)

  probs <- irf_fmp(theta = theta, bmat = bmat, maxncat = maxncat,
                   returncat = 0:(maxncat - 1), cvec = cvec, dvec = dvec)

  cumprobs <- apply(probs, c(1, 2), cumsum)

  dat <- apply(cumprobs, c(2, 3), function(x) sum(runif(1) > x, na.rm = TRUE))

  dat
}
