#' Find Theta Surrogates
#'
#' Compute surrogate theta values as the
#' normalized first principal component scores.
#'
#' @param dat Matrix of binary item responses.
#'
#' @return Vector of surrogate theta values.
#'
#' @description Compute surrogate theta values as the set of normalized first
#' principal component scores.
#'
#' @examples
#'
#' set.seed(2342)
#' bmat <- sim_bmat(n_items = 5, k = 2)$bmat
#'
#' theta <- rnorm(50)
#' dat <- sim_data(bmat = bmat, theta = theta)
#'
#' tsur <- get_surrogates(dat)
#'
#' @references Liang, L., & Browne, M. W. (2015). A quasi-parametric method for
#' fitting flexible item response functions. \emph{Journal of Educational and
#' Behavioral Statistics}, \emph{40}, 5--34. \doi{10.3102/1076998614556816}
#' @export


get_surrogates <- function(dat) {
  dev_scores <- scale(dat, scale = FALSE)
  svd_dev <- svd(dev_scores)
  pc1 <- svd_dev$u[, 1]

  ## reverse principal component scores if necessary
  if (sum(svd_dev$v[, 1]) < 0)
      pc1 <- -pc1

  qnorm(rank(pc1) / (length(pc1) + 1))
}
