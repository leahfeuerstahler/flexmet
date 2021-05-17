#' Root Integrated Mean Squared Difference Between FMP IRFs
#'
#' Compute the root integrated mean squared error (RIMSE) between two FMP IRFs.
#'
#' @param bvec1 Either a vector of FMP item parameters or a function
#' corresponding to a non-FMP IRF. Functions should have exactly one argument,
#' corresponding to the latent trait.
#' @param bvec2 Either a vector of FMP item parameters or a function
#' corresponding to a non-FMP IRF. Functions should have exactly one argument,
#' corresponding to the latent trait.
#' @param c1 Lower asymptote parameter for bvec1.
#' Ignored if bvec1 is a function.
#' @param d1 Upper asymptote parameter for bvec1.
#' Ignored if bvec1 is a function.
#' @param c2 Lower asymptote parameter for bvec2.
#' Ignored if bvec2 is a function.
#' @param d2 Upper asymptote parameter for bvec2.
#' Ignored if bvec2 is a function.
#' @param int Matrix with two columns, used for numerical integration.
#' Column 1 is a grid of theta values, column 2 are normalized densities
#' associated with the column 1 values
#'
#' @return Root integrated mean squared difference between two IRFs.
#'
#' @examples
#'
#' set.seed(2342)
#' bmat <- sim_bmat(n_items = 1, k = 2)$bmat
#'
#' theta <- rnorm(500)
#' dat <- sim_data(bmat = bmat, theta = theta)
#'
#' # k = 0
#' fmp0 <- fmp_1(dat = dat, k = 0, tsur = theta)
#'
#' # k = 1
#' fmp1 <- fmp_1(dat = dat, k = 1, tsur = theta)
#'
#' ## compare estimated curves to the data-generating curve
#' rimse(fmp0$bmat, bmat)
#' rimse(fmp1$bmat, bmat)
#'
#'
#' @references Ramsay, J. O. (1991). Kernel smoothing approaches to
#' nonparametric item characteristic curve estimation. \emph{Psychometrika},
#' \emph{56}, 611--630. \doi{10.1007/BF02294494}
#'
#' @export


rimse <- function(bvec1, bvec2, ncat = 2, 
                  c1 = NULL, d1 = NULL, c2 = NULL, d2 = NULL,
                  int = int_mat()) {

    # find model-predicted probabilities for item 1
    if (is.function(bvec1)) {
      irf1 <- bvec1(int[, 1])
    }  else{
      irf1 <- irf_fmp(theta = int[, 1], bmat = bvec1, maxncat = ncat, 
                      cvec = c1, dvec = d1, returncat = 0:(ncat - 1))
    }
  
    irf1 <- apply(irf1, c(1, 2), function(x) 
      sum(x * 0:(ncat - 1), na.rm = TRUE))

    # find model-predicted probabilities for item 2
    if (is.function(bvec2)) {
        irf2 <- bvec2(int[, 1])
    } else {
        irf2 <- irf_fmp(theta = int[, 1], bmat = bvec2, maxncat = ncat, 
                        cvec = c2, dvec = d2, returncat = 0:(ncat - 1))
    }
    
    irf2 <- apply(irf2, c(1, 2), function(x) 
      sum(x * 0:(ncat - 1), na.rm = TRUE))

    # compute the root integrated mean squared difference between the two IRFs or TRFs
    sqrt(sum( (irf1 - irf2) ^ 2 * int[, 2]))
}
