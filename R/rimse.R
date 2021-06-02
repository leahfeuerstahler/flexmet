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
#' @param ncat Number of response categories (first ncat - 1 elemnts of bvec1
#' and bvec2 are intercepts)
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
#' @return Root integrated mean squared difference between two IRFs
#' (dichotomous items) or expected item scores (polytomous items).
#'
#' @examples
#'
#' set.seed(2342)
#' bmat <- sim_bmat(n_items = 2, k = 2, ncat = c(2, 5))$bmat
#'
#' theta <- rnorm(500)
#' dat <- sim_data(bmat = bmat, theta = theta, maxncat = 5)
#'
#' # k = 0
#' fmp0a <- fmp_1(dat = dat[, 1], k = 0, tsur = theta)
#' fmp0b <- fmp_1(dat = dat[, 2], k = 0, tsur = theta)
#'
#' # k = 1
#' fmp1a <- fmp_1(dat = dat[, 1], k = 1, tsur = theta)
#' fmp1b <- fmp_1(dat = dat[, 2], k = 1, tsur = theta)
#'
#' ## compare estimated curves to the data-generating curve
#' rimse(fmp0a$bmat, bmat[1, -c(2:4)])
#' rimse(fmp1a$bmat, bmat[1, -c(2:4)])
#' rimse(fmp0b$bmat, bmat[2, ], ncat = 5)
#' rimse(fmp1b$bmat, bmat[2, ], ncat = 5)
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

    # compute the root integrated mean squared difference b/n two IRFs or TRFs
    sqrt(sum((irf1 - irf2) ^ 2 * int[, 2]))
}
