#' Linear and Nonlinear Item Parameter Linking
#'
#' Link two sets of FMP item parameters using linear or nonlinear
#' transformations of the latent trait.
#'
#' @name linking
#' @aliases sl_link
#' @aliases hb_link
#'
#' @param bmat1 FMP item parameters on an anchor test.
#' @param bmat2 FMP item parameters to be rescaled.
#' @param maxncat Maximum number of response categories (the first maxncat - 1
#' columns of bmat1 and bmat2 are intercepts)
#' @param cvec1 Vector of lower asymptote parameters for the anchor test.
#' @param cvec2 Vector of lower asymptote parameters corresponding to the
#' rescaled item parameters.
#' @param dvec1 Vector of upper asymptote parameters for the anchor test.
#' @param dvec2 Vector of upper asymptote parameters corresponding to the
#' rescaled item parameters.
#' @param k_theta Complexity of the latent trait transformation (k_theta = 0 is
#' linear, k_theta > 0 is nonlinear).
#' @param int Matrix with two columns, used for numerical integration. Column 1
#' is a grid of theta values, column 2 are normalized densities associated with
#' the column 1 values.
#' @param \dots Additional arguments passed to optim.
#'
#' @details The goal of item parameter linking is to find a metric
#' transformation such that the fitted parameters for one test can be
#' transformed to the same metric as those for the other test. In the Haebara
#' approach, the overall sum of squared differences between the original and
#' transformed individual item response functions is minimized. In the
#' Stocking-Lord approach, the sum of squared differences between the original
#' and transformed test response functions is minimized. See
#' Feuerstahler (2016, 2019) for details on linking with the FMP model.
#'
#' @return
#' \item{par}{(Greek-letter) parameters estimated by optim.}
#' \item{value}{Value of the minimized criterion function.}
#' \item{counts}{Number of function counts in optim.}
#' \item{convergence}{Convergence criterion given by optim.}
#' \item{message}{Message given by optim.}
#' \item{tvec}{Vector of theta transformation coefficients
#' \eqn{(t = t0,....,t(2k_\theta+1))}}
#' \item{bmat}{Transformed bmat2 item parameters.}
#'
#' @examples
#'
#' set.seed(2342)
#' bmat <- sim_bmat(n_items = 20, k = 2)$bmat
#'
#' theta1 <- rnorm(100)
#' theta2 <- rnorm(100, mean = -1)
#'
#' dat1 <- sim_data(bmat = bmat, theta = theta1)
#' dat2 <- sim_data(bmat = bmat, theta = theta2)
#'
#' # estimate each model with fixed-effects and k = 0
#' fmp0_1 <- fmp(dat = dat1, k = 0, em = FALSE)
#' fmp0_2 <- fmp(dat = dat2, k = 0, em = FALSE)
#'
#' # Stocking-Lord linking
#'
#' \donttest{
#' sl_res <- sl_link(bmat1 = fmp0_1$bmat[1:5, ],
#'                   bmat2 = fmp0_2$bmat[1:5, ],
#'                   k_theta = 0)
#'
#' 
#' hb_res <- hb_link(bmat1 = fmp0_1$bmat[1:5, ],
#'                   bmat2 = fmp0_2$bmat[1:5, ],
#'                   k_theta = 0)
#' }
#'
# # check
# i <- 1 #item index
#
# rimse(fmp0_1$bmat[i, ], fmp0_2$bmat[i, ]) # before transformation
# rimse(fmp0_1$bmat[i, ], sl_res$bmat[i, ]) # after transformation
#
#'
#'
#' @references
#'
#' Feuerstahler, L. M. (2016). \emph{Exploring alternate latent trait metrics
#' with the filtered monotonic polynomial IRT model} (Unpublished dissertation).
#' University of Minnesota, Minneapolis, MN.
#' \url{http://hdl.handle.net/11299/182267}
#'
#' Feuerstahler, L. M. (2019). Metric Transformations and the Filtered
#' Monotonic Polynomial Item Response Model. \emph{Psychometrika}, \emph{84},
#' 105--123. \doi{10.1007/s11336-018-9642-9}
#'
#' Haebara, T. (1980). Equating logistic ability scales by a weighted least
#' squares method. \emph{Japanese Psychological Research}, \emph{22}, 144--149.
#' \doi{10.4992/psycholres1954.22.144}
#'
#' Stocking, M. L., & Lord, F. M. (1983). Developing a common metric in item
#' response theory. \emph{Applied Psychological Measurement}, \emph{7},
#' 201--210. \doi{10.1002/j.2333-8504.1982.tb01311.x}
#'
#' @importFrom stats dnorm optim
#' @export
#'



sl_link <- function(bmat1, bmat2, maxncat = 2,
                    cvec1 = NULL, cvec2 = NULL,
                    dvec1 = NULL, dvec2 = NULL, k_theta,
                    int = int_mat(), ...) {

    # function to find the square root of the sum of squared
    #   differences between FMP IRFs
    rimse_dif <- function(greek_tvec, bmat1, bmat2, maxncat,
                          cvec1, cvec2, dvec1, dvec2, int) {

        # find xi, omega, alpha, and tau from greek_tvec
        xi <- greek_tvec[1]
        omega <- greek_tvec[2]

        if (length(greek_tvec) == 2) {
            alpha <- NULL
            tau <- NULL
        } else {
            alpha <- greek_tvec[seq(3, length(greek_tvec) - 1, by = 2)]
            tau <- greek_tvec[seq(4, length(greek_tvec), by = 2)]
        }
        tvec <- greek2b(xi, omega, alpha, tau)

        theta1 <- int[, 1]

        bmat2star <- t(apply(bmat2, 1, transform_b, tvec = tvec,
                             ncat = maxncat))

        irf1 <- irf_fmp(theta = theta1,
                        bmat = bmat1,
                        maxncat = maxncat,
                        cvec = cvec1,
                        dvec = dvec1,
                        returncat = 0:(maxncat - 1))
        irf1 <- apply(irf1, c(1, 2), function(x)
            sum(x * 0:(maxncat - 1), na.rm = TRUE))

        irf2 <- irf_fmp(theta = theta1,
                        bmat = bmat2star,
                        maxncat = maxncat,
                        cvec = cvec2,
                        dvec = dvec2,
                        returncat = 0:(maxncat - 1))
        irf2 <- apply(irf2, c(1, 2), function(x)
            sum(x * 0:(maxncat - 1), na.rm = TRUE))

        sqrt(sum(int[, 2] * (irf1 - irf2) ^ 2))
    }

    out <- optim(par = c(0, 1, rep(0, 2 * k_theta)), fn = rimse_dif,
                 bmat1 = bmat1, bmat2 = bmat2, maxncat = maxncat,
                 cvec1 = cvec1, cvec2 = cvec2,
                 dvec1 = dvec1, dvec2 = dvec2,
                 int = int, ...)

    greekvec <- out$par
    ifelse(k_theta == 0,

           # if k_theta = 0
           tvec <- greek2b(xi = greekvec[1:(maxncat - 1)],
                           omega = greekvec[maxncat]),

           # if k_theta > 0
           tvec <- greek2b(xi = greekvec[1:(maxncat - 1)],
                           omega = greekvec[maxncat],
                           alpha = greekvec[c(maxncat - 1 + 2 * (1:k_theta))],
                           tau = greekvec[c(maxncat + 2 * (1:k_theta))]))

    out$tvec <- tvec
    out$bmat <- t(apply(bmat2, 1, transform_b, tvec = tvec, ncat = maxncat))
    out
}

#' @rdname linking
#' @export


hb_link <- function(bmat1, bmat2, maxncat = 2,
                    cvec1 = NULL, cvec2 = NULL,
                    dvec1 = NULL, dvec2 = NULL, k_theta,
                    int = int_mat(), ...) {

    # function to find the root integrated mean squared difference
    #   between two FMP TRFs
    rimse_dif <- function(greek_tvec, bmat1, bmat2, maxncat,
                          cvec1, cvec2, dvec1, dvec2, int) {
        xi <- greek_tvec[1]
        omega <- greek_tvec[2]

        if (length(greek_tvec) == 2) {
            alpha <- NULL
            tau <- NULL
        } else {
            alpha <- greek_tvec[seq(3, length(greek_tvec) - 1, by = 2)]
            tau <- greek_tvec[seq(4, length(greek_tvec), by = 2)]
        }
        tvec <- greek2b(xi, omega, alpha, tau)

        theta1 <- int[, 1]

        bmat2star <- t(apply(bmat2, 1, transform_b, tvec = tvec,
                             ncat = maxncat))

        irf1 <- irf_fmp(theta = theta1,
                        bmat = bmat1,
                        maxncat = maxncat,
                        cvec = cvec1,
                        dvec = dvec1,
                        returncat = 0:(maxncat - 1))
        irf1 <- apply(irf1, c(1, 2), function(x)
            sum(x * 0:(maxncat - 1), na.rm = TRUE))

        irf2 <- irf_fmp(theta = theta1,
                        bmat = bmat2star,
                        maxncat = maxncat,
                        cvec = cvec2,
                        dvec = dvec2,
                        returncat = 0:(maxncat - 1))
        irf2 <- apply(irf2, c(1, 2), function(x)
            sum(x * 0:(maxncat - 1), na.rm = TRUE))

        sqrt(sum(int[, 2] * (rowSums(irf1) - rowSums(irf2)) ^ 2))
    }

    out <- optim(par = c(0, 1, rep(0, 2 * k_theta)), fn = rimse_dif,
                 bmat1 = bmat1, bmat2 = bmat2, maxncat = maxncat,
                 cvec1 = cvec1, cvec2 = cvec2,
                 dvec1 = dvec1, dvec2 = dvec2,
                 int = int, ...)

    greekvec <- out$par
    ifelse(k_theta == 0,

           # if k_theta = 0
           tvec <- greek2b(xi = greekvec[1:(maxncat - 1)],
                           omega = greekvec[maxncat]),

           # if k_theta > 0
           tvec <- greek2b(xi = greekvec[1:(maxncat - 1)],
                           omega = greekvec[maxncat],
                           alpha = greekvec[c(maxncat - 1 + 2 * (1:k_theta))],
                           tau = greekvec[c(maxncat + 2 * (1:k_theta))]))

    out$tvec <- tvec
    out$bmat <- t(apply(bmat2, 1, transform_b, tvec = tvec, ncat = maxncat))

    out
}
