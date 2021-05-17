#' Latent Trait Estimation
#'
#' Compute latent trait estimates using either maximum likelihood (ML) or
#' expected a posteriori (EAP) trait estimation.
#'
#' @name th_est_ml
#' @param dat Data matrix of binary item responses with one column for each
#' item. Alternatively, a vector of binary item responses for one person.
#' @param bmat Matrix of FMP item parameters, one row for each item.
#' @param cvec Vector of lower asymptote parameters, one element for each item.
#' @param dvec Vector of upper asymptote parameters, one element for each item.
#' @param lb Lower bound at which to truncate ML estimates.
#' @param ub Upper bound at which to truncate ML estimates.
#' @param int Matrix with two columns used for numerical integration in EAP.
#' Column 1 contains the x coordinates and Column 2 contains the densities.
#'
#' @return Matrix with two columns: est and either sem or psd
#'    \item{est}{Latent trait estimate}
#'    \item{sem}{Standard error of measurement (mle estimates)}
#'    \item{psd}{Posterior standard deviation (eap estimates)}
#'
#' @examples
#'
#' set.seed(3453)
#' bmat <- sim_bmat(n_items = 20, k = 0)$bmat
#'
#' theta <- rnorm(10)
#' dat <- sim_data(bmat = bmat, theta = theta)
#'
#' ## mle estimates
#' mles <- th_est_ml(dat = dat, bmat = bmat)
#'
#' ## eap estimates
#' eaps <- th_est_eap(dat = dat, bmat = bmat)
#'
#' cor(mles[,1], eaps[,1])
#' # 0.9967317
#'
#' @importFrom stats optimize
#'
#' @export

th_est_ml <- function(dat, bmat, maxncat = 2,
                      cvec = NULL, dvec = NULL,
                      lb = -4, ub = 4) {

    # log likelihood function
    loglik <- function(theta, resp, maxncat, bmat, cvec, dvec) {
        P <- irf_fmp(theta = theta, bmat = bmat, 
                     maxncat = maxncat,
                     cvec = cvec, dvec = dvec,
                     returncat = 0:(maxncat-1))
        logP <- sapply(1:dim(P)[2], function(i) log(P[, i, resp[i] + 1]))
        sum(logP)
    }

    # find mles
    mles <- apply(dat, 1, function(resp)
      optimize(f = loglik, interval = c(lb, ub), resp = resp,
               bmat = bmat, maxncat = maxncat, cvec = cvec,
               dvec = dvec, maximum = TRUE)$maximum)

    # find standard errors
    infos <- iif_fmp(theta = mles, bmat = bmat, maxncat = maxncat,
                     cvec = cvec, dvec = dvec)

    sems <- 1 / sqrt(rowSums(infos))

    # prepare output
    mles <- cbind(mles, sems)

    # nice column names
    colnames(mles) <- c("est", "sem")

    # output estimates and standard errors
    mles
}

#' @rdname th_est_ml
#' @export

th_est_eap <- function(dat, bmat, maxncat = 2,
                       cvec = NULL, dvec = NULL,
                       int = int_mat(npts = 33)) {

  P <- irf_fmp(theta = int[, 1], bmat = bmat,
               maxncat = maxncat, cvec = cvec, 
               dvec = dvec, returncat = 0:(maxncat-1))

  # compute eaps and posterior standard deviations
  eaps <- apply(dat, 1, function(resp) {

    # dat <- matrix(dat, ncol = length(dat), nrow = nrow(p), byrow = TRUE)
    
    lik <- sapply(1:dim(P)[2], function(i) P[, i, resp[i] + 1])
    lik <- apply(lik, 1, prod)

    eap <- sum(int[, 1] * lik * int[, 2]) / sum(lik * int[, 2])

    psd <- sqrt(sum( (int[, 1] - eap) ^ 2 * lik * int[, 2]) /
                  sum(lik * int[, 2]))

    c(eap, psd)
    }
  )

  eaps <- t(eaps)

  # nice column names
  colnames(eaps) <- c("est", "psd")

  # output estimates and standard errors (posterior standard deviations)
  eaps
}
