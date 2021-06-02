#' FMP Item Response Function
#'
#' Find FMP item response probabilities
#' for user-supplied item and person parameters.
#'
#' @param theta Vector of latent trait parameters.
#' @param bmat Items x parameters matrix of FMP item parameters (or a vector of
#' FMP item parameters for a single item).
#' @param maxncat Maximum number of response categories (the first maxncat - 1
#' columns of bmat are intercepts).
#' @param returncat Response categories for which probabilities should be
#' returned, 0,..., maxncat - 1.
#' @param cvec Optional vector of lower asymptote parameters. If cvec = NULL,
#' then all lower asymptotes set to 0.
#' @param dvec Optional vector of upper asymptote parameters. If dvec = NULL,
#' then all upper asymptotes set to 1.
#'
#' @return Matrix of item response probabilities.
#'
#' @examples
#'
#' # plot the IRF for an item with 4 response categories and k = 2
#'
#' set.seed(2342)
#' bmat <- sim_bmat(n_items = 1, ncat = 4, k = 2)$bmat
#'
#' theta <- seq(-3, 3, by = .01)
#'
#' probability <- irf_fmp(theta = theta, bmat = bmat,
#'                        maxncat = 4, returncat = 0:3)
#' 
#' plot(theta, probability[, , 1], type = 'l', ylab = "probability")
#' points(theta, probability[, , 2], type = 'l')
#' points(theta, probability[, , 3], type = 'l')
#' points(theta, probability[, , 4], type = 'l')
#'
#' @export


irf_fmp <- function(theta, bmat, maxncat = 2, returncat = NA,
                    cvec = NULL, dvec = NULL) {

    if (any(is.na(returncat))) {
        if (maxncat == 2)
            returncat <- 1 else{
                returncat <- 0:(maxncat - 1)
            }
    }

    if (maxncat > 2 & (!is.null(cvec) | !is.null(dvec)))
        message("Beware! Asymptote parameters only available with maxncat = 2!")

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

    ntheta <- length(theta)

    if (maxncat > 2) {

        theta <- t(sapply(theta, function(x)
            x ^ (1:(ncol(bmat) - maxncat + 1))))
        if (nrow(theta) != ntheta) theta <- t(theta)

        b0 <- bmat[, 1:(maxncat - 1), drop = FALSE]
        whichbinary <- which(apply(b0, 1, function(x) sum(!is.na(x))) == 1)
        bm <- bmat[, maxncat:(ncol(bmat)), drop = FALSE]

        xis <- as.matrix(apply(b0, 1, cumsum))

        out <- array(NA, dim = c(ntheta, nrow(bmat), maxncat))
        out[, , 1] <- 1
        for (i in 2:maxncat) {
            out[, , i] <- exp(theta %*% t(bm) * (i - 1) +
                                  rep(1, ntheta) %*% t(xis[i - 1, ]))
        }
        out[is.infinite(out)] <- 1e+308
        den <- apply(out, c(1, 2), sum, na.rm = TRUE)
        for (i in 1:maxncat)
            out[, , i] <- out[, , i] / den

        ## ifelse statement for categories 0 and 1 for binary items
        out[, whichbinary, 2] <- cvec[whichbinary] +
            (dvec[whichbinary] - cvec[whichbinary]) *
            out[, whichbinary, 2]
        out[, whichbinary, 1] <- 1 - out[, whichbinary, 2]

    } else{ # if maxncat = 2
        theta <- t(sapply(theta, function(x) x ^ (0:(ncol(bmat) - 1))))
        out <- array(NA, dim = c(ntheta, nrow(bmat), maxncat))
        out[, , 2] <- cvec + (dvec - cvec) / (1 + exp(-theta %*% t(bmat)))
        out[, , 1] <- 1 - out[, , 2]
    }

    out <- out[, , returncat + 1, drop = FALSE]

    if (dim(out)[3] == 1) out <- matrix(out, nrow = dim(out)[1])

    out
}
