#' Polynomial Functions
#'
#' Evaluate a forward or inverse (monotonic) polynomial function.
#'
#' @name inv_poly
#' @param x Scalar polynomial function input.
#' @param y Scalar polynomial function output.
#' @param coefs Vector of coefficients that define a monotonic polynomial,
#' see details.
#' @param lb Lower bound of the search interval.
#' @param ub Upper bound of the search interval.
#'
#' @details
#' \deqn{x = t_0 + t_1y + t_2y^2 + ...}{x = t0 + t1*y + t2*y^2 + ...}
#' Then, for coefs = \eqn{(t_0, t_1, t_2, ...)^\prime}{(t0, t1, t2, ...)'},
#' this function finds the corresponding
#' \eqn{y} value (inv_poly) or \eqn{x} value (fw_poly).
#'
#' @importFrom stats uniroot
#' @export

inv_poly <- function(x, coefs, lb = -1000, ub = 1000) {

    coefs <- as.numeric(coefs)

    func <- function(y, x, coefs) {
        fw_poly(y = y, coefs = coefs) - x
    }

    out <- uniroot(f = func, interval = c(lb, ub), x = x, coefs = coefs)

    y <- out$root

    y

}

#' @rdname inv_poly
#' @export

fw_poly <- function(y, coefs){

  y <- y ^ (0:(length(coefs) - 1))

  t(y) %*% coefs
}



#' Find the Greek-Letter Parameterization corresponding to a b Vector of
#' Item Parameters
#'
#' Convert the b vector of item parameters (polynomial coefficients) to the
#' corresponding Greek-letter parameterization (used to ensure monotonicitiy).
#'
#' @param bvec b vector of item parameters (i.e., polynomial coefficients).
#' @param eps Convergence tolerance.
#'
#' @details See \link{greek2b} for more information about the b (polynomial
#' coefficient) and Greek-letter parameterizations of the FMP model.
#'
#' @return A vector of item parameters in the Greek-letter parameterization.
#'
#'
#' @examples
#'
#' (bvec <- greek2b(xi = 0, omega = 1, alpha = c(.1, .1), tau = c(-2, -2)))
#' ## 0.00000000  2.71828183 -0.54365637  0.29961860 -0.03950623  0.01148330
#'
#' (b2greek(bvec))
#' ##  0.0  1.0  0.1 -2.0  0.1 -2.0
#'
#' @references Liang, L., & Browne, M. W. (2015). A quasi-parametric method for
#' fitting flexible item response functions. \emph{Journal of Educational and
#' Behavioral Statistics}, \emph{40}, 5--34. \doi{10.3102/1076998614556816}
#'
#' @importFrom stats optim
#'
#' @seealso \link{greek2b}
#'
#' @export

b2greek <- function(bvec, eps = 1e-08) {

    # retain only non-zero elements in bvec
    check <- bvec == 0
    while (!all(check)) check <- check[-c(1:2)]
    bvec <- bvec[1:(length(bvec) - length(check))]

    # save the number of zeros
    zeros <- numeric(length(check))


    # find k value
    k <- (length(bvec) - 2) / 2

    # extract xi and omega
    xi <- bvec[1]
    omega <- log(bvec[2])

    # k = 0
    if (k == 0) {
        out <- c(xi, omega)
    } else{

      # k = 1
      if (k == 1) {
        alpha <- -bvec[3] / bvec[2]
        tau <- log(3 * bvec[4] / bvec[2] - alpha ^ 2)
        out <- c(xi, omega, alpha, tau)
      } else{

        # k = 2 (use numerical methods to find bvec)
        if (k > 1) {
          func <- function(alphatau, xi, omega, bvec) {
            alphatau <- matrix(alphatau, ncol = 2, byrow = TRUE)
            alpha <- alphatau[, 1]
            tau <- alphatau[, 2]

            bcand <- greek2b(xi = xi, omega = omega,
                             alpha = alpha, tau = tau)

            # sum of squared differences between target and candidate bvecs
            diffvec <- bcand - bvec
            t(diffvec) %*% diffvec
          }

          startvals <- rep(c(.1, -2), k)

          # minimize func
          out <- optim(startvals, fn = func,
                       method = "BFGS", xi = xi,
                       omega = omega, bvec = bvec)

          # check the solution
          if (out$value < eps)
            out <- c(xi, omega, out$par) else out <- NA
        }
      }
    }
    # output a vector of the same length as the input vector
    c(out, zeros)
}


#' Find the b Vector from a Greek-Letter Parameterization of Item Parameters.
#'
#' Convert the Greek-letter parameterization of item parameters (used to ensure
#' monotonicitiy) to the b-vector parameterization (polynomial coefficients).
#'
#' @param xi see details
#' @param omega see details
#' @param alpha see details, vector of length k, set to NULL if k = 0
#' @param tau see details, vector of length k, set to NULL if k = 0
#'
#' @details   For
#' \deqn{m(\theta) = b_{0} + b_{1}\theta + b_{2}\theta^2 + \cdots +
#' b_{2k+1}\theta^{2k+1}}{m(\theta) = b0 + b1\theta + b2\theta^2 + ... +
#' b(2k+1)\theta^{2k+1}}
#' to be a monotonic function, a necessary and sufficient condition is that its
#' first derivative,
#' \deqn{p(\theta) = a_{0} + a_{1}\theta + ... + a_{2k}\theta^{2k},}{p(\theta)
#' = a0 + a1\theta + ... + a(2k)\theta^{2k},}
#' is nonnegative at all theta. Here, let
#' \deqn{b_{0} = \xi}{b0 = \xi}
#' be the constant of integration and
#' \deqn{b_{s} = a_{s-1}/s}{b(s) = a(s-1)/s}
#' for \eqn{s = 1, 2, ..., 2k+1}.
#' Notice that \eqn{p(\theta)} is a polynomial function of degree \eqn{2k}.
#' A nonnegative polynomial of an even degree can be re-expressed as the
#' product of k quadratic functions.
#'
#' If \eqn{k \geq 1}{k >= 1}:
#' \deqn{p(\theta) =  \exp{\omega} \Pi_{s=1}^{k}[1 - 2\alpha_{s}\theta +
#' (\alpha_{s}^2+ \exp(\tau_{s}))\theta^2]}{p(\theta) = exp{\omega}
#' \Pi_{s=1}^{k}[1 - 2\alpha(s)\theta + (\alpha(s)^2+ exp(\tau(s)))\theta^2]}
#'
#' If \eqn{k = 0}:
#' \deqn{p(\theta) = 0.}
#'
#' @return A vector of item parameters in the b parameterization.
#'
#' @examples
#'
#' (bvec <- greek2b(xi = 0, omega = 1, alpha = .1, tau = -1))
#' ## 0.0000000  2.7182818 -0.2718282  0.3423943
#'
#' (b2greek(bvec))
#' ##  0.0  1.0  0.1 -1.0
#'
#' @references Liang, L., & Browne, M. W. (2015). A quasi-parametric method for
#' fitting flexible item response functions. \emph{Journal of Educational and
#' Behavioral Statistics}, \emph{40}, 5--34. \doi{10.3102/1076998614556816}
#'
#' @seealso \link{b2greek}
#'
#' @export

greek2b <- function(xi, omega, alpha = NULL, tau = NULL) {

    a <- exp(omega)

    # if k > 0, find higher-order polynomial coefficients
    if (!(is.null(alpha) | is.null(tau))) {

        # use T-matrices to find a coefficients
        Tlist <- lapply(1:length(tau), find_t,
                        alpha = alpha, tau = tau)

        for (i in 1:length(Tlist)) {
            a <- Tlist[[i]] %*% a
        }
    }

    a <- as.numeric(a)

    # find b coefficients
    b <- c(xi, a / 1:length(a))

    b
}


#' Transform FMP Item Parameters
#'
#' Given FMP item parameters for a single item and the polynomial coefficients
#' defining a latent trait transformation, find the transformed FMP
#' item parameters.
#'
#' @name transform_b
#' @param bvec Vector of item parameters on the \eqn{\theta} metric: (b0,
#' b1, b2, b3, ...).
#' @param bstarvec Vector of item parameters on the \eqn{\theta^{\star}}{\theta*} metric:
#' (b*0, b*1, b*2, b*3, ...)
#' @param tvec Vector of theta transformation polynomial coefficients: (t0, t1,
#' t2, t3, ...)
#'
#' @return Vector of transformed FMP item parameters.
#'
#' @details Equivalent item response models can be written
#' \deqn{P(\theta) = b_0 + b_1\theta + b_2\theta^2 + \cdots +
#' b_{2k+1}\theta^{2k+1}}{P(\theta) = b0 + b1\theta + b2\theta^2 + ... +
#' b(2k+1)\theta^{2k+1}}
#'
#' and
#'
#' \deqn{P(\theta^\star) = b^\star_0 + b^\star_1\theta^\star +
#' b^\star_2\theta^{\star2}+\cdots + b^\star_{2k^\star+1}\theta^{2k^\star+1}}{
#' P(\theta*) = b*0+b*1\theta* + b*2\theta*^2 + ... + b*_{2k*+1}\theta^{2k*+1}}
#'
#' where
#'
#' \deqn{\theta = t_0 + t_1\theta^\star + t_2\theta^{\star 2} + \cdots +
#' t_{2k_\theta+1}\theta^{\star2k_\theta+1}}{\theta = t0 + t1\theta* +
#' t2\theta*^2 + ... + t_{2k_\theta+1}\theta*^{2k_\theta+1}.}
#'
#' @examples
#'
#' ## example parameters from Table 7 of Reise & Waller (2003)
#' ## goal: transform IRT model to sum score metric
#'
#' a <- c(0.57, 0.68, 0.76, 0.72, 0.69, 0.57, 0.53, 0.64,
#'        0.45, 1.01, 1.05, 0.50, 0.58, 0.58, 0.60, 0.59,
#'        1.03, 0.52, 0.59, 0.99, 0.95, 0.39, 0.50)
#' b <- c(0.87, 1.02, 0.87, 0.81, 0.75, -0.22, 0.14, 0.56,
#'        1.69, 0.37, 0.68, 0.56, 1.70, 1.20, 1.04, 1.69,
#'        0.76, 1.51, 1.89, 1.77, 0.39, 0.08, 2.02)
#'
#' ## convert from difficulties and discriminations to FMP parameters
#'
#' b1 <- 1.702 * a
#' b0 <- - 1.702 * a * b
#' bmat <- cbind(b0, b1)
#'
#' ## theta transformation vector (k_theta = 3)
#' ##  see vignette for details about how to find tvec
#'
#' tvec <- c(-3.80789e+00, 2.14164e+00, -6.47773e-01, 1.17182e-01,
#'           -1.20807e-02, 7.02295e-04, -2.13809e-05, 2.65177e-07)
#'
#' ## transform bmat
#' bstarmat <- t(apply(bmat, 1, transform_b, tvec = tvec))
#'
#' ## inspect transformed parameters
#' signif(head(bstarmat), 2)
#'
#' ## plot test response function
#' ##  should be a straight line if transformation worked
#'
#' curve(rowSums(irf_fmp(x, bmat = bstarmat)), xlim = c(0, 23),
#'       ylim = c(0, 23), xlab = expression(paste(theta,"*")),
#'       ylab = "Expected Sum Score")
#' abline(0, 1, col = 2)
#'
#' @export

transform_b <- function(bvec, tvec) {

    k <- (length(bvec) - 2) / 2

    w <- find_w(tvec = tvec, k = k)
    bvec <- as.numeric(bvec)
    bstarvec <- w %*% bvec

    bstarvec
}

#' @rdname transform_b
#' @export

inv_transform_b <- function(bstarvec, tvec) {

    kstar <- (length(bstarvec) - 2) / 2
    ktheta <- (length(tvec) - 2) / 2

    k <- (kstar - ktheta) / (2 * ktheta + 1)

    stopifnot(k %% 1 == 0)

    w <- find_w(tvec = tvec, k = k)

    bvec <- solve(t(w) %*% w) %*% t(w) %*% bstarvec

    bvec
}
