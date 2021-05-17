#' Randomly Generate FMP Parameters
#'
#' Generate monotonic polynomial coefficients for user-specified item
#' complexities and prior distributions.
#'
#' @param n_items Number of items for which to simulate item parameters.
#' @param k Either a scalar for the item complexity of all items or a
#'  vector of length n_items if different items have different item
#'  complexities.
#' @param xi_dist Vector of two elements indicating the lower and upper
#'   bounds of the uniform distribution from which to draw xsi parameters.
#' @param omega_dist Vector of two elements indicating the lower and upper
#'   bounds of the uniform distribution from which to draw omega parameters.
#' @param alpha_dist Vector of two elements indicating the lower and upper
#'   bounds of the uniform distribution from which to draw alpha parameters.
#'   Ignored if all k = 0.
#' @param tau_dist Vector of two elements indicating the lower and upper
#'   bounds of the uniform distribution from which to draw tau parameters.
#'   Ignored if all k = 0.
#'
#' @return
#' \item{bmat}{Item parameters in the b parameterization (polynomial
#' coefficients).}
#' \item{greekmat}{Item parameters in the Greek-letter parameterization}
#'
#' @details Randomly generate FMP item parameters for a given k value.
#'
#' @examples
#' ## generate FMP item parameters for 5 items all with k = 2
#' set.seed(2342)
#' pars <- sim_bmat(n_items = 5, k = 2)
#' pars$bmat
#'
#' ## generate FMP item parameters for 5 items with varying k values
#' set.seed(2432)
#' pars <- sim_bmat(n_items = 5, k = c(1, 2, 0, 0, 2))
#' pars$bmat
#'
#' @importFrom stats runif
#' @export

sim_bmat <- function(n_items, k,
                     ncat = 2,
                     xi_dist = c(-1, 1),
                     omega_dist = c(-1, 1),
                     alpha_dist = c(-1, .5),
                     tau_dist = c(-7, -1)){

  maxk <- max(k)
  
  maxncat <- max(ncat)

  bmat <- matrix(0, nrow = n_items, ncol = 2 * maxk + maxncat)

  if (length(k) == 1) k <- rep(k, n_items)
  
  if (length(k) != n_items) stop("k must either have 1 or n_items elements")
  
  if (length(ncat) == 1) ncat <- rep(ncat, n_items)
  
  if (length(ncat) != n_items) stop("ncat must either have 1 or n_items elements")

  # randomly draw xi and omega parameters
  xi <- matrix(runif(n_items * (maxncat - 1), xi_dist[1], xi_dist[2]), ncol = maxncat - 1)
  omega <- runif(n_items, omega_dist[1], omega_dist[2])

  # randomly draw alpha and tau parameters
  alpha <- matrix(0, nrow = n_items, ncol = maxk)
  tau <- matrix(-Inf, nrow = n_items, ncol = maxk)

  for (i in 1:n_items){
    if (ncat[i] < maxncat) xi[i, ncat[i]:ncol(xi)] <- NA
    if (k[i] != 0){
      alpha[i, 1:k[i]] <- runif(k[i], alpha_dist[1], alpha_dist[2])
      tau[i, 1:k[i]] <- runif(k[i], tau_dist[1], tau_dist[2])
      bmat[i, ] <- greek2b(xi = xi[i, ], omega = omega[i],
                           alpha = alpha[i, ], tau = tau[i, ])
    } else bmat[i, 1:maxncat] <- greek2b(xi = xi[i, ], omega = omega[i])
  }

  # bind together greekmat
  greekmat <- matrix(rbind(alpha, tau), nrow = n_items)
  greekmat <- cbind(xi, omega, greekmat)

  # make nice column names
  if (maxk == 0) colnames(greekmat) <- c(paste0("xi", 1:(maxncat - 1)), "omega") else
    colnames(greekmat) <- c(paste0("xi", 1:(maxncat - 1)), "omega", 
                            paste0(rep(c("alpha", "tau"), maxk),
                                   rep(1:maxk, each = 2)))
  if(maxncat > 2) 
    colnames(bmat) <- c(paste0("b0_", 1:(maxncat - 1)), paste0("b", 1:(2*maxk + 1))) else
      colnames(bmat) <- paste0("b", 0:(ncol(bmat) - 1))

  # output both bmat and greekmat
  list(bmat = bmat, greekmat = greekmat)
}
