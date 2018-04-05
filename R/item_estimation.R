#' Estimate FMP Item Parameters
#'
#' Estimate FMP item parameters for a single item using user-specified
#' theta values (fixed-effects) using fmp_1, or estimate FMP item parameters
#' for multiple items using fixed-effects or random-effects with fmp.
#'
#' @name fmp
#' @aliases fmp
#' @aliases fmp_1
#'
#' @param dat Vector of 0/1 item responses for N (# subjects) examinees.
#' @param k Vector of item complexities for each item, see details.
#' @param tsur Vector of N (# subjects) surrogate theta values.
#' @param start_vals Start values, For fmp_1, a vector of length 2k+2 in the
#' following order:
#'
#' If k = 0: (xi, omega)
#'
#' If k = 1: (xi, omega, alpha1, tau1)
#'
#' If k = 2: (xi, omega, alpha1, tau1, alpha2, tau2)
#'
#' and so forth. For fmp, add start values for item 1, followed by those for
#' item 2, and so forth. For further help, first fit the model without start
#' values, then inspect the outputted parmat data frame.
#'
#' @param em Logical, use random-effects estimation using the EM algorith? If
#' FALSE, fixed effects estimation is used with theta surrogates.
#' @param eps Covergence tolerance for the EM algorithm. The EM algorithm is
#' said to converge is the maximum absolute difference between parameter
#' estimates for successive iterations is less than eps. Ignored if em = FALSE.
#' @param n_quad Number of quadrature points for EM integration. Ignored if
#' em = FALSE
#' @param max_em Maximum number of EM iterations.
#'
#' @param method Optimization method passed to optim.
#' @param \ldots Additional arguments passed to optim.
#'
#' @return
#' \item{bmat}{Matrix of estimated b-matrix parameters, each row corresponds
#' to an item, and contains b0, b1, ...b(max(k)).}
#' \item{parmat}{Data frame of parameter estimation information, including the
#' Greek-letter parameterization, starting value, and parameter estimate.}
#' \item{k}{Vector of item complexities chosen for each item.}
#' \item{log_lik}{Model log likelihood.}
#' \item{mod}{Optimization information, including output from optim.}
#' \item{AIC}{Model AIC.}
#' \item{BIC}{Model BIC.}
#'
#' @details The FMP item response function for a single item is specified using
#' the composite function,
#' \deqn{P(\theta)=[1+\exp(-m(\theta))]^{-1},}{P(\theta)=[1+exp(-m(\theta))]^{-1},}
#' where \eqn{m(\theta)} is an unbounded and monotonically
#' increasing polynomial function of the latent trait \eqn{\theta}.
#'
#' The item complexity parameter \eqn{k} controls the degree of the polynomial:
#' \deqn{m(\theta)=b_0+b_1\theta+b_2\theta^{2}+...+b_{2k+1}
#' \theta^{2k+1},}{m(\theta)=b0+b1\theta+b2\theta^{2}+...+b(2k+1)\theta^{2k+1},}
#' where \eqn{2k+1} equals the order of the polynomial,
#' \eqn{k} is a nonnegative integer, and
#' \deqn{b=(b0,b1,...,b(2k+1))'}
#' are item parameters that define the location and shape of the IRF. The
#' vector \eqn{b} is called the b-vector parameterization of the FMP Model.
#' When \eqn{k=0}, the FMP IRF equals
#' \deqn{P(\theta)= [1+\exp(-(b_0+b_1\theta))] ^{-1},}{P(\theta)=
#' [1+exp(-(b0+b1\theta))] ^{-1},}
#' and is equivalent to the slope-threshold parameterization of the
#' two-parameter item response model.
#'
#' For \eqn{m(\theta)} to be a monotonic function, the FMP IRF can also be
#' expressed as a function of the vector
#'
#' \deqn{\gamma = (\xi, \omega, \alpha_1, \tau_1, \alpha_2, \tau_2,
#' \cdots \alpha_k,\tau_k)'.}{\gamma = (\xi, \omega, \alpha1, \tau1, \alpha2,
#' \tau2, ... \alpha_k, \tau_k)'.}
#'
#' The \eqn{\gamma} vector is called the Greek-letter parameterization of the
#' FMP model. See Feuerstahler (2016) or Liang & Browne (2015) for details
#' about the relationship between the b-vector and Greek-letter
#' parameterizations.
#'
#'
#' @examples
#'
#' set.seed(2342)
#' bmat <- sim_bmat(n_items = 5, k = 2)$bmat
#'
#' theta <- rnorm(50)
#' dat <- sim_data(bmat = bmat, theta = theta)
#'
#' ## fixed-effects estimation for item 1
#'
#' tsur <- get_surrogates(dat)
#'
#' # k = 0
#' fmp0_it_1 <- fmp_1(dat = dat[, 1], k = 0, tsur = tsur)
#'
#' # k = 1
#' fmp1_it_1 <- fmp_1(dat = dat[, 1], k = 1, tsur = tsur)
#'
#' # k = 2
#' fmp2_it_1 <- fmp_1(dat = dat[, 1], k = 2, tsur = tsur)
#'
#'
#' ## fixed-effects estimation for all items
#'
#' fmp0_fixed <- fmp(dat = dat, k = 0, em = FALSE)
#'
#' ## random-effects estimation for all items
#'
#' fmp0_random <- fmp(dat = dat, k = 0, em = TRUE)
#'
#' @references
#'
#' Elphinstone, C. D. (1983). A target distribution model for nonparametric
#' density estimation. \emph{Communication in Statistics--Theory
#' and Methods}, \emph{12}, 161--198. \doi{10.1080/03610928308828450}
#'
#' Elphinstone, C. D. (1985). \emph{A method of distribution and density
#' estimation} (Unpublished dissertation). University of South Africa,
#' Pretoria, South Africa. \doi{20.500.11892/132832}
#'
#' Falk, C. F., & Cai, L. (2016a). Maximum marginal likelihood estimation of a
#' monotonic polynomial generalized partial credit model with applications to
#' multiple group analysis. \emph{Psychometrika}, \emph{81}, 434--460.
#' \doi{10.1007/s11336-014-9428-7}
#'
#' Falk, C. F., & Cai, L. (2016b). Semiparametric item response functions in
#' the context of guessing. \emph{Journal of Educational Measurement},
#' \emph{53}, 229--247. \doi{10.1111/jedm.12111}
#'
#' Feuerstahler, L. M. (2016). \emph{Exploring alternate latent trait metrics
#' with the filtered monotonic polynomial IRT model} (Unpublished dissertation).
#' University of Minnesota, Minneapolis, MN. \url{http://hdl.handle.net/11299/182267}
#'
#' Liang, L. (2007). \emph{A semi-parametric approach to estimating item
#' response functions} (Unpublished dissertation). The Ohio
#' State University, Columbus, OH. Retrieved from https://etd.ohiolink.edu/
#'
#' Liang, L., & Browne, M. W. (2015). A quasi-parametric method for
#' fitting flexible item response functions. \emph{Journal of Educational
#' and Behavioral Statistics}, \emph{40}, 5--34. \doi{10.3102/1076998614556816}
#'
#' @importFrom stats dnorm qnorm optimize
#' @export


fmp_1 <- function(dat, k, tsur, start_vals = NULL, method = "BFGS", ...){

  parmat <- data.frame("item" = rep(1, 2 * k + 2))

  if (k == 0) parnames <- c("xi", "omega") else
    parnames <- c("xi", "omega",
                  sapply(1:k, function(x)
                    paste(c("alpha", "tau"), x, sep = "")))

  parmat$name <- parnames
  parmat$est <- TRUE


  if (is.null(start_vals)){
    parmat$value[grep("xi", parmat$name)] <- qnorm(mean(dat))
    parmat$value[grep("omega", parmat$name)] <- log(1)
    parmat$value[grepl("alpha", parmat$name) & parmat$est] <- .1
    parmat$value[grepl("tau", parmat$name) & parmat$est] <- log(.1)
  } else
    parmat$value <- start_vals

  ## estimate the model
  mod <- optim(par = parmat$value,
               fn = logl, gr = gr_logl,
               method = "BFGS",
               dat = dat, theta = tsur,
               parmat = parmat,
               control = ...)

  ## save the Greek parameters
  parmat$estimate <- mod$par

  ## find the standard errors
  info <- try(solve(hess_logl(parvec = parmat$estimate,
                                   dat = dat, thetas = tsur, parmat = parmat)))

  parmat$error <- suppressWarnings(try(sqrt(diag(info))))


  xi <- parmat$estimate[grep("xi", parmat$name)]
  omega <- parmat$estimate[grep("omega", parmat$name)]
  alpha <- parmat$estimate[grep("alpha", parmat$name)]
  tau <- parmat$estimate[grep("tau", parmat$name)]


  if (k == 0){
    bmat <- c(xi, exp(omega))
  } else{
    bmat <- greek2b(xi = xi, omega = omega, alpha = alpha, tau = tau)
  }


  out <- list(bmat = bmat, parmat = parmat,
              k = k, log_lik = mod$value, mod = mod,
              AIC = 2 * mod$value + 2 * sum(parmat$est),
              BIC = 2 * mod$value + sum(parmat$est) * log(length(dat)))

  out

}

#' @rdname fmp
#' @export

fmp <- function(dat, k, start_vals = NULL,
                em = TRUE, eps = 1e-04, n_quad = 49,
                method = "BFGS", max_em = 500, ...){

  n_items <- ncol(dat)
  n_subj <- nrow(dat)

  maxk <- max(k) # highest order item

  ### create a matrix of parameter information

  ## index the estimated parameters
  parmat <- data.frame("item" = c(rep(1:n_items, each = 2 * maxk + 2)))

  ## name the estimated parameters

  # item parameters (theta metric)
  if (maxk == 0) parnames <- c("xi", "omega") else
    parnames <- c("xi", "omega",
                  sapply(1:maxk, function(x)
                    paste(c("alpha", "tau"), x, sep = "")))

  parmat$name <- rep(parnames, n_items)

  # indicate whether parameters are estimated
  estmat <- sapply(k, function(x) c(rep(TRUE, 2 + 2 * x),
                                    rep(FALSE, 2 * (maxk - x))))
  parmat$est <- as.logical(estmat)

  # alpha/tau should equal zero/-Inf if not estimated
  # don't need to set alpha = -Inf since all alphas are estimated
  parmat$value <- 0
  parmat$value[grep("tau", parmat$name)] <- -Inf

  # add start values
  if (!is.null(start_vals)) parmat$value[parmat$est] <- start_vals
  if (is.null(start_vals)){
    parmat$value[grep("xi", parmat$name)] <- qnorm(colMeans(dat))
    parmat$value[grep("omega", parmat$name)] <- log(1)
    parmat$value[grepl("alpha", parmat$name) & parmat$est] <- .1
    parmat$value[grepl("tau", parmat$name) & parmat$est] <- log(.1)
  }

  if (!em){

    tsur <- get_surrogates(dat)

    mods <- lapply(1:n_items, function(it){
      parmati <- parmat[parmat$item == it, ]

      optim(par = parmati$value[parmati$est],
            fn = logl, gr = gr_logl,
            method = method,
            dat = dat[, it], thetas = tsur,
            parmat = parmati, ...)
    })

    mod <- list()
    mod$values <- sapply(mods, function(x) x$value)
    mod$value <- sum(sapply(mods, function(x) x$value))
    mod$convergence <- sapply(mods, function(x) x$convergence)
    mod$counts <- t(sapply(mods, function(x) x$counts))

    mod$AICs <- 2 * mod$values + 2 * tapply(parmat$est, parmat$item, sum)
    mod$BICs <- 2 * mod$values +
                tapply(parmat$est, parmat$item, sum) * log(n_subj)

    ## save the Greek parameters
    parmat$estimate <- parmat$value
    parmat$estimate[parmat$est] <- unlist(lapply(mods, function(x) x$par))

    greek_mat <- matrix(parmat$estimate,
                       nrow = n_items, ncol = 2 * maxk + 2,
                       byrow = TRUE)

    if (maxk == 0){
      bmat <- greek_mat
      bmat[, 2] <- exp(bmat[, 2])
    } else{
      bmat <- t(apply(greek_mat, 1, function(x){
        xi <- x[1]
        omega <- x[2]
        alpha <- x[c(1 + 2 * (1:maxk))]
        tau <- x[c(2 + 2 * (1:maxk))]
        greek2b(xi, omega, alpha, tau)
      }))
    }

    ## find the standard errors
    info <- try(solve(hess_logl(parvec = parmat$estimate,
                                dat = dat, thetas = tsur, parmat = parmat)))

    parmat$error <- suppressWarnings(try(sqrt(diag(info))))

  }

  if (em){

    em_out <- em_alg(dat = dat,
                   eps = eps, n_quad = n_quad,
                   method = method, parmat = parmat, max_em = max_em, ...)

    mod <- em_out$mod
    mod$iter <- em_out$iter

    ## save the Greek parameters
    parmat$estimate <- 0
    parmat$estimate[parmat$est] <- mod$par

    greek_mat <- matrix(parmat$estimate,
                       nrow = n_items, ncol = 2 * maxk + 2,
                       byrow = TRUE)

    if (maxk == 0){
      bmat <- greek_mat
      bmat[, 2] <- exp(bmat[, 2])
    } else{
      bmat <- t(apply(greek_mat, 1, function(x){
        xi <- x[1]
        omega <- x[2]
        alpha <- x[c(1 + 2 * (1:maxk))]
        tau <- x[c(2 + 2 * (1:maxk))]
        greek2b(xi, omega, alpha, tau)
      }))
    }
  }


  out <- list(bmat = bmat, parmat = parmat,
              k = k, log_lik = mod$value, mod = mod,
              AIC = 2 * mod$value + 2 * sum(parmat$est),
              BIC = 2 * mod$value + sum(parmat$est) * log(n_subj))


  out

}
