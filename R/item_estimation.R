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
#' @param dat Vector of item responses for N (# subjects) examinees. Binary 
#' data should be coded 0/1, and polytomous data should be coded 0, 1, 2, etc.
#' @param k Vector of item complexities for each item, see details. If
#' k < ncol(dat), k's will be recycled.
#' @param tsur Vector of N (# subjects) surrogate theta values.
#' @param start_vals Start values, For fmp_1, a vector of length 2k+2 in the
#' following order:
#'
#' If k = 0: (xi_1, ..., x_{C_i - 1}, omega)
#'
#' If k = 1: (xi_1, ..., x_{C_i - 1}, omega, alpha1, tau1)
#'
#' If k = 2: (xi_1, ..., x_{C_i - 1}, omega, alpha1, tau1, alpha2, tau2)
#'
#' and so forth. For fmp, add start values for item 1, followed by those for
#' item 2, and so forth. For further help, first fit the model without start
#' values, then inspect the outputted parmat data frame.
#'
#' @param em If "mirt", use the mirt (Chalmers, 2012) package to estimate
#' item parameters. If TRUE, random-effects estimation is used via the EM
#' algorithm. If FALSE, fixed effects estimation is used with theta
#' surrogates.
#' @param eps Covergence tolerance for the EM algorithm. The EM algorithm is
#' said to converge is the maximum absolute difference between parameter
#' estimates for successive iterations is less than eps. Ignored if em = FALSE.
#' @param n_quad Number of quadrature points for EM integration. Ignored if
#' em = FALSE
#' @param max_em Maximum number of EM iterations (for em = TRUE only).
#' @param method Optimization method passed to optim.
#' @param priors List of prior information used to estimate the item parameters.
#' The list should have up to 4 elements named xi, omega, alpha, tau. Each list
#' should be a vector of length 3: the name of the prior distribution ("norm" or
#' "none"), the first parameter of the prior distribution, and the second
#' parameter of the prior distribution. Currently, "norm" and 'none" are the
#' only available prior distributions.
#' @param \ldots Additional arguments passed to optim (if em != "mirt") or mirt
#' (if em == "mirt").
#'
#' @return
#' \item{bmat}{Matrix of estimated b-matrix parameters, each row corresponds
#' to an item, and contains b0, b1, ...b(max(k)).}
#' \item{parmat}{Data frame of parameter estimation information, including the
#' Greek-letter parameterization, starting value, and parameter estimate.}
#' \item{k}{Vector of item complexities chosen for each item.}
#' \item{log_lik}{Model log likelihood.}
#' \item{mod}{If em == "mirt", the mirt object. Otherwise, optimization
#' information, including output from optim.}
#' \item{AIC}{Model AIC.}
#' \item{BIC}{Model BIC.}
#'
#' @details The FMP item response function for a single item \eqn{i} with 
#' responses in categories \eqn{c = 0, ..., C_i - 1} is specified using the 
#' composite function,
#' \deqn{P(X_i = c | \theta) = exp(\sum_{v=0}^c(b_0i_{v}} + m_i(\theta))) / 
#' (\sum_{u=0}^{C_i - 1} exp(\sum_{v=0}^u(b_{0i_{v}} + m_i(\theta)))) }
#' where \eqn{m(\theta)} is an unbounded and monotonically increasing polynomial
#' function of the latent trait \eqn{\theta}, excluding the intercept (s). 
#'
#' The item complexity parameter \eqn{k} controls the degree of the polynomial:
#' \deqn{m(\theta)=b_1\theta+b_2\theta^{2}+...+b_{2k+1}
#' \theta^{2k+1},}{m(\theta)=b1\theta+b2\theta^{2}+...+b(2k+1)\theta^{2k+1},}
#' where \eqn{2k+1} equals the order of the polynomial,
#' \eqn{k} is a nonnegative integer, and
#' \deqn{b=(b1,...,b(2k+1))'}
#' are item parameters that define the location and shape of the IRF. The
#' vector \eqn{b} is called the b-vector parameterization of the FMP Model.
#' When \eqn{k=0}, the FMP IRF equals either the slope-threshold 
#' parameterization of the two-parameter item response model (if maxncat = 2) 
#' or Muraki's (1992) generalized partial credit model (if maxncat > 2).
#'
#' For \eqn{m(\theta)} to be a monotonic function, the FMP IRF can also be
#' expressed as a function of the vector
#'
#' \deqn{\gamma = (\xi, \omega, \alpha_1, \tau_1, \alpha_2, \tau_2,
#' \cdots \alpha_k,\tau_k)'.}{\gamma = (\xi, \omega, \alpha1, \tau1, \alpha2,
#' \tau2, ... \alpha_k, \tau_k)'.}
#'
#' The \eqn{\gamma} vector is called the Greek-letter parameterization of the
#' FMP model. See Falk & Cai (2016a), Feuerstahler (2016), or Liang & Browne 
#' (2015) for details about the relationship between the b-vector and 
#' Greek-letter parameterizations.
#'
#'
#' @examples
#'
#' set.seed(2345)
#' bmat <- sim_bmat(n_items = 5, k = 2, ncat = 4)$bmat
#'
#' theta <- rnorm(50)
#' dat <- sim_data(bmat = bmat, theta = theta, maxncat = 4)
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
#' ## random-effects estimation
#'
#' fmp0_random <- fmp(dat = dat, k = 0, em = TRUE)
#' 
#' ## random-effects estimation using mirt's estimation engine
#' 
#' fmp0_mirt <- fmp(dat = dat, k = 0, em = "mirt")
#'
#' @references
#'
#' Chalmers, R. P. (2012). mirt: A multidimensional item response theory
#' package for the R environment. \emph{Journal of Statistical Software},
#' \emph{48}, 1--29. \doi{10.18637/jss.v048.i06}
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
#' University of Minnesota, Minneapolis, MN.
#' \url{http://hdl.handle.net/11299/182267}
#'
#' Feuerstahler, L. M. (2019). Metric Transformations and the Filtered
#' Monotonic Polynomial Item Response Model. \emph{Psychometrika}, \emph{84},
#' 105--123. \doi{10.1007/s11336-018-9642-9}
#'
#' Liang, L. (2007). \emph{A semi-parametric approach to estimating item
#' response functions} (Unpublished dissertation). The Ohio
#' State University, Columbus, OH. Retrieved from https://etd.ohiolink.edu/
#'
#' Liang, L., & Browne, M. W. (2015). A quasi-parametric method for
#' fitting flexible item response functions. \emph{Journal of Educational
#' and Behavioral Statistics}, \emph{40}, 5--34. \doi{10.3102/1076998614556816}
#' 
#' Muraki, E. (1992). A generalized partial credit model: Application of an EM
#' algorithm. \emph{Applied Psychological Measurement}, \emph{16}, 159--176. 
#' \doi{10.1177/014662169201600206}
#'
#' @importFrom stats dnorm qnorm optimize
#' @export


fmp_1 <- function(dat, k, tsur, start_vals = NULL, method = "CG",
                  priors = list(xi = c("none", NaN, NaN),
                                omega = c("none", NaN, NaN),
                                alpha = c("none", NaN, NaN),
                                tau = c("none", NaN, NaN)), ...) {

  missing <- is.na(dat)
  if (any(missing)) {
    dat <- dat[!missing]
    tsur <- tsur[!missing]
    message(paste(sum(missing), "values removed due to missing data"))
  }

  ncat <- max(dat) + 1

  parmat <- data.frame("item" = rep(1, 2 * k + ncat))

  if (k == 0) parnames <- c(paste0("xi", 1:(ncat - 1)), "omega") else
    parnames <- c(paste0("xi", 1:(ncat - 1)), "omega",
                  sapply(1:k, function(x)
                    paste(c("alpha", "tau"), x, sep = "")))

  parmat$name <- parnames
  parmat$est <- TRUE
  parmat$value <- 0

  if (is.null(start_vals)) {
    parmat$value[grep("xi", parmat$name)] <-
      qnorm(cumsum(table(dat))[1:(ncat - 1)] / length(dat))
    parmat$value[grep("omega", parmat$name)] <- log(1)
    parmat$value[grepl("alpha", parmat$name) & parmat$est] <- .1
    parmat$value[grepl("tau", parmat$name) & parmat$est] <- log(.1)
  } else
    parmat$value <- start_vals

  if (is.null(priors$xi)) priors$xi <- c("none", NaN, NaN)
  if (is.null(priors$omega)) priors$omega <- c("none", NaN, NaN)
  if (is.null(priors$alpha)) priors$alpha <- c("none", NaN, NaN)
  if (is.null(priors$tau)) priors$tau <- c("none", NaN, NaN)

  parmat$prior_type <- "none"
  parmat$prior_2 <- parmat$prior_1 <- NaN
  parmat$prior_type[grep("xi", parmat$name)] <- priors$xi[1]
  parmat$prior_1[grep("xi", parmat$name)] <- as.numeric(priors$xi[2])
  parmat$prior_2[grep("xi", parmat$name)] <- as.numeric(priors$xi[3])
  parmat$prior_type[grep("omega", parmat$name)] <- priors$omega[1]
  parmat$prior_1[grep("omega", parmat$name)] <- as.numeric(priors$omega[2])
  parmat$prior_2[grep("omega", parmat$name)] <- as.numeric(priors$omega[3])
  parmat$prior_type[grep("alpha", parmat$name)] <- priors$alpha[1]
  parmat$prior_1[grep("alpha", parmat$name)] <- as.numeric(priors$alpha[2])
  parmat$prior_2[grep("alpha", parmat$name)] <- as.numeric(priors$alpha[3])
  parmat$prior_type[grep("tau", parmat$name)] <- priors$tau[1]
  parmat$prior_1[grep("tau", parmat$name)] <- as.numeric(priors$tau[2])
  parmat$prior_2[grep("tau", parmat$name)] <- as.numeric(priors$tau[3])

  ## estimate the model
  mod <- optim(par = parmat$value,
               fn = logl, gr = gr_logl,
               method = method,
               dat = dat, thetas = tsur,
               maxncat = ncat,
               parmat = parmat,
               control = ...)

  ## save the Greek parameters
  parmat$estimate <- mod$par

  ## find the standard errors
  info <- try(solve(hess_logl(parvec = parmat$estimate, maxncat = ncat,
                                   dat = dat, thetas = tsur, parmat = parmat)))

  parmat$error <- suppressWarnings(try(sqrt(diag(info))))


  xi <- parmat$estimate[grep("xi", parmat$name)]
  omega <- parmat$estimate[grep("omega", parmat$name)]
  alpha <- parmat$estimate[grep("alpha", parmat$name)]
  tau <- parmat$estimate[grep("tau", parmat$name)]


  if (k == 0) {
    bmat <- c(xi, exp(omega))
  } else{
    bmat <- greek2b(xi = xi, omega = omega, alpha = alpha, tau = tau)
  }

  out <- list(bmat = bmat, parmat = parmat,
              k = k, ncat = ncat, log_lik = mod$value, mod = mod,
              AIC = 2 * mod$value + 2 * sum(parmat$est),
              BIC = 2 * mod$value + sum(parmat$est) * log(length(dat)))

  out
}

#' @rdname fmp
#' @export

fmp <- function(dat, k, start_vals = NULL,
                em = TRUE, eps = 1e-04, n_quad = 49,
                method = "BFGS", max_em = 500,
                priors = list(xi = c("none", NaN, NaN),
                              omega = c("none", NaN, NaN),
                              alpha = c("none", NaN, NaN),
                              tau = c("none", NaN, NaN)), ...) {

  missing <- apply(dat, 1, function(d) all(is.na(d)))
  if (any(missing)) {
    dat <- subset(dat, !missing)
    message(paste(sum(missing), "rows removed due to missing data"))
  }

  ncat <- apply(dat, 2, max, na.rm = TRUE) + 1
  maxncat <- max(ncat)

  n_items <- ncol(dat)
  n_subj <- nrow(dat)

  maxk <- max(k) # highest order item

  if (length(k) != n_items) k <- rep_len(k, n_items)

  ### create a matrix of parameter information

  ## index the estimated parameters
  parmat <- data.frame("item" = c(rep(1:n_items, each = 2 * maxk + maxncat)))

  ## name the estimated parameters

  # item parameters (theta metric)
  if (maxk == 0) parnames <- c(paste0("xi", 1:(maxncat - 1)), "omega") else
    parnames <- c(paste0("xi", 1:(maxncat - 1)), "omega",
                  sapply(1:maxk, function(x)
                    paste(c("alpha", "tau"), x, sep = "")))

  parmat$name <- rep(parnames, n_items)

  # indicate whether parameters are estimated
  estmat <- sapply(1:n_items, function(i) c(rep(TRUE, ncat[i] - 1),
                                            rep(FALSE, maxncat - ncat[i]),
                                            rep(TRUE, 1 + 2 * k[i]),
                                            rep(FALSE, 2 * (maxk - k[i]))))
  parmat$est <- as.logical(estmat)

  # alpha/tau should equal zero/-Inf if not estimated
  # don't need to set alpha = -Inf since all alphas are estimated
  parmat$value <- 0
  parmat$value[grep("tau", parmat$name)] <- -Inf

  # add start values
  if (!is.null(start_vals)) parmat$value[parmat$est] <- start_vals
  if (is.null(start_vals) & em != "mirt") {
    for (i in 1:n_items) {
      parmat$value[grepl("xi", parmat$name) & parmat$item == i & parmat$est] <-
        qnorm(cumsum(table(dat[, i]))[1:(ncat[i] - 1)] / length(dat[, i]))
    }
    parmat$value[grep("omega", parmat$name)] <- log(1)
    parmat$value[grepl("alpha", parmat$name) & parmat$est] <- 0
    parmat$value[grepl("tau", parmat$name) & parmat$est] <- log(.01)
  }

  if (is.null(priors$xi)) priors$xi <- c("none", NaN, NaN)
  if (is.null(priors$omega)) priors$omega <- c("none", NaN, NaN)
  if (is.null(priors$alpha)) priors$alpha <- c("none", NaN, NaN)
  if (is.null(priors$tau)) priors$tau <- c("none", NaN, NaN)

  parmat$prior_type <- "none"
  parmat$prior_2 <- parmat$prior_1 <- NaN
  parmat$prior_type[grep("xi", parmat$name)] <- priors$xi[1]
  parmat$prior_1[grep("xi", parmat$name)] <- as.numeric(priors$xi[2])
  parmat$prior_2[grep("xi", parmat$name)] <- as.numeric(priors$xi[3])
  parmat$prior_type[grep("omega", parmat$name)] <- priors$omega[1]
  parmat$prior_1[grep("omega", parmat$name)] <- as.numeric(priors$omega[2])
  parmat$prior_2[grep("omega", parmat$name)] <- as.numeric(priors$omega[3])
  parmat$prior_type[grep("alpha", parmat$name)] <- priors$alpha[1]
  parmat$prior_1[grep("alpha", parmat$name)] <- as.numeric(priors$alpha[2])
  parmat$prior_2[grep("alpha", parmat$name)] <- as.numeric(priors$alpha[3])
  parmat$prior_type[grep("tau", parmat$name)] <- priors$tau[1]
  parmat$prior_1[grep("tau", parmat$name)] <- as.numeric(priors$tau[2])
  parmat$prior_2[grep("tau", parmat$name)] <- as.numeric(priors$tau[3])

  if (em == "mirt") {
    itemtype <- ifelse(k == 0, ifelse(ncat == 2, "2PL", "gpcm"), "monopoly")

    pars <- mirt::mirt(dat = as.data.frame(dat), model = 1,
                       itemtype = itemtype, monopoly.k = k,
                       pars = "values", ...)

    itemnames <- unique(pars$item)

    # if not given start values, record the ones used by mirt
    if (is.null(start_vals)) {
      for (i in 1:n_items) {
        subpars <- pars[pars$item == itemnames[i], ]
        vals <- parmat$value[parmat$item == i]
        if (k[i] == 0) {
          vals[1] <- subpars$value[2] # xi
          vals[2] <- log(subpars$value[1]) #omega
          parmat$value[parmat$item == i] <- vals
        } else { # if monopoly
          vals[1] <- subpars$value[2] # xi
          vals[2] <- subpars$value[1] #omega
          vals[3:nrow(subpars)] <- subpars$value[-c(1, 2)]
          parmat$value[parmat$item == i] <- vals
        }
      }
    } else{ # if given start values, add to pars
      for (i in 1:n_items) {
        vals <- parmat$value[parmat$est & parmat$item == i]
        if (k[i] == 0) {
          pars$value[pars$item == itemnames[i]] <-
            c(exp(vals[2]), vals[1], 0, 1) # a1, d, g = 0, u = 1
        } else { # if monopoly
          pars$value[pars$item == itemnames[i]] <-
            vals[c(2, 1, 3:length(vals))]
        }
      }
    }

    # priors
    pars$prior.type[grep("xi", pars$name)] <- priors$xi[1]
    pars$prior_1[grep("xi", pars$name)] <- as.numeric(priors$xi[2])
    pars$prior_2[grep("xi", pars$name)] <- as.numeric(priors$xi[3])
    pars$prior.type[grep("omega", pars$name)] <- priors$omega[1]
    pars$prior_1[grep("omega", pars$name)] <- as.numeric(priors$omega[2])
    pars$prior_2[grep("omega", pars$name)] <- as.numeric(priors$omega[3])
    pars$prior.type[grep("alpha", pars$name)] <- priors$alpha[1]
    pars$prior_1[grep("alpha", pars$name)] <- as.numeric(priors$alpha[2])
    pars$prior_2[grep("alpha", pars$name)] <- as.numeric(priors$alpha[3])
    pars$prior.type[grep("tau", pars$name)] <- priors$tau[1]
    pars$prior_1[grep("tau", pars$name)] <- as.numeric(priors$tau[2])
    pars$prior_2[grep("tau", pars$name)] <- as.numeric(priors$tau[3])

    mod <- mirt::mirt(dat = as.data.frame(dat), model = 1, itemtype = itemtype,
                      monopoly.k = k, pars = pars, TOL = eps,
                      quadpts = n_quad, ...)

    ## update bmat, parmat
    mirtests <- mirt::extract.mirt(mod, "parvec")
    mirtests <- tapply(mirtests, parmat$item[parmat$est], function(x) {
      x[c(1, 2:maxncat)] <- x[c(2:maxncat, 1)]
      if (length(x) == maxncat){
        x[maxncat] <- log(x[maxncat])
        if(maxncat > 2){
          for(i in 3:maxncat) x[i - 1] <- x[i - 1] - sum(x[1:(i - 2)])
        }
      } 
      x
    })

    parmat$estimate <- parmat$value
    parmat$estimate[parmat$est] <- unlist(mirtests)

    log_lik <- mirt::extract.mirt(mod, "logLik")
    aic <- mirt::extract.mirt(mod, "AIC")
    bic <- mirt::extract.mirt(mod, "BIC")

  } else{

    if (!em) {

      tsur <- get_surrogates(dat)

      mods <- lapply(1:n_items, function(it) {
        parmati <- parmat[parmat$item == it, ]
        parmati$item <- 1

        optim(par = parmati$value[parmati$est],
              fn = logl, gr = gr_logl,
              method = method,
              dat = dat[, it], thetas = tsur,
              maxncat = maxncat,
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

      ## find the standard errors
      info <- try(solve(hess_logl(parvec = parmat$estimate, dat = dat,
                                  maxncat = maxncat, thetas = tsur,
                                  parmat = parmat)))

      parmat$error <- suppressWarnings(try(sqrt(diag(info))))

      log_lik <- mod$value
      aic <- 2 * mod$value + 2 * sum(parmat$est)
      bic <- 2 * mod$value + sum(parmat$est) * log(n_subj)

  }

    if (em) {

      em_out <- em_alg(dat = dat,
                       eps = eps, n_quad = n_quad, maxncat = maxncat,
                       method = method, parmat = parmat, max_em = max_em, ...)

      mod <- em_out$mod
      mod$iter <- em_out$iter

      ## save the Greek parameters
      parmat$estimate <- parmat$value
      parmat$estimate[parmat$est] <- mod$par

      log_lik <- mod$value
      aic <- 2 * mod$value + 2 * sum(parmat$est)
      bic <- 2 * mod$value + sum(parmat$est) * log(n_subj)

    }

  }

  if (maxk > 0) {
    bmat <- t(sapply(1:n_items, function(i) {
      subparmat <- subset(parmat, parmat$item == i)
      greek2b(xi = subparmat$estimate[grep("xi", subparmat$name)],
              omega = subparmat$estimate[grep("omega", subparmat$name)],
              alpha = subparmat$estimate[grep("alpha", subparmat$name)],
              tau = subparmat$estimate[grep("tau", subparmat$name)])
    }))
  } else{
    bmat <- t(sapply(1:n_items, function(i) {
      subparmat <- subset(parmat, parmat$item == i)
      greek2b(xi = subparmat$estimate[grep("xi", subparmat$name)],
              omega = subparmat$estimate[grep("omega", subparmat$name)],
              alpha = NULL,
              tau = NULL)
    }))
  }

  if (maxncat > 2)
    colnames(bmat) <- c(paste0("b0_", 1:(maxncat - 1)),
                        paste0("b", 1:(2 * maxk + 1))) else
      colnames(bmat) <- paste0("b", 0:(ncol(bmat) - 1))

  out <- list(bmat = bmat, parmat = parmat,
              k = k, log_lik = log_lik, mod = mod,
              maxncat = maxncat, AIC = aic, BIC = bic)

  out

}
