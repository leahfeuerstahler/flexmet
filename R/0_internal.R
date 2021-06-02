partial_m <- function(parveci, theta, maxk, ncat = 2) {

  j <- 1:(2 * maxk + 1)

  # L&B Eq B10
  nu <- 1 / j * theta ^ j

  omegai <- parveci[ncat]

  if (maxk == 0) {

    out <- c(rep(1, ncat - 1), nu * exp(omegai))

  } else {

    # extract alpha and tau parameters
    alphai <- parveci[c(ncat - 1 + 2 * (1:maxk))]
    taui <- parveci[c(ncat + 2 * (1:maxk))]

    # find all T matrices
    t_list <- lapply(1:maxk, find_t, alpha = alphai, tau = taui)

    # compute partial w.r.t. omega - L&B Eq B11 or F&Ca Appendix 2
    m_omega <- t_list[[1]] * exp(omegai)

    if (maxk > 1) {

      for (j in 2:maxk) {
        m_omega <- t_list[[j]] %*% m_omega
      }
    }
    m_omega <- t(nu) %*% m_omega

    # compute partials w.r.t. alpha & tau - L&B Eq B11 or F&Ca Appendix 2
    m_alphatau <- matrix(nrow = 2, ncol = maxk)

    for (j in 1:maxk) {
      t_list_alpha <- t_list_tau <- t_list
      t_list_alpha[[j]] <- find_t_general(j = j, el1 = 0, el2 = -2,
                                         el3 = 2 * alphai[j])
      t_list_tau[[j]] <- find_t_general(j = j, el1 = 0, el2 = 0,
                                       el3 = exp(taui[j]))
      tmp_alpha <- tmp_tau <- exp(omegai)
      for (i in 1:length(t_list)) {
        tmp_alpha <- t_list_alpha[[i]] %*% tmp_alpha
        tmp_tau <- t_list_tau[[i]] %*% tmp_tau
      }

      m_alphatau[1, j] <- t(nu) %*% tmp_alpha
      m_alphatau[2, j] <- t(nu) %*% tmp_tau
    }

    # output vector of first derivatives
    out <- c(rep(1, ncat - 1), m_omega, as.numeric(m_alphatau))
  }
  out
}

partial2_m <- function(parveci, theta, maxk, ncat = 2) {

  j <- 1:(2 * maxk + 1)

  # L&B Eq B10
  nu <- 1 / j * theta ^ j

  omegai <- parveci[ncat]

  if (maxk == 0) {

    out <- matrix(0, nrow = 2, ncol = 2)
    out[2, 2] <- nu * exp(omegai)

  } else {

    out <- matrix(0, nrow = 2 * maxk + ncat, ncol = 2 * maxk + ncat)

    # extract alpha and tau parameters
    alphai <- parveci[c(ncat - 1 + 2 * (1:maxk))]
    taui <- parveci[c(ncat + 2 * (1:maxk))]

    # find T matrices
    t_list <- lapply(1:maxk, find_t, alpha = alphai, tau = taui)

    # w.r.t. omega
    m_omega <- t_list[[1]] * exp(omegai)

    if (maxk > 1) {
      for (j in 2:maxk) {
        m_omega <- t_list[[j]] %*% m_omega
      }
    }
    out[ncat, ncat] <- t(nu) %*% m_omega

    # get more T matrices
    for (j1 in 1:maxk) {
      t_list_alpha1 <- t_list_tau1 <- t_list_alpha_alpha <-
        t_list_tau_tau <- t_list

      t_list_alpha1[[j1]] <- find_t_general(j = j1, el1 = 0, el2 = -2,
                                           el3 = 2 * alphai[j1])

      t_list_tau1[[j1]] <- find_t_general(j = j1, el1 = 0, el2 = 0,
                                         el3 = exp(taui[j1]))

      t_list_alpha_alpha[[j1]] <- find_t_general(j = j1, el1 = 0,
                                                el2 = 0, el3 = 2)

      t_list_tau_tau[[j1]] <- find_t_general(j = j1, el1 = 0,
                                            el2 = 0, el3 = exp(taui[j1]))

      # w.r.t. alpha, omega
      tmp_alpha_omega <- exp(omegai)

      # w.r.t. tau, omega
      tmp_tau_omega <- exp(omegai)

      # w.r.t. alpha^2
      tmp_alpha_alpha <- exp(omegai)

      # w.r.t. tau^2
      tmp_tau_tau <- exp(omegai)

      # w.r.t. alpha, tau (same index) = 0

      for (i in 1:length(t_list_alpha1)) {
        tmp_alpha_omega <- t_list_alpha1[[i]] %*% tmp_alpha_omega
        tmp_tau_omega <- t_list_tau1[[i]] %*% tmp_tau_omega
        tmp_alpha_alpha <- t_list_alpha_alpha[[i]] %*% tmp_alpha_alpha
        tmp_tau_tau <- t_list_tau_tau[[i]] %*% tmp_tau_tau
      }

      out[ncat, ncat - 1 + 2 * j1] <- t(nu) %*% tmp_alpha_omega
      out[ncat, ncat + 2 * j1] <- t(nu) %*% tmp_tau_omega
      out[ncat - 1 + 2 * j1, ncat - 1 + 2 * j1] <- t(nu) %*% tmp_alpha_alpha
      out[ncat + 2 * j1, ncat + 2 * j1] <- t(nu) %*% tmp_tau_tau

      for (j2 in 1:maxk) {

        if (j1 < j2) {
          t_list_alpha2 <- t_list_alpha1
          t_list_alpha_tau1 <- t_list_alpha1
          t_list_alpha_tau2 <- t_list_tau1
          t_list_tau_tau <- t_list_tau1

          t_list_alpha2[[j2]] <- find_t_general(j = j2, el1 = 0, el2 = -2,
                                               el3 = 2 * alphai[j2])
          t_list_alpha_tau1[[j2]] <- find_t_general(j = j2, el1 = 0, el2 = 0,
                                                   el3 = exp(taui[[j2]]))
          t_list_alpha_tau2[[j2]] <- find_t_general(j = j2, el1 = 0, el2 = -2,
                                                   el3 = 2 * alphai[j2])
          t_list_tau_tau[[j2]] <- find_t_general(j = j2, el1 = 0, el2 = 0,
                                                el3 = exp(taui[[j2]]))


          # w.r.t. alpha, alpha (different indices)
          tmp_alpha_alpha <- exp(omegai)

          # w.r.t. alpha, tau (diferent indices)
          tmp_alpha_tau1 <- exp(omegai)
          tmp_alpha_tau2 <- exp(omegai)

          # w.r.t. tau, tau (different indices)
          tmp_tau_tau <- exp(omegai)

          for (i in 1:length(t_list_alpha1)) {
            tmp_alpha_alpha <- t_list_alpha2[[i]] %*% tmp_alpha_alpha
            tmp_alpha_tau1 <- t_list_alpha_tau1[[i]] %*% tmp_alpha_tau1
            tmp_alpha_tau2 <- t_list_alpha_tau2[[i]] %*% tmp_alpha_tau2
            tmp_tau_tau <- t_list_tau_tau[[i]] %*% tmp_tau_tau
          }

          out[ncat - 1 + 2 * j1, ncat - 1 + 2 * j2] <- t(nu) %*% tmp_alpha_alpha
          out[ncat - 1 + 2 * j1, ncat + 2 * j2] <- t(nu) %*% tmp_alpha_tau1
          out[ncat + 2 * j1, ncat - 1 + 2 * j2] <- t(nu) %*% tmp_alpha_tau2
          out[ncat + 2 * j1, ncat + 2 * j2] <- t(nu) %*% tmp_tau_tau
        }
      }
    }

  }
  t(out)[lower.tri(out, diag = TRUE)]
}

logl <- function(parvec, dat, thetas, maxncat = 2, parmat) {

  dat <- as.matrix(dat)

  n_items <- ncol(dat)

  parmat$value[parmat$est] <- parvec
  pars <- matrix(parmat$value, nrow = n_items, byrow = TRUE)

  maxk <- (ncol(pars) - maxncat) / 2

  xi <- pars[, 1:(maxncat - 1), drop = FALSE]
  omega <- pars[, maxncat]

  if (maxk > 0) {
    alpha <- as.matrix(pars[, c(maxncat - 1 + 2 * (1:maxk))]) ## 3, 5, 7, ...
    tau <- as.matrix(pars[, c(maxncat + 2 * (1:maxk))]) ## 4, 6, 8, ...

  if (n_items == 1) {
      alpha <- matrix(alpha, nrow = 1)
      tau <- matrix(tau, nrow = 1)
    }
  }

  # convert greek parameters to b parameters

  bmat <- matrix(nrow = n_items, ncol = 2 * maxk + maxncat)

  if (maxk == 0) {
    for (it in 1:n_items) {
      bmat[it, ] <- greek2b(xi = xi[it, ], omega = omega[it],
                           alpha = NULL, tau = NULL)
    }
  }

  if (maxk > 0) {
    for (it in 1:n_items) {
      bmat[it, ] <- greek2b(xi = xi[it, ], omega = omega[it],
                           alpha = alpha[it, ], tau = tau[it, ])
    }
  }

  # find probabilities
  probs <- irf_fmp(theta = thetas, bmat = bmat,
                   maxncat = maxncat, returncat = 0:(maxncat - 1))

  # avoid probs of 1 or 0
  probs[probs == 1] <- 1 - 1e-16
  probs[probs == 0] <- 1e-16

  logprobs <- sapply(1:dim(probs)[2], function(i)
    log(probs[, i, ][cbind(1:dim(probs)[1], dat[, i] + 1)]))

  # calculate log priors
  priormat <- subset(parmat, parmat$est & parmat$prior_type != "none")
  priorlogl <- ifelse(nrow(priormat) > 0,
                      log(get(paste0("d", priormat$prior_type))(
                        priormat$value, priormat$prior_1, priormat$prior_2)),
                      0)

  -sum(logprobs, na.rm = TRUE) - sum(priorlogl)
}

gr_logl <- function(parvec, dat, thetas, maxncat, parmat) {

  dat <- as.matrix(dat)

  n_items <- ncol(dat)
  ntheta <- nrow(dat)

  parmat$value[parmat$est] <- parvec
  
  pars <- matrix(parmat$value, nrow = n_items, byrow = TRUE)

  maxk <- (ncol(pars) - maxncat) / 2

  xi <- pars[, 1:(maxncat - 1), drop = FALSE]
  omega <- pars[, maxncat]

  if (maxk > 0) {
    alpha <- as.matrix(pars[, c(maxncat - 1 + 2 * (1:maxk))]) ## 3, 5, 7, ...
    tau <- as.matrix(pars[, c(maxncat + 2 * (1:maxk))]) ## 4, 6, 8, ...

    if (n_items == 1) {
      alpha <- matrix(alpha, nrow = 1)
      tau <- matrix(tau, nrow = 1)
    }
  }

  bmat <- matrix(nrow = n_items, ncol = 2 * maxk + maxncat)

  if (maxk == 0) {
    for (it in 1:n_items) {
      bmat[it, ] <- greek2b(xi = xi[it, ], omega = omega[it],
                           alpha = NULL, tau = NULL)
    }
  }

  if (maxk > 0) {
    for (it in 1:n_items) {
      bmat[it, ] <- greek2b(xi = xi[it, ], omega = omega[it],
                           alpha = alpha[it, ], tau = tau[it, ])
    }
  }

  # find probabilities
  probs <- irf_fmp(theta = thetas, bmat = bmat,
                   maxncat = maxncat, returncat = 0:(maxncat - 1))

  # avoid probs of 1 or 0
  probs[probs == 1] <- 1 - 1e-16
  probs[probs == 0] <- 1e-16

  out <- matrix(nrow = n_items, ncol = ncol(pars))

  for (it in 1:n_items) {
    # N x (2 * maxk + maxncat)
    partial_m_mat <- t(sapply(thetas, partial_m, ncat = maxncat,
                             parveci = pars[it, ], maxk = maxk))
    probs_it <- as.matrix(probs[, it, ]) # N x maxcat
    if (nrow(probs_it) != ntheta) probs_it <- t(probs_it)
    g_it <- 1 / probs_it[, 1, drop = FALSE] # denominator N x 1
    f_it <- as.matrix(apply(probs_it, 2, "*", g_it)) # numerator N x maxcat
    if (nrow(f_it) != ntheta) f_it <- t(f_it)

    f1_it <- array(NA, dim = c(dim(partial_m_mat), maxncat))
    f1_it[, , 1] <- 0
    for (c in 2:maxncat) {
      f1_it[, 1:(maxncat - 1), c] <- 0
      f1_it[, 1:(c - 1), c] <- partial_m_mat[, 1:(c - 1)] * f_it[, c]
      f1_it[, maxncat:(dim(partial_m_mat)[2]), c] <-
        partial_m_mat[, maxncat:(dim(partial_m_mat)[2])] * (c - 1) * f_it[, c]
    }
    f1_it[f1_it == Inf] <- 1e308
    f1_it[f1_it == -Inf] <- -1e308

    g1_it <- apply(f1_it, c(1, 2), sum, na.rm = TRUE) # N x (2 * maxk + maxncat)

    priorinfo <- subset(parmat, parmat$item == it)

    for (l in 1:dim(partial_m_mat)[2]) {

      # can currently only handle normal priors
      gr_prior <- ifelse(priorinfo$est[l] & priorinfo$prior_type[l] == "norm",
                         - (priorinfo$value[l] - priorinfo$prior_1[l]) /
                           priorinfo$prior_2[l]^2, 0)

      f1_it_l <- as.matrix(f1_it[, l, ])
      if (nrow(f1_it_l) != ntheta) f1_it_l <- t(f1_it_l)

      out[it, l] <- -sum(f1_it_l[cbind(1:ntheta, dat[, it] + 1)] /
                           f_it[cbind(1:ntheta, dat[, it] + 1)] -
                           g1_it[, l] / g_it) - gr_prior
    }
  }

  as.numeric(t(out))[parmat$est]

}

hess_logl <- function(parvec, dat, thetas, maxncat, parmat) {

  dat <- as.matrix(dat)

  n_items <- ncol(dat)

  parmat$value[parmat$est] <- parvec
  pars <- matrix(parmat$value, nrow = n_items, byrow = TRUE)

  maxk <- (ncol(pars) - 2) / 2

  # use analytic hessian if maxk = 0, else do cross-product of gradient
  if (maxncat == 2) {
    xi <- pars[, 1]
    omega <- pars[, 2]

    if (maxk > 0) {
      alpha <- as.matrix(pars[, c(1 + 2 * (1:maxk))]) ## 3, 5, 7, ...
      tau <- as.matrix(pars[, c(2 + 2 * (1:maxk))]) ## 4, 6, 8, ...

      if (n_items == 1) {
        alpha <- matrix(alpha, nrow = 1)
        tau <- matrix(tau, nrow = 1)
      }
    }

    # convert greek parameters to b parameters
    n_items <- 1

    bmat <- matrix(nrow = n_items, ncol = 2 * maxk + 2)

    if (maxk == 0) {
      for (it in 1:n_items) {
        bmat[it, ] <- greek2b(xi = xi[it], omega = omega[it],
                              alpha = NULL, tau = NULL)
      }
    }

    if (maxk > 0) {
      for (it in 1:n_items) {
        bmat[it, ] <- greek2b(xi = xi[it], omega = omega[it],
                              alpha = alpha[it, ], tau = tau[it, ])
      }
    }

    # find probabilities
    probs <- irf_fmp(theta = thetas, bmat = bmat)

    # avoid probs of 1 or 0
    probs[probs == 1] <- 1 - 1e-16
    probs[probs == 0] <- 1e-16

    # compute hessian--block matrix

    out <- matrix(0, nrow = n_items * (2 * maxk + 2),
                  ncol = n_items * (2 * maxk + 2))
    for (it in 1:n_items) {
      partial_m_mat <- sapply(thetas, partial_m,
                              parveci = pars[it, ], maxk = maxk)
      partial2_m_mat <- sapply(thetas, partial2_m,
                               parveci = pars[it, ], maxk = maxk)

      out_it <- idx <- matrix(0, nrow = 2 * maxk + 2, ncol = 2 * maxk + 2)
      idx[lower.tri(t(idx), diag = TRUE)] <- 1:((ncol(pars) + 1) *
                                                  ncol(pars) / 2)
      idx <- t(idx)

      for (l in 1:(2 * maxk + 2)) {
        for (t in 1:(2 * maxk + 2)) {
          if (l <= t) {
            out_it[l, t] <- out[t, l] <- sum(probs[, it] * (1 - probs[, it]) *
                                             partial_m_mat[l, ] *
                                             partial_m_mat[t, ] -
                                             (dat[, it] - probs[, it]) *
                                             partial2_m_mat[idx[l, t], ],
                                             na.rm = TRUE)
            out[((it - 1) * (2 * maxk + 2) + 1):(it * (2 * maxk + 2)),
                ((it - 1) * (2 * maxk + 2) + 1):(it * (2 * maxk + 2))] <- out_it
          }
        }
      }
    }
  } else{

    crossprods <- lapply(1:length(thetas), function(i) {
      gr <- gr_logl(parvec = parmat$value, dat = dat[i, ], thetas = thetas[i],
                    maxncat = maxncat, parmat = parmat)
      gr %*% t(gr)
    })

    out <- Reduce("+", crossprods)
  }

  out
}

em_alg <- function(dat, eps = 1e-04, maxncat = maxncat,
                   n_quad = 49, method, parmat, max_em, ...) {

  quad_nodes <- seq(-6, 6, length = n_quad)
  quad_wts <- dnorm(quad_nodes)
  quad_wts <- quad_wts / sum(quad_wts)

  n_subj <- nrow(dat)
  n_items <- ncol(dat)

  # start iterations
  iter <- 1
  estep_out <- e_step(dat = dat, parvec = parmat$value[parmat$est],
                     n_quad = n_quad, quad_nodes = quad_nodes,
                     quad_wts = quad_wts, n_subj = n_subj,
                     n_items = n_items, maxncat = maxncat, parmat = parmat)
  mstep_out <- m_step(estep_out,
                     parvec = estep_out$parvec,
                     quad_nodes = quad_nodes, quad_wts = quad_wts,
                     n_subj = n_subj, n_items = n_items, maxncat = maxncat,
                     method = method, parmat = parmat, ...)

  pars1 <- mstep_out$par
  maxchange <- NA

  cat("\riter: ", iter,
      " M step conv:", mstep_out$convergence,
      " maxchange = ", round(maxchange, 6),
      " -logl = ", mstep_out$value)
  maxchange <- 10

  # iterate
  while ((iter < 5 | abs(maxchange) > eps) & iter < max_em) {
    iter <- iter + 1
    pars <- pars1

    estep_out <- e_step(dat = dat, parvec = pars,
                       n_quad = n_quad, quad_nodes = quad_nodes,
                       quad_wts = quad_wts, n_subj = n_subj,
                       n_items = n_items, maxncat = maxncat, parmat = parmat)

    mstep_out <- m_step(estep_out,
                       parvec = estep_out$parvec,
                       quad_nodes = quad_nodes, quad_wts = quad_wts,
                       n_subj = n_subj, n_items = n_items, method = method,
                       maxncat = maxncat, parmat = parmat, ...)

    pars1 <- mstep_out$par
    maxchange <- max(abs(pars1 - pars))

    cat("\riter: ", iter,
        " M step conv:", mstep_out$convergence,
        " maxchange = ", round(maxchange, 6),
        " -logl = ", mstep_out$value)

  }

  cat("\n")

  mstep_out$n_bar <- estep_out$n_bar
  mstep_out$r_bar <- estep_out$r_bar
  mstep_out$quad_nodes <- quad_nodes
  mstep_out$quad_wts <- quad_wts

  list(mod = mstep_out, iter = iter, maxchange = maxchange)
}

e_step <- function(dat, parvec, n_quad, quad_nodes, quad_wts,
                   n_subj, n_items, maxncat, parmat) {

  dat <- as.matrix(dat)


  parmat$value[parmat$est] <- parvec
  pars <- matrix(parmat$value, nrow = n_items, byrow = TRUE)

  maxk <- (ncol(pars) - maxncat) / 2

  xi <- pars[, 1:(maxncat - 1), drop = FALSE]
  omega <- pars[, maxncat]


  if (maxk > 0) {
    alpha <- as.matrix(pars[, c(maxncat - 1 + 2 * (1:maxk))]) # 3, 5, 7, ...
    tau <- as.matrix(pars[, c(maxncat + 2 * (1:maxk))]) # 4, 6, 8, ...

    if (n_items == 1) {
      alpha <- matrix(alpha, nrow = 1)
      tau <- matrix(tau, nrow = 1)
    }
  }

  # convert greek parameters to b parameters
  bmat <- matrix(nrow = n_items, ncol = 2 * maxk + maxncat)
  if (maxk == 0) {
    for (it in 1:n_items) {
      bmat[it, ] <- greek2b(xi = xi[it, ], omega = omega[it],
                           alpha = NULL, tau = NULL)
    }
  }

  if (maxk > 0) {
    for (it in 1:n_items) {
      bmat[it, ] <- greek2b(xi = xi[it, ], omega = omega[it],
                           alpha = alpha[it, ], tau = tau[it, ])
    }
  }

  probs <- irf_fmp(theta = quad_nodes, bmat = bmat, maxncat = maxncat,
                   returncat = 0:(maxncat - 1)) # n_quad x n_items x maxncat

# lxr is n_subj x n_quad
  lxr <- t(apply(dat, 1, function(x)
    apply(probs, 1, function(y)
      prod(y[cbind(1:n_items, x + 1)], na.rm = TRUE))))

  # vector of denominators of length n_subj
  denoms <- apply(lxr, 1, function(x) t(x) %*% quad_wts)

  # expected item score:
  r_bar <- array(dim = c(n_quad, n_items, maxncat))
  for (c in 0:(maxncat - 1)) {
    r_bar[, , c + 1] <- sapply(1:n_items, function(j) # n_quad x n_items
      sapply(1:n_quad, function(r)
        sum((dat[, j] == c) * lxr[, r] * quad_wts[r] / denoms, na.rm = TRUE)))
  }

  list(r_bar = r_bar, parvec = parvec)
}

m_step <- function(estep_out, parvec, quad_nodes, quad_wts,
                   n_subj, n_items, maxncat, method, parmat, ...) {

  optim(par = parvec,
        fn = logl_em, gr = gr_logl_em,
        method = method,
        r_bar = estep_out$r_bar,
        quad_nodes = quad_nodes, quad_wts = quad_wts,
        n_subj = n_subj, n_items = n_items,
        maxncat = maxncat, parmat = parmat, ...)
}

logl_em <- function(parvec,
                    r_bar, quad_nodes, quad_wts,
                    n_subj, n_items, maxncat, parmat) {

  parmat$value[parmat$est] <- parvec
  pars <- matrix(parmat$value, nrow = n_items, byrow = TRUE)

  maxk <- (ncol(pars) - maxncat) / 2

  xi <- pars[, 1:(maxncat - 1), drop = FALSE]
  omega <- pars[, maxncat]

  if (maxk > 0) {
    alpha <- as.matrix(pars[, c(maxncat - 1 + 2 * (1:maxk))]) # 3, 5, 7, ...
    tau <- as.matrix(pars[, c(maxncat + 2 * (1:maxk))]) # 4, 6, 8, ...

    if (n_items == 1) {
      alpha <- matrix(alpha, nrow = 1)
      tau <- matrix(tau, nrow = 1)
    }
  }

  bmat <- matrix(nrow = n_items, ncol = 2 * maxk + maxncat)

  if (maxk == 0) {
    for (it in 1:n_items) {
      bmat[it, ] <- greek2b(xi = xi[it, ], omega = omega[it],
                           alpha = NULL, tau = NULL)
    }
  }

  if (maxk > 0) {
    for (it in 1:n_items) {
      bmat[it, ] <- greek2b(xi = xi[it, ], omega = omega[it],
                           alpha = alpha[it, ], tau = tau[it, ])
    }
  }
  probs <- irf_fmp(theta = quad_nodes, bmat = bmat,
                   maxncat = maxncat, returncat = 0:(maxncat - 1))

  # avoid probs of 1 or 0
  probs[probs == 1] <- 1 - 2e-16
  probs[probs == 0] <- 2e-16

  priormat <- subset(parmat, parmat$est & parmat$prior_type != "none")
  priorlogl <- ifelse(nrow(priormat) > 0,
                      log(get(paste0("d", priormat$prior_type))(
                        priormat$value, priormat$prior_1, priormat$prior_2)),
                      0)

  -sum(r_bar * log(probs)) - sum(priorlogl)
}

gr_logl_em <- function(parvec,
                       r_bar, quad_nodes, quad_wts,
                       maxncat = maxncat,
                       n_subj, n_items, parmat) {

  parmat$value[parmat$est] <- parvec
  pars <- matrix(parmat$value, nrow = n_items, byrow = TRUE)

  maxk <- (ncol(pars) - maxncat) / 2

  xi <- pars[, 1:(maxncat - 1), drop = FALSE]
  omega <- pars[, maxncat]

  if (maxk > 0) {
    alpha <- as.matrix(pars[, c(maxncat - 1 + 2 * (1:maxk))]) # 3, 5, 7, ...
    tau <- as.matrix(pars[, c(maxncat + 2 * (1:maxk))]) # 4, 6, 8, ...

    if (n_items == 1) {
      alpha <- matrix(alpha, nrow = 1)
      tau <- matrix(tau, nrow = 1)
    }
  }

  bmat <- matrix(nrow = n_items, ncol = 2 * maxk + maxncat)

  if (maxk == 0) {
    for (it in 1:n_items) {
      bmat[it, ] <- greek2b(xi = xi[it, ], omega = omega[it],
                           alpha = NULL, tau = NULL)
    }
  }

  if (maxk > 0) {
    for (it in 1:n_items) {
      bmat[it, ] <- greek2b(xi = xi[it, ], omega = omega[it],
                           alpha = alpha[it, ], tau = tau[it, ])
    }
  }

  probs <- irf_fmp(theta = quad_nodes, bmat = bmat,
                   maxncat = maxncat, returncat = 0:(maxncat - 1))

  ## avoid probs of 1 or 0
  probs[probs == 1] <- 1 - 1e-16
  probs[probs == 0] <- 1e-16

  out <- matrix(nrow = n_items, ncol = ncol(pars))

  for (it in 1:n_items) {
    partial_m_mat <- t(sapply(quad_nodes, partial_m, ncat = maxncat,
                              parveci = pars[it, ], maxk = maxk))
    probs_it <- as.matrix(probs[, it, ]) # Q x maxcat
    if (nrow(probs_it) != length(quad_nodes)) probs_it <- t(probs_it)
    g_it <- 1 / probs_it[, 1, drop = FALSE] # denominator Q x 1
    f_it <- as.matrix(apply(probs_it, 2, "*", g_it)) # numerator Q x maxcat
    if (nrow(f_it) != length(quad_nodes)) f_it <- t(f_it)

    f1_it <- array(NA, dim = c(dim(partial_m_mat), maxncat))
    f1_it[, , 1] <- 0
    for (c in 2:maxncat) {
      f1_it[, 1:(maxncat - 1), c] <- 0
      f1_it[, 1:(c - 1), c] <- partial_m_mat[, 1:(c - 1)] * f_it[, c]
      f1_it[, maxncat:(dim(partial_m_mat)[2]), c] <-
        partial_m_mat[, maxncat:(dim(partial_m_mat)[2])] * (c - 1) * f_it[, c]
    }
    f1_it[f1_it == Inf] <- 1e308
    f1_it[f1_it == -Inf] <- -1e308

    g1_it <- apply(f1_it, c(1, 2), sum, na.rm = TRUE) # Q x (2 * maxk + maxncat)

    priorinfo <- subset(parmat, parmat$item == it)

    for (l in 1:dim(partial_m_mat)[2]) {
      # can currently only handle normal priors
      gr_prior <- ifelse(priorinfo$est[l] & priorinfo$prior_type[l] == "norm",
                         - (priorinfo$value[l] - priorinfo$prior_1[l]) /
                           priorinfo$prior_2[l]^2, 0)

      f1_it_l <- as.matrix(f1_it[, l, ])
      if (nrow(f1_it_l) != length(quad_nodes)) f1_it_l <- t(f1_it_l)

      out[it, l] <- -sum(r_bar[, it, ] *
                           (f1_it_l / f_it - (g1_it[, l] / g_it) %*%
                              t(rep(1, maxncat))), na.rm = TRUE) - gr_prior
    }
  }

  estmat <- matrix(parmat$est, nrow = n_items, byrow = TRUE)

  as.numeric(t(out)[t(estmat)])
}

find_t <- function(j, alpha, tau) {
  out <- matrix(0, nrow = 2 * j + 1, ncol = 2 * j - 1)
  diag(out) <- 1

  if (j > 1) {
    diag(out[-1, ]) <- -2 * alpha[j]
    diag(out[-c(1, 2), ]) <- alpha[j] ^ 2 + exp(tau[j])
  }

  if (j == 1) {
    out[2, ] <- -2 * alpha[j]
    out[3, ] <- alpha[j] ^ 2 + exp(tau[j])
  }
  out
}

find_t_general <- function(j, el1, el2, el3) {
  out <- matrix(0, nrow = 2 * j + 1, ncol = 2 * j - 1)
  diag(out) <- el1

  if (j > 1) {
    diag(out[-1, ]) <- el2
    diag(out[-c(1, 2), ]) <- el3
  }

  if (j == 1) {
    out[2, ] <- el2
    out[3, ] <- el3
  }

  out
}

find_w <- function(tvec, k) {

  k_theta <- (length(tvec) - 2) / 2
  k_star <- 2 * k * k_theta + k + k_theta

  w <- matrix(0, nrow = 2 * k_star + 2, ncol = 2 * k + 2)

  wlist <- list()
  wlist[[1]] <- 1
  for (s in 2:(2 * k + 2)) {
    vs <- find_v(wlist[[s - 1]], k_theta)
    wlist[[s]] <- vs %*% tvec
  }

  for (s in 1:length(wlist)) {
    w[1:length(wlist[[s]]), s] <- wlist[[s]]
  }
  w
}

find_v <- function(wvec, k_theta) {
  v <- matrix(0, nrow = length(wvec) + 2 * k_theta + 1,
                 ncol = 2 * k_theta + 2)

  for (i in seq_len(ncol(v))) v[i:(i + length(wvec) - 1), i] <- wvec

  v
}
