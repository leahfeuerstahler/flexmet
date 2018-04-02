partial_m <- function(parveci, theta, maxk){

  j <- 1:(2 * maxk + 1)

  # L&B Eq B10
  nu <- 1 / j * theta ^ j

  omegai <- parveci[2]

  if (maxk == 0){

    out <- c(1, nu * exp(omegai))

  } else {

    # extract alpha and tau parameters
    alphai <- parveci[c(1 + 2 * (1:maxk))]
    taui <- parveci[c(2 + 2 * (1:maxk))]

    # find all T matrices
    Tlist <- lapply(1:maxk, find_t, alpha = alphai, tau = taui)

    # compute partial w.r.t. omega - L&B Eq B11 or F&Ca Appendix 2
    m.omega <- Tlist[[1]] * exp(omegai)

    if (maxk > 1){

      for (j in 2:maxk){
        m.omega <- Tlist[[j]] %*% m.omega
      }
    }
    m.omega <- t(nu) %*% m.omega

    # compute partials w.r.t. alpha & tau - L&B Eq B11 or F&Ca Appendix 2
    m.alphatau <- matrix(nrow = 2, ncol = maxk)

    for (j in 1:maxk){
      Tlist.alpha <- Tlist.tau <- Tlist
      Tlist.alpha[[j]] <- find_t_general(j = j, el1 = 0, el2 = -2,
                                         el3 = 2 * alphai[j])
      Tlist.tau[[j]] <- find_t_general(j = j, el1 = 0, el2 = 0,
                                       el3 = exp(taui[j]))
      tmp.alpha <- tmp.tau <- exp(omegai)
      for (i in 1:length(Tlist)){
        tmp.alpha <- Tlist.alpha[[i]] %*% tmp.alpha
        tmp.tau <- Tlist.tau[[i]] %*% tmp.tau
      }

      m.alphatau[1, j] <- t(nu) %*% tmp.alpha
      m.alphatau[2, j] <- t(nu) %*% tmp.tau
    }

    # output vector of first derivatives
    out <- c(1, m.omega, as.numeric(m.alphatau))
  }
  out
}

partial2_m <- function(parveci, theta, maxk){

  j <- 1:(2 * maxk + 1)

  # L&B Eq B10
  nu <- 1 / j * theta ^ j

  omegai <- parveci[2]

  if (maxk == 0){

    out <- matrix(0, nrow = 2, ncol = 2)
    out[2, 2] <- nu * exp(omegai)

  } else {

    out <- matrix(0, nrow = 2 * maxk + 2, ncol = 2 * maxk + 2)

    # extract alpha and tau parameters
    alphai <- parveci[c(1 + 2 * (1:maxk))]
    taui <- parveci[c(2 + 2 * (1:maxk))]

    # find T matrices
    Tlist <- lapply(1:maxk, find_t, alpha = alphai, tau = taui)

    # w.r.t. omega
    m.omega <- Tlist[[1]] * exp(omegai)

    if (maxk > 1){
      for (j in 2:maxk){
        m.omega <- Tlist[[j]] %*% m.omega
      }
    }
    out[2, 2] <- t(nu) %*% m.omega

    # get more T matrices
    for (j1 in 1:maxk){
      Tlist_alpha1 <- Tlist_tau1 <- Tlist_alpha_alpha <- Tlist_tau_tau <- Tlist

      Tlist_alpha1[[j1]] <- find_t_general(j = j1, el1 = 0, el2 = -2,
                                           el3 = 2 * alphai[j1])

      Tlist_tau1[[j1]] <- find_t_general(j = j1, el1 = 0, el2 = 0,
                                         el3 = exp(taui[j1]))

      Tlist_alpha_alpha[[j1]] <- find_t_general(j = j1, el1 = 0,
                                                el2 = 0, el3 = 2)

      Tlist_tau_tau[[j1]] <- find_t_general(j = j1, el1 = 0,
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

      for (i in 1:length(Tlist_alpha1)){
        tmp_alpha_omega <- Tlist_alpha1[[i]] %*% tmp_alpha_omega
        tmp_tau_omega <- Tlist_tau1[[i]] %*% tmp_tau_omega
        tmp_alpha_alpha <- Tlist_alpha_alpha[[i]] %*% tmp_alpha_alpha
        tmp_tau_tau <- Tlist_tau_tau[[i]] %*% tmp_tau_tau
      }

      out[2, 1 + 2 * j1] <- t(nu) %*% tmp_alpha_omega
      out[2, 2 + 2 * j1] <- t(nu) %*% tmp_tau_omega
      out[1 + 2 * j1, 1 + 2 * j1] <- t(nu) %*% tmp_alpha_alpha
      out[2 + 2 * j1, 2 + 2 * j1] <- t(nu) %*% tmp_tau_tau

      for (j2 in 1:maxk){

        if (j1 < j2){
          Tlist_alpha2 <- Tlist_alpha1
          Tlist_alpha_tau1 <- Tlist_alpha1
          Tlist_alpha_tau2 <- Tlist_tau1
          Tlist_tau_tau <- Tlist_tau1

          Tlist_alpha2[[j2]] <- find_t_general(j = j2, el1 = 0, el2 = -2,
                                               el3 = 2 * alphai[j2])
          Tlist_alpha_tau1[[j2]] <- find_t_general(j = j2, el1 = 0, el2 = 0,
                                                   el3 = exp(taui[[j2]]))
          Tlist_alpha_tau2[[j2]] <- find_t_general(j = j2, el1 = 0, el2 = -2,
                                                   el3 = 2 * alphai[j2])
          Tlist_tau_tau[[j2]] <- find_t_general(j = j2, el1 = 0, el2 = 0,
                                                el3 = exp(taui[[j2]]))


          # w.r.t. alpha, alpha (different indices)
          tmp_alpha_alpha <- exp(omegai)

          # w.r.t. alpha, tau (diferent indices)
          tmp_alpha_tau1 <- exp(omegai)
          tmp_alpha_tau2 <- exp(omegai)

          # w.r.t. tau, tau (different indices)
          tmp_tau_tau <- exp(omegai)

          for (i in 1:length(Tlist_alpha1)){
            tmp_alpha_alpha <- Tlist_alpha2[[i]] %*% tmp_alpha_alpha
            tmp_alpha_tau1 <- Tlist_alpha_tau1[[i]] %*% tmp_alpha_tau1
            tmp_alpha_tau2 <- Tlist_alpha_tau2[[i]] %*% tmp_alpha_tau2
            tmp_tau_tau <- Tlist_tau_tau[[i]] %*% tmp_tau_tau
          }

          out[1 + 2 * j1, 1 + 2 * j2] <- t(nu) %*% tmp_alpha_alpha
          out[1 + 2 * j1, 2 + 2 * j2] <- t(nu) %*% tmp_alpha_tau1
          out[2 + 2 * j1, 1 + 2 * j2] <- t(nu) %*% tmp_alpha_tau2
          out[2 + 2 * j1, 2 + 2 * j2] <- t(nu) %*% tmp_tau_tau
        }
      }
    }

  }
  t(out)[lower.tri(out, diag = TRUE)]
  #out
}

logl <- function(parvec, dat, thetas, parmat){

  dat <- as.matrix(dat)

  n_items <- ncol(dat)

  parmat$value[parmat$est] <- parvec
  pars <- matrix(parmat$value, nrow = n_items, byrow = TRUE)

  maxk <- (ncol(pars) - 2) / 2

  xi <- pars[, 1]
  omega <- pars[, 2]

  if (maxk > 0){
    alpha <- as.matrix(pars[, c(1 + 2 * (1:maxk))]) ## 3, 5, 7, ...
    tau <- as.matrix(pars[, c(2 + 2 * (1:maxk))]) ## 4, 6, 8, ...

    if (n_items == 1){
      alpha <- matrix(alpha, nrow = 1)
      tau <- matrix(tau, nrow = 1)
    }
  }

  # convert greek parameters to b parameters

  bmat <- matrix(nrow = n_items, ncol = 2 * maxk + 2)

  if (maxk == 0){
    for (it in 1:n_items){
      bmat[it, ] <- greek2b(xi = xi[it], omega = omega[it],
                           alpha = NULL, tau = NULL)
    }
  }

  if (maxk > 0){
    for (it in 1:n_items){
      bmat[it, ] <- greek2b(xi = xi[it], omega = omega[it],
                           alpha = alpha[it, ], tau = tau[it, ])
    }
  }

  # find probabilities
  probs <- irf_fmp(theta = thetas, bmat = bmat)

  # avoid probs of 1 or 0
  probs[probs == 1] <- 1 - 1e-16
  probs[probs == 0] <- 1e-16

  -sum(dat * log(probs) + (1 - dat) * log(1 - probs))
}

gr_logl <- function(parvec, dat, thetas, parmat){

  dat <- as.matrix(dat)

  n_items <- ncol(dat)

  parmat$value[parmat$est] <- parvec
  pars <- matrix(parmat$value, nrow = n_items, byrow = TRUE)

  maxk <- (ncol(pars) - 2) / 2

  xi <- pars[, 1]
  omega <- pars[, 2]

  if (maxk > 0){
    alpha <- as.matrix(pars[, c(1 + 2 * (1:maxk))]) ## 3, 5, 7, ...
    tau <- as.matrix(pars[, c(2 + 2 * (1:maxk))]) ## 4, 6, 8, ...

    if (n_items == 1){
      alpha <- matrix(alpha, nrow = 1)
      tau <- matrix(tau, nrow = 1)
    }
  }

  bmat <- matrix(nrow = n_items, ncol = 2 * maxk + 2)

  if (maxk == 0){
    for (it in 1:n_items){
      bmat[it, ] <- greek2b(xi = xi[it], omega = omega[it],
                           alpha = NULL, tau = NULL)
    }
  }

  if (maxk > 0){
    for (it in 1:n_items){
      bmat[it, ] <- greek2b(xi = xi[it], omega = omega[it],
                           alpha = alpha[it, ], tau = tau[it, ])
    }
  }

  probs <- irf_fmp(theta = thetas, bmat = bmat)

  # avoid probs of 1 or 0
  probs[probs == 1] <- 1 - 1e-16
  probs[probs == 0] <- 1e-16

  out <- matrix(nrow = n_items, ncol = ncol(pars))

  for (it in 1:n_items){
    partial_m_mat <- sapply(thetas, partial_m,
                            parveci = pars[it, ], maxk = maxk)

    for (l in 1:ncol(pars)){
      out[it, l] <- -sum( (dat[, it] - probs[, it]) * partial_m_mat[l, ])
    }
  }

  as.numeric(t(out))[parmat$est]

}

hess_logl <- function(parvec, dat, thetas, parmat){

  dat <- as.matrix(dat)

  n_items <- ncol(dat)

  parmat$value[parmat$est] <- parvec
  pars <- matrix(parmat$value, nrow = n_items, byrow = TRUE)

  maxk <- (ncol(pars) - 2) / 2

  xi <- pars[, 1]
  omega <- pars[, 2]

  if (maxk > 0){
    alpha <- as.matrix(pars[, c(1 + 2 * (1:maxk))]) ## 3, 5, 7, ...
    tau <- as.matrix(pars[, c(2 + 2 * (1:maxk))]) ## 4, 6, 8, ...

    if (n_items == 1){
      alpha <- matrix(alpha, nrow = 1)
      tau <- matrix(tau, nrow = 1)
    }
  }

  # convert greek parameters to b parameters
  n_items <- 1

  bmat <- matrix(nrow = n_items, ncol = 2 * maxk + 2)

  if (maxk == 0){
    for (it in 1:n_items){
      bmat[it, ] <- greek2b(xi = xi[it], omega = omega[it],
                            alpha = NULL, tau = NULL)
    }
  }

  if (maxk > 0){
    for (it in 1:n_items){
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

  OUT <- matrix(0, nrow = n_items * (2 * maxk + 2),
                ncol = n_items * (2 * maxk + 2))
  for (it in 1:n_items){
    partial_m_mat <- sapply(thetas, partial_m,
                            parveci = pars[it, ], maxk = maxk)
    partial2_m_mat <- sapply(thetas, partial2_m,
                             parveci = pars[it, ], maxk = maxk)

    out <- idx <- matrix(0, nrow = 2 * maxk + 2, ncol = 2 * maxk + 2)
    idx[lower.tri(t(idx), diag = TRUE)] <- 1:( (ncol(pars) + 1) *
                                                 ncol(pars) / 2)
    idx <- t(idx)

    for (l in 1:(2 * maxk + 2)){
      for (t in 1:(2 * maxk + 2)){
        if (l <= t){
          out[l, t] <- out[t, l] <- sum(probs[, it] * (1 - probs[, it]) *
                                        partial_m_mat[l, ] *
                                        partial_m_mat[t, ] -
                                        (dat[, it] - probs[, it]) *
                                        partial2_m_mat[idx[l, t], ])
          OUT[ ( (it - 1) * (2 * maxk + 2) + 1):(it * (2 * maxk + 2)),
              ( (it - 1) * (2 * maxk + 2) + 1):(it * (2 * maxk + 2))] <- out
        }
      }
    }
  }

  OUT
}

em_alg <- function(dat, eps = 1e-04,
                  n_quad = 50, method, parmat, max_em, ...){

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
                     n_items = n_items, parmat = parmat)
  mstep_out <- m_step(estep_out,
                     parvec = estep_out$parvec,
                     quad_nodes = quad_nodes, quad_wts = quad_wts,
                     n_subj = n_subj, n_items = n_items, method = method,
                     parmat = parmat, ...)

  pars1 <- mstep_out$par
  maxchange <- NA

  cat("\riter: ", iter,
      " M step conv:", mstep_out$convergence,
      " maxchange = ", round(maxchange, 6),
      " -logl = ", mstep_out$value)
  maxchange <- 10

  # iterate
  while (abs(maxchange) > eps & iter < max_em){
    iter <- iter + 1
    pars <- pars1

    estep_out <- e_step(dat = dat, parvec = pars,
                       n_quad = n_quad, quad_nodes = quad_nodes,
                       quad_wts = quad_wts, n_subj = n_subj,
                       n_items = n_items, parmat = parmat)

    mstep_out <- m_step(estep_out,
                       parvec = estep_out$parvec,
                       quad_nodes = quad_nodes, quad_wts = quad_wts,
                       n_subj = n_subj, n_items = n_items, method = method,
                       parmat = parmat, ...)

    pars1 <- mstep_out$par
    maxchange <- max(abs(pars1 - pars))
    # cat(paste(rep("\b", 75), collapse = ""),
    #     "iter: ", iter,
    #     " M step conv:", mstep_out$convergence,
    #     " maxchange = ", round(maxchange, 6),
    #     " -logl = ", mstep_out$value)

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
                   n_subj, n_items, parmat){

  dat <- as.matrix(dat)


  parmat$value[parmat$est] <- parvec
  pars <- matrix(parmat$value, nrow = n_items, byrow = TRUE)

  maxk <- (ncol(pars) - 2) / 2

  xi <- pars[, 1]
  omega <- pars[, 2]


  if (maxk > 0){
    alpha <- as.matrix(pars[, c(1 + 2 * (1:maxk))]) # 3, 5, 7, ...
    tau <- as.matrix(pars[, c(2 + 2 * (1:maxk))]) # 4, 6, 8, ...

    if (n_items == 1){
      alpha <- matrix(alpha, nrow = 1)
      tau <- matrix(tau, nrow = 1)
    }
  }

  # convert greek parameters to b parameters
  bmat <- matrix(nrow = n_items, ncol = 2 * maxk + 2)
  if (maxk == 0){
    for (it in 1:n_items){
      bmat[it, ] <- greek2b(xi = xi[it], omega = omega[it],
                           alpha = NULL, tau = NULL)
    }
  }

  if (maxk > 0){
    for (it in 1:n_items){
      bmat[it, ] <- greek2b(xi = xi[it], omega = omega[it],
                           alpha = alpha[it, ], tau = tau[it, ])
    }
  }

  probs <- irf_fmp(theta = quad_nodes, bmat = bmat) # n_quad x n_items

# lxr is n_subj x n_quad
  lxr <- t(apply(dat, 1, function(x) # dat is n_subj x n_items
    apply(probs, 1, function(y) # probs is n_quad x n_items
      prod(y ^ x * (1 - y) ^ (1 - x)))))

  # vector of denominators of length n_subj
  denoms <- apply(lxr, 1, function(x) t(x) %*% quad_wts)

  # expected item score:
  r_bar <- sapply(1:n_items, function(j) # n_quad x n_items
    sapply(1:n_quad, function(r)
      sum(dat[, j] * lxr[, r] * quad_wts[r] / denoms)))

  # expected number of persons at each quadpt
  n_bar <- sapply(1:n_quad, function(r) sum(lxr[, r] * quad_wts[r] / denoms))

  list(r_bar = r_bar, n_bar = n_bar,
       parvec = parvec)
}

m_step <- function(estep_out, parvec, quad_nodes, quad_wts,
                   n_subj, n_items, method, parmat, ...){

  optim(par = parvec,
        fn = logl_em, gr = gr_logl_em,
        method = method,
        n_bar = estep_out$n_bar, r_bar = estep_out$r_bar,
        quad_nodes = quad_nodes, quad_wts = quad_wts,
        n_subj = n_subj, n_items = n_items, parmat = parmat, ...)
}

logl_em <- function(parvec, n_bar, r_bar, quad_nodes, quad_wts,
                    n_subj, n_items, parmat){

  n_bar_mat <- n_bar %*% t(rep(1, n_items))

  parmat$value[parmat$est] <- parvec
  pars <- matrix(parmat$value, nrow = n_items, byrow = TRUE)

  maxk <- (ncol(pars) - 2) / 2

  xi <- pars[, 1]
  omega <- pars[, 2]

  if (maxk > 0){
    alpha <- as.matrix(pars[, c(1 + 2 * (1:maxk))]) # 3, 5, 7, ...
    tau <- as.matrix(pars[, c(2 + 2 * (1:maxk))]) # 4, 6, 8, ...

    if (n_items == 1){
      alpha <- matrix(alpha, nrow = 1)
      tau <- matrix(tau, nrow = 1)
    }
  }

  bmat <- matrix(nrow = n_items, ncol = 2 * maxk + 2)

  if (maxk == 0){
    for (it in 1:n_items){
      bmat[it, ] <- greek2b(xi = xi[it], omega = omega[it],
                           alpha = NULL, tau = NULL)
    }
  }

  if (maxk > 0){
    for (it in 1:n_items){
      bmat[it, ] <- greek2b(xi = xi[it], omega = omega[it],
                           alpha = alpha[it, ], tau = tau[it, ])
    }
  }
  probs <- irf_fmp(theta = quad_nodes, bmat = bmat)

  # avoid probs of 1 or 0
  probs[probs == 1] <- 1 - 2e-16
  probs[probs == 0] <- 2e-16

  -sum(r_bar * log(probs) + (n_bar_mat - r_bar) * log(1 - probs))
}

gr_logl_em <- function(parvec, n_bar, r_bar, quad_nodes, quad_wts,
                       n_subj, n_items, parmat){

  r_bar <- matrix(r_bar, ncol = n_items)

  parmat$value[parmat$est] <- parvec
  pars <- matrix(parmat$value, nrow = n_items, byrow = TRUE)

  maxk <- (ncol(pars) - 2) / 2

  xi <- pars[, 1]
  omega <- pars[, 2]

  if (maxk > 0){
    alpha <- as.matrix(pars[, c(1 + 2 * (1:maxk))]) # 3, 5, 7, ...
    tau <- as.matrix(pars[, c(2 + 2 * (1:maxk))]) # 4, 6, 8, ...

    if (n_items == 1){
      alpha <- matrix(alpha, nrow = 1)
      tau <- matrix(tau, nrow = 1)
    }
  }

  bmat <- matrix(nrow = n_items, ncol = 2 * maxk + 2)

  if (maxk == 0){
    for (it in 1:n_items){
      bmat[it, ] <- greek2b(xi = xi[it], omega = omega[it],
                           alpha = NULL, tau = NULL)
    }
  }

  if (maxk > 0){
    for (it in 1:n_items){
      bmat[it, ] <- greek2b(xi = xi[it], omega = omega[it],
                           alpha = alpha[it, ], tau = tau[it, ])
    }
  }

  probs <- irf_fmp(theta = quad_nodes, bmat = bmat)

  ## avoid probs of 1 or 0
  probs[probs == 1] <- 1 - 1e-16
  probs[probs == 0] <- 1e-16

  out <- matrix(nrow = n_items, ncol = ncol(pars))

  for (it in 1:n_items){
    part1 <- (r_bar[, it] - n_bar * probs[, it]) # quad_nodes
    part2 <- t(partial_m_mat <- sapply(quad_nodes, partial_m,
                                       parveci = pars[it, ], maxk = maxk))
    out[it, ] <- -t(part1) %*% part2
  }

  estmat <- matrix(parmat$est, nrow = n_items, byrow = TRUE)

  as.numeric(t(out)[t(estmat)])
}

find_t <- function(j, alpha, tau){
  out <- matrix(0, nrow = 2 * j + 1, ncol = 2 * j - 1)
  diag(out) <- 1

  if (j > 1){
    diag(out[-1, ]) <- -2 * alpha[j]
    diag(out[-c(1, 2), ]) <- alpha[j] ^ 2 + exp(tau[j])
  }

  if (j == 1){
    out[2, ] <- -2 * alpha[j]
    out[3, ] <- alpha[j] ^ 2 + exp(tau[j])
  }

  out
}

find_t_general <- function(j, el1, el2, el3){
  out <- matrix(0, nrow = 2 * j + 1, ncol = 2 * j - 1)
  diag(out) <- el1

  if (j > 1){
    diag(out[-1, ]) <- el2
    diag(out[-c(1, 2), ]) <- el3
  }

  if (j == 1){
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

  for (i in 1:ncol(v)) v[i:(i + length(wvec) - 1), i] <- wvec

  v
}
