### fitting function
fast3PL <- function(y, weights = NULL, impact = NULL, start = NULL,
  control = list(optimizer = c("Newton", "BFGS", "L-BFGS", "CG", "GD", "NM", "DE", "PSO"),
    accelerator = c("none", "Ramsay", "SQUAREM", "Zhou"), maxit = 500L, reltol = 1e-4, Q = 61L,
    global = FALSE, criterion = c("ll", "l2", "l2_itemwise")),
  algo_settings = list(err_tol = 1e-8, iter_max = 2000L, lbfgs_par_M = 10L,
    cg_method = c("FR", "PR", "FR-PR", "HS", "DY", "HZ"), cg_restart_threshold = 0.1,
    gd_method = c("Basic", "Momentum", "NAG", "AdaGrad", "RMSprop", "AdaDelta", "Adam/AdaMax", "Nadam/NadaMax"),
    nm_par_alpha = 1.0, nm_par_beta = 0.5, nm_par_gamma = 2.0, nm_par_delta = 0.5,
    gd_step_size = 0.1, gd_momentum_par = 0.9, gd_norm_term = 10e-08, gd_ada_rho = 0.9,
    gd_adam_beta_1 = 0.9, gd_adam_beta_2 = 0.999, gd_ada_max = FALSE, de_n_pop = 200L,
    de_n_gen = 1000L, de_check_freq = -1L, de_mutation_method = c("rand", "best"), de_par_F = 0.8,
    de_par_CR = 0.9, pso_n_pop = 100L, pso_n_gen = 1000L, pso_check_freq = -1L,
    pso_inertia_method = c("ld", "dampening"), pso_par_w_min = 0.1, pso_par_w_max = 0.99,
    pso_par_w_damp = 0.99, pso_velocity_method = c("fw", "ld"), pso_par_c_cog = 2.0,
    pso_par_c_soc = 2.0, pso_par_initial_c_cog = 2.5, pso_par_final_c_cog = 0.5,
    pso_par_initial_c_soc = 0.5, pso_par_final_c_soc = 2.5))
{
  cl <- match.call()
  y <- as.matrix(y)
  M <- dim(y)[1]
  N <- dim(y)[2]

  weights <- if(is.null(weights)) {
    rep(1, M)
  } else {
    weights
  }

  impact <- if(is.null(impact)) {
    as.factor(rep("all", M))
  } else {
    as.factor(impact)
  }

  ### FIXME: start as vector
  start <- if(is.null(start)) {
    as <- rep(0.851, N)
    ds <- qnorm(colMeans(y)) * 1.95
    gs <- rep(0.15, N)
    gs <- log(gs / (1 - gs))
    tsuf <- rowSums(y * matrix(as, M, N, byrow = TRUE))
    xs <- aggregate(tsuf ~ impact, FUN = mean)$tsuf
    vs <- aggregate(tsuf ~ impact, FUN = var)$tsuf
    xs <- xs - xs[1]
    vs <- vs / vs[1]
    list(c(as, ds, gs), xs, vs)
  } else {
    start
  }

  fit <- fit(y, M, N, weights, impact, start, "3PL", control, algo_settings)
  if(!fit$converged) warning("EM algorithm did not converge.")
  fit$y <- y
  fit$weights <- weights
  fit$impact <- impact
  fit$call <- cl
  class(fit) <- "threepl"

  return(fit)
}

### methods
coef.threepl <- function(object, ...)
{
  pars <- as.vector(object$ipars)
  N <- length(pars) / 3L
  names(pars) <- c(paste0("a_", 1L:N), paste0("d_", 1L:N), paste0("g_", 1L:N))
  G <- dim(object$mu)[1L]
  if(G > 1L) {
    gpars <- c(object$mu[-1L], object$sg[-1L])
    names(gpars) <- paste0(c("mu_", "sg_"), levels(object$impact)[-1L])
    pars <- c(pars, gpars)
  }

  return(pars)
}

vcov.threepl <- function(object, ...)
{
  sc <- estfun.threepl(object, weights = sqrt(weights(object)))
  vc <- chol2inv(chol(crossprod(sc)))
  rownames(vc) <- colnames(vc) <- names(coef(object))

  return(vc)
}

logLik.threepl <- function(object, ...)
{
  return(structure(object$ll, df = length(coef(object)), class = "logLik"))
}

weights.threepl <- function(object, ...)
{
  return(object$weights)
}

print.threepl <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("3PL model parameters:\n")
  print(coef(x), digits = digits)
  invisible(x)
}

summary.threepl <- function(object, vcov. = NULL, ...)
{
  cf <- coef(object)
  if(is.null(vcov.))
    vc <- vcov(object)
  else {
    if(is.function(vcov.)) vc <- vcov.(object) else vc <- vcov.
  }
  cf <- cbind(cf, sqrt(diag(vc)))
  colnames(cf) <- c("Estimate", "Std. Error")

  object$coefficients <- cf
  class(object) <- "summary.threepl"

  return(object)
}

print.summary.threepl <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("\nCall:\n")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat("3PL model parameters and standard errors:\n")
  print(x$coefficients, digits = digits)

  cat("\nLog-likelihood:", format(signif(x$ll, digits)),
    "(df =", paste(length(x$coefficients[, 1]), ")", sep = ""), "\n")
  cat("Converged:", if(x$converged == 1) "yes" else "no", "\n")
  cat("Termination criterion value:", format(signif(x$critval, digits)), "\n")
  cat("Number of EM cycles:", x$iter, "\n")
  cat("Elapsed estimation time (s):", format(signif(x$time[length(x$time)] * 1e-9, digits)), "\n")

  invisible(x)
}

nobs.threepl <- function (object, ...)
{
  return(dim(object$y)[1])
}

### FIXME: push this in C++
estfun.threepl <- function(x, weights = x$weights, ...)
{
  dat <- x$y
  impact <- x$impact
  G <- length(levels(impact))
  impact <- as.numeric(impact)
  coefs <- coef(x)
  estnms <- names(coefs)
  X <- as.vector(x$X)
  Q <- length(X)
  AX <- x$AX
  M <- dim(dat)[1]
  N <- dim(dat)[2]
  p <- rep.int(2L, N)
  K <- 2L
  means <- x$mu
  vars <- x$sg ^ 2
  aest <- coefs[1:N]
  dest <- coefs[(N + 1):(2 * N)]
  gest <- plogis(coefs[((2 * N) + 1):(3 * N)])
  pX <- LX <- matrix(0, M, Q)
  scores_a <- scores_d <- scores_g <- vector("list", G)
  nest <- 3L * N
  scores <- matrix(0, M, nest)
  for(g in 1L:G) {
    px_tmp_ <- matrix(1, sum(impact == g), Q)
    scores_a[[g]] <- scores_d[[g]] <- scores_g[[g]] <- vector("list", N)
    M_tmp <- sum(impact == g)
    for(j in 1L:N) {
      pat <- dat[impact == g, j] + 1
      aest_tmp <- aest[j]
      dest_tmp <- dest[j]
      gest_tmp <- gest[j]
      exp_tmp <- exp(aest_tmp * X + dest_tmp)
      px_tmp <- gest_tmp + (1 - gest_tmp) * (exp_tmp / (1 + exp_tmp))
      if(any(px_tmp == 1)) {
        px_tmp[px_tmp == 1] <- 1 - .Machine$double.neg.eps
      }
      if(any(px_tmp == 0)) {
        px_tmp[px_tmp == 0] <- .Machine$double.neg.eps
      }
      px_tmp <- rbind(1 - px_tmp, px_tmp)
      px_tmp_sel <- px_tmp[pat, ]
      ex <- is.na(px_tmp_sel)
      if(any(ex)) {
        px_tmp_sel[which(ex, arr.ind = TRUE)] <- 1
      }
      px_tmp_ <- px_tmp_ * px_tmp_sel
      gdiffu <- gest_tmp - 1
      expone <- 1 + exp_tmp
      exponesq <- expone ^ 2
      px_da_tmp <- -(gdiffu * X * exp_tmp) / exponesq
      px_da_tmp <- rbind(-px_da_tmp, px_da_tmp)
      scores_a[[g]][[j]] <- px_da_tmp[pat, ] / px_tmp[pat, ]
      px_dd_tmp <- -(gdiffu * exp_tmp) / exponesq
      px_dd_tmp <- rbind(-px_dd_tmp, px_dd_tmp)
      scores_d[[g]][[j]] <- px_dd_tmp[pat, ] / px_tmp[pat, ]
      gest_tmp_logit <- qlogis(gest_tmp)
      px_dg_tmp <- exp(gest_tmp_logit) / (((1 + exp(gest_tmp_logit)) ^ 2) * expone)
      px_dg_tmp <- rbind(-px_dg_tmp, px_dg_tmp)
      scores_g[[g]][[j]] <- px_dg_tmp[pat, ] / px_tmp[pat, ]
    }
    LX[impact == g, ] <- px_tmp_
    pX[impact == g, ] <- px_tmp_ * matrix(AX[, g], M_tmp, Q, TRUE)
    pX[impact == g, ] <- pX[impact == g, ] / rowSums(pX[impact == g, , drop = FALSE])
    scores_a[[g]] <-
    lapply(scores_a[[g]], function(sc_a) {
      rowSums(sc_a * pX[impact == g, , drop = FALSE])
    })
    scores_d[[g]] <-
    lapply(scores_d[[g]], function(sc_d) {
      rowSums(sc_d * pX[impact == g, , drop = FALSE])
    })
    scores_g[[g]] <-
     lapply(scores_g[[g]], function(sc_g) {
       rowSums(sc_g * pX[impact == g, , drop = FALSE])
    })
    text <- paste0("cbind(", paste0(
      c(sapply(1L:N, function(j) {
        gsub("item", j, "scores_a[[g]][[item]]")
      }),
      sapply(1L:N, function(j) {
        gsub("item", j, "scores_d[[g]][[item]]")
      }),
      sapply(1L:N, function(j) {
        gsub("item", j, "scores_g[[g]][[item]]")
      })), collapse = ","), ")"
    )
    scores[impact == g, ] <- eval(parse(text = text))
  }
  if(G > 1L) {
    PX <- numeric(M)
    for(g in 2L:G) {
      PX[impact == g] <- rowSums(LX[impact == g, ] * matrix(AX[, g], sum(impact == g), Q, TRUE))
    }
    for(g in 2L:G) {
      scores <- cbind(scores, matrix(0, M, 2L))
      m <- means[g]
      v <- vars[g]
      scores[impact == g, dim(scores)[2L] - 1L] <- ((v ^ -1) * (PX[impact == g] ^ -1) *
        colSums(matrix(X - m, Q, sum(impact == g)) * t(LX[impact == g, , drop = FALSE]) *
          matrix(AX[, g], Q, sum(impact == g))))
      scores[impact == g, dim(scores)[2L]] <- (- (1 / 2) * (PX[impact == g] ^ -1) *
        colSums(matrix((v ^ -1) - ((X - m) ^ 2) * (v ^ -2), Q, sum(impact == g)) *
          t(LX[impact == g, , drop = FALSE]) * matrix(AX[, g], Q, sum(impact == g))))
    }
  }
  scores <- scores * weights
  scores[which(is.na(scores), arr.ind = TRUE)] <- 0
  colnames(scores) <- estnms

  return(scores)
}
