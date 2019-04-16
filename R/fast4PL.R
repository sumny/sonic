### fitting function
fast4PL <- function(y, weights = NULL, impact = NULL, start = NULL,
  control = NULL, algo_settings = NULL)
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
    us <- rep(0.85, N)
    us <- log(us / (1 - us))
    tsuf <- rowSums(y * matrix(as, M, N, byrow = TRUE))
    xs <- aggregate(tsuf ~ impact, FUN = mean)$tsuf
    vs <- aggregate(tsuf ~ impact, FUN = var)$tsuf
    xs <- xs - xs[1]
    vs <- vs / vs[1]
    list(c(as, ds, gs, us), xs, vs)
  } else {
    start
  }

  fit <- fit(y, M, N, weights, impact, start, "4PL", control, algo_settings)
  if(!fit$converged) warning("EM algorithm did not converge.")
  fit$y <- y
  fit$weights <- weights
  fit$impact <- impact
  fit$call <- cl
  class(fit) <- "fourpl"

  return(fit)
}

### methods
coef.fourpl <- function(object, ...)
{
  pars <- as.vector(object$ipars)
  N <- length(pars) / 4L
  names(pars) <- c(paste0("a_", 1L:N), paste0("d_", 1L:N), paste0("g_", 1L:N), paste0("u_", 1L:N))
  G <- dim(object$mu)[1L]
  if(G > 1L) {
    gpars <- c(object$mu[-1L], object$sg[-1L])
    names(gpars) <- paste0(c("mu_", "sg_"), levels(object$impact)[-1L])
    pars <- c(pars, gpars)
  }

  return(pars)
}

vcov.fourpl <- function(object, ...)
{
  sc <- estfun.fourpl(object, weights = sqrt(weights(object)))
  vc <- chol2inv(chol(crossprod(sc)))
  rownames(vc) <- colnames(vc) <- names(coef(object))

  return(vc)
}

logLik.fourpl <- function(object, ...)
{
  return(structure(object$ll, df = length(coef(object)), class = "logLik"))
}

weights.fourpl <- function(object, ...)
{
  return(object$weights)
}

print.fourpl <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("4PL model parameters:\n")
  print(coef(x), digits = digits)
  invisible(x)
}

summary.fourpl <- function(object, vcov. = NULL, ...)
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
  class(object) <- "summary.fourpl"

  return(object)
}

print.summary.fourpl <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("\nCall:\n")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat("4PL model parameters and standard errors:\n")
  print(x$coefficients, digits = digits)

  cat("\nLog-likelihood:", format(signif(x$ll, digits)),
    "(df =", paste(length(x$coefficients[, 1]), ")", sep = ""), "\n")
  cat("Converged:", if(x$converged == 1) "yes" else "no", "\n")
  cat("Termination criterion value:", format(signif(x$critval, digits)), "\n")
  cat("Number of EM cycles:", x$iter, "\n")
  cat("Elapsed estimation time (s):", format(signif(x$time[length(x$time)] * 1e-9, digits)), "\n")

  invisible(x)
}

nobs.fourpl <- function (object, ...)
{
  return(dim(object$y)[1])
}

### FIXME: push this in C++
estfun.fourpl <- function(x, weights = x$weights, ...)
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
  uest <- plogis(coefs[((3 * N) + 1):(4 * N)])
  pX <- LX <- matrix(0, M, Q)
  scores_a <- scores_d <- scores_g <- scores_u <- vector("list", G)
  nest <- 4L * N
  scores <- matrix(0, M, nest)
  for(g in 1L:G) {
    px_tmp_ <- matrix(1, sum(impact == g), Q)
    scores_a[[g]] <- scores_d[[g]] <- scores_g[[g]] <- scores_u[[g]] <- vector("list", N)
    M_tmp <- sum(impact == g)
    for(j in 1L:N) {
      pat <- dat[impact == g, j] + 1
      aest_tmp <- aest[j]
      dest_tmp <- dest[j]
      gest_tmp <- gest[j]
      uest_tmp <- uest[j]
      exp_tmp <- exp(aest_tmp * X + dest_tmp)
      px_tmp <- gest_tmp + (uest_tmp - gest_tmp) * (exp_tmp / (1 + exp_tmp))
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
      uest_tmp_logit <- qlogis(uest_tmp)
      px_du_tmp <- exp(aest_tmp * X + dest_tmp + uest_tmp_logit) /
        (((1 + exp(uest_tmp_logit)) ^ 2) * expone)
      px_du_tmp <- rbind(-px_du_tmp, px_du_tmp)
      scores_u[[g]][[j]] <- px_du_tmp[pat, ] / px_tmp[pat, ]
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
    scores_u[[g]] <-
     lapply(scores_u[[g]], function(sc_u) {
       rowSums(sc_u * pX[impact == g, , drop = FALSE])
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
      }),
      sapply(1L:N, function(j) {
        gsub("item", j, "scores_u[[g]][[item]]")
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
