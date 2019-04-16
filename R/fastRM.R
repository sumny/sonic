### fitting function
fastRM <- function(y, weights = NULL, impact = NULL, start = NULL,
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
    ds <- qnorm(colMeans(y)) * 1.95
    tsuf <- rowSums(y)
    xs <- aggregate(tsuf ~ impact, FUN = mean)$tsuf
    vs <- aggregate(tsuf ~ impact, FUN = var)$tsuf
    xs <- xs - xs[1]
    vs <- vs / vs[1]
    list(ds, xs, vs)
  } else {
    start
  }

  fit <- fit(y, M, N, weights, impact, start, "RM", control, algo_settings)
  if(!fit$converged) warning("EM algorithm did not converge.")
  fit$y <- y
  fit$weights <- weights
  fit$impact <- impact
  fit$call <- cl
  class(fit) <- "rm"

  return(fit)
}

### methods
coef.rm <- function(object, ...)
{
  pars <- as.vector(object$ipars)
  N <- length(pars)
  pars <- c(pars, object$sg[1L])
  names(pars) <- c(paste0("d_", 1L:N), paste0("sg_", levels(object$impact)[1L]))
  G <- dim(object$mu)[1L]
  if(G > 1L) {
    gpars <- c(object$mu[-1L], object$sg[-1L])
    names(gpars) <- paste0(c("mu_", "sg_"), levels(object$impact)[-1L])
    pars <- c(pars, gpars)
  }

  return(pars)
}

vcov.rm <- function(object, ...)
{
  sc <- estfun.rm(object, weights = sqrt(weights(object)))
  vc <- chol2inv(chol(crossprod(sc)))
  rownames(vc) <- colnames(vc) <- names(coef(object))

  return(vc)
}

logLik.rm <- function(object, ...)
{
  return(structure(object$ll, df = length(coef(object)), class = "logLik"))
}

weights.rm <- function(object, ...)
{
  return(object$weights)
}

print.rm <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("RM model parameters:\n")
  print(coef(x), digits = digits)
  invisible(x)
}

summary.rm <- function(object, vcov. = NULL, ...)
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
  class(object) <- "summary.rm"

  return(object)
}

print.summary.rm <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("\nCall:\n")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat("RM model parameters and standard errors:\n")
  print(x$coefficients, digits = digits)

  cat("\nLog-likelihood:", format(signif(x$ll, digits)),
    "(df =", paste(length(x$coefficients[, 1]), ")", sep = ""), "\n")
  cat("Converged:", if(x$converged == 1) "yes" else "no", "\n")
  cat("Termination criterion value:", format(signif(x$critval, digits)), "\n")
  cat("Number of EM cycles:", x$iter, "\n")
  cat("Elapsed estimation time (s):", format(signif(x$time[length(x$time)] * 1e-9, digits)), "\n")

  invisible(x)
}

nobs.rm <- function (object, ...)
{
  return(dim(object$y)[1])
}

### FIXME: push this in C++
estfun.rm <- function(x, weights = x$weights, ...)
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
  dest <- coefs[1:N]
  pX <- LX <- matrix(0, M, Q)
  scores_d <- vector("list", G)
  nest <- N
  scores <- matrix(0, M, nest)
  for(g in 1L:G) {
    px_tmp_ <- matrix(1, sum(impact == g), Q)
    scores_d[[g]] <- vector("list", N)
    M_tmp <- sum(impact == g)
    for(j in 1L:N) {
      pat <- dat[impact == g, j] + 1
      dest_tmp <- dest[j]
      exp_tmp <- exp(X + dest_tmp)
      px_tmp <- exp_tmp / (1 + exp_tmp)
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
      expone <- 1 + exp_tmp
      exponesq <- expone ^ 2
      px_dd_tmp <- exp_tmp / exponesq
      px_dd_tmp <- rbind(-px_dd_tmp, px_dd_tmp)
      scores_d[[g]][[j]] <- px_dd_tmp[pat, ] / px_tmp[pat, ]
    }
    LX[impact == g, ] <- px_tmp_
    pX[impact == g, ] <- px_tmp_ * matrix(AX[, g], M_tmp, Q, TRUE)
    pX[impact == g, ] <- pX[impact == g, ] / rowSums(pX[impact == g, , drop = FALSE])
    scores_d[[g]] <-
    lapply(scores_d[[g]], function(sc_d) {
      rowSums(sc_d * pX[impact == g, , drop = FALSE])
    })
    text <- paste0("cbind(", paste0(
      sapply(1L:N, function(j) {
        gsub("item", j, "scores_d[[g]][[item]]")
      }), collapse = ","), ")"
    )
    scores[impact == g, ] <- eval(parse(text = text))
  }
  PX <- numeric(M)
  for(g in 1L:G) {
    PX[impact == g] <- rowSums(LX[impact == g, ] * matrix(AX[, g], sum(impact == g), Q, TRUE))
  }
  for(g in 1L:G) {
    if(g == 1) {
      scores <- cbind(scores, matrix(0, M, 1L))
      m <- means[g]
      v <- vars[g]
      scores[impact == g, dim(scores)[2L]] <- (- (1 / 2) * (PX[impact == g] ^ -1) *
        colSums(matrix((v ^ -1) - ((X - m) ^ 2) * (v ^ -2), Q, sum(impact == g)) *
          t(LX[impact == g, , drop = FALSE]) * matrix(AX[, g], Q, sum(impact == g))))
    } else {
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
