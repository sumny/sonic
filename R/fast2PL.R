fast2PL <- function(y, weights = NULL, impact = NULL, start = NULL,
  control = list(optimizer = c("BFGS", "L-BFGS-B", "CG", "GD", "NM", "DE", "PSO"), accelerator = c("none", "Ramsay", "SQUAREM"), maxit = 500L, reltol = 1e-4, Q = 61L, global = FALSE, criterium = c("logL", "L2")),
  algo_settings = list(err_tol = 1e-8, iter_max = 2000L, lbfgs_par_M = 10L,
    cg_method = c("FR", "PR", "FR-PR", "HS", "DY", "HZ"), cg_restart_threshold = 0.1,
    gd_method = c("Basic", "Momentum", "NAG", "AdaGrad", "RMSprop", "AdaDelta", "Adam/AdaMax", "Nadam/NadaMax"),
    nm_par_alpha = 1.0, nm_par_beta  = 0.5, nm_par_gamma = 2.0, nm_par_delta = 0.5,
    de_n_pop = 200L, de_n_gen = 1000L, de_check_freq = -1L, de_mutation_method = c("rand", "best"), de_par_F = 0.8, de_par_CR = 0.9,
    pso_n_pop = 100L, pso_n_gen = 1000L, pso_check_freq = -1L, pso_inertia_method = c("ld", "dampening"), pso_par_w_min = 0.1, pso_par_w_max = 0.99, pso_par_w_damp = 0.99,
    pso_velocity_method = c("fw", "ld"), pso_par_c_cog = 2.0, pso_par_c_soc = 2.0, pso_par_initial_c_cog = 2.5, pso_par_final_c_cog = 0.5, pso_par_initial_c_soc = 0.5, pso_par_final_c_soc = 2.5),
  gd_settings = list(step_size = 0.1, momentum_par = 0.9, norm_term = 10e-08, ada_rho = 0.9, adam_beta_1 = 0.9, adam_beta_2 = 0.999, ada_max = FALSE))
{
  ## FIXME: push everything in C++ afterwards

  y <- as.matrix(y)
  M <- dim(y)[1]
  N <- dim(y)[2]

  ## weights
  weights <- if(is.null(weights)) {
    rep(1, M)
  } else {
    weights
  }

  ## impact
  impact <- if(is.null(impact)) {
    as.factor(rep("all", M))
  } else {
    impact
  }

  ## start
  start <- if(is.null(start)) {
    as <- rep(0.851, N)
    ds <- qnorm(colMeans(y)) * 1.95
    tsuf <- rowSums(y * matrix(as, M, N, byrow = TRUE))
    xs <- aggregate(tsuf ~ impact, FUN = mean)$tsuf
    vs <- aggregate(tsuf ~ impact, FUN = var)$tsuf
    xs <- xs - xs[1]
    vs <- vs / vs[1]
    list(c(as, ds), xs, vs)
  } else {
    start
  }

  ## EM and general control settings
  Rcontrol <- list(optimizer = 0L, accelerator = 0L, maxit = 500L, reltol = 1e-4, Q = 61L, global = FALSE, criterium = 0L)

  ## optimizer (starting from 0)
  control$optimizer <- as.integer(factor(match.arg(control$optimizer,
    c("BFGS", "L-BFGS-B", "CG", "GD", "NM", "DE", "PSO")),
    levels = c("BFGS", "L-BFGS-B", "CG", "GD", "NM", "DE", "PSO"))) - 1L

  ## accelerator (starting from 0)
  control$accelerator <- as.integer(factor(match.arg(control$accelerator,
    c("none", "Ramsay", "SQUAREM")),
    levels = c("none", "Ramsay", "SQUAREM"))) - 1L

  ## maxit, reltol, Q, global

  ## criterium (starting from 0)
  control$criterium <- as.integer(factor(match.arg(control$criterium,
    c("logL", "L2")),
    levels = c("logL", "L2"))) - 1L

  Rcontrol[names(control)] <- control

  ## optim algo settings
  Ralgo_settings <- list(err_tol = 1e-8, iter_max = 2000L, lbfgs_par_M = 10L, cg_method = 1L, cg_restart_threshold = 0.1,
    gd_method = 0L, nm_par_alpha = 1.0, nm_par_beta  = 0.5, nm_par_gamma = 2.0, nm_par_delta = 0.5,
    de_n_pop = 200L, de_n_gen = 1000L, de_check_freq = -1L, de_mutation_method = 1L, de_par_F = 0.8, de_par_CR = 0.9,
    pso_n_pop = 100L, pso_n_gen = 1000L, pso_check_freq = -1L, pso_inertia_method = 1L, pso_par_w_min = 0.1, pso_par_w_max = 0.99, pso_par_w_damp = 0.99,
    pso_velocity_method = 1L, pso_par_c_cog = 2.0, pso_par_c_soc = 2.0, pso_par_initial_c_cog = 2.5, pso_par_final_c_cog = 0.5, pso_par_initial_c_soc = 0.5, pso_par_final_c_soc = 2.5)
  ## CG method (starting from 1)
  algo_settings$cg_method <- as.integer(factor(match.arg(algo_settings$cg_method,
    c("FR", "PR", "FR-PR", "HS", "DY", "HZ")),
    levels = c("FR", "PR", "FR-PR", "HS", "DY", "HZ")))

  ## GD method (starting from 0)
  algo_settings$gd_method <- as.integer(factor(match.arg(algo_settings$gd_method,
    c("Basic", "Momentum", "NAG", "AdaGrad", "RMSprop", "AdaDelta", "Adam/AdaMax", "Nadam/NadaMax")),
    levels = c("Basic", "Momentum", "NAG", "AdaGrad", "RMSprop", "AdaDelta", "Adam/AdaMax", "Nadam/NadaMax"))) - 1L

  ## DE mutation method (starting from 1)
  algo_settings$de_mutation_method <- as.integer(factor(match.arg(algo_settings$de_mutation_method,
    c("rand", "best")),
    levels = c("rand", "best")))

  ## PSO inertia method (starting from 1)
  algo_settings$pso_inertia_method <- as.integer(factor(match.arg(algo_settings$pso_inertia_method,
    c("ld", "dampening")),
    levels = c("ld", "dampening")))

  ## PSO velocity method (starting from 1)
  algo_settings$pso_velocity_method <- as.integer(factor(match.arg(algo_settings$pso_velocity_method,
    c("fw", "ld")),
    levels = c("fw", "ld")))

  Ralgo_settings[names(algo_settings)] <- algo_settings

  ## GD settings
  ## FIXME: put this in algo_settings
  Rgd_settings <- list(step_size = 0.1, momentum_par = 0.9, norm_term = 10e-08, ada_rho = 0.9, adam_beta_1 = 0.9, adam_beta_2 = 0.999, ada_max = FALSE)

  Rgd_settings[names(gd_settings)] <- gd_settings

  ## fit
  fit <- .Call("fit", y, weights, as.integer(impact) - 1, start, Rcontrol, Ralgo_settings, Rgd_settings)

  class(fit) <- "fast2PL"

  return(fit)
}

