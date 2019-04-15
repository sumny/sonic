### FIXME: check defaults
fit <- function(y, M, N, weights, impact, start, model, control, algo_settings)
{

  ### model (starting from 0)
  model <- as.integer(factor(match.arg(model,
    c("2PL", "RM", "3PL", "3PLu", "4PL")),
    levels = c("2PL", "RM", "3PL", "3PLu", "4PL"))) - 1L

  ### EM and general control settings
  Rcontrol <- list(optimizer = 1L, accelerator = 3L, maxit = 500L, reltol = 1e-4, Q = 61L, global = FALSE, criterion = 0L)

  ### optimizer (starting from 0)
  control$optimizer <- as.integer(factor(match.arg(control$optimizer,
    c("Newton", "BFGS", "L-BFGS", "CG", "GD", "NM", "DE", "PSO")),
    levels = c("Newton", "BFGS", "L-BFGS", "CG", "GD", "NM", "DE", "PSO"))) - 1L

  ### accelerator (starting from 0)
  control$accelerator <- as.integer(factor(match.arg(control$accelerator,
    c("none", "Ramsay", "SQUAREM", "Zhou")),
    levels = c("none", "Ramsay", "SQUAREM", "Zhou"))) - 1L

  ### maxit, reltol, Q, global

  ### criterion (starting from 0)
  control$criterion <- as.integer(factor(match.arg(control$criterion,
    c("ll", "l2", "l2_itemwise")),
    levels = c("ll", "l2", "l2_itemwise"))) - 1L

  ### Newton currently only 2PL or RM
  Rcontrol[names(control)] <- control
  if(((model != 0) && (model != 1)) && Rcontrol$optimizer == 0) {
    warning("'optimizer = Newton' currently only supported for 2PL or RM. Switched to 'BFGS'.")
    Rcontrol$optimizer <- 1
  }

  ### global currently only 2PL - always worse than local so drop this soon
  if((model != 0) && (Rcontrol$global == TRUE)) {
    warning("'global = TRUE' currently only supported for 2PL. Set to 'FALSE'.")
    Rcontrol$global <- FALSE
  }

  if(control$criterion == 2) {
    warning("'criterion = l2_itemwise' requires 'global = FALSE'.")
    Rcontrol$global <- FALSE
  }

  ### optim algo settings
  Ralgo_settings <- list(err_tol = 1e-8, iter_max = 25L, lbfgs_par_M = 10L,
    cg_method = c("FR", "PR", "FR-PR", "HS", "DY", "HZ"), cg_restart_threshold = 0.1,
    gd_method = c("Basic", "Momentum", "NAG", "AdaGrad", "RMSprop", "AdaDelta", "Adam/AdaMax", "Nadam/NadaMax"),
    nm_par_alpha = 1.0, nm_par_beta = 0.5, nm_par_gamma = 2.0, nm_par_delta = 0.5,
    gd_step_size = 0.1, gd_momentum_par = 0.9, gd_norm_term = 10e-08, gd_ada_rho = 0.9, gd_adam_beta_1 = 0.9, gd_adam_beta_2 = 0.999, gd_ada_max = FALSE,
    de_n_pop = 200L, de_n_gen = 1000L, de_check_freq = -1L, de_mutation_method = c("rand", "best"), de_par_F = 0.8, de_par_CR = 0.9,
    pso_n_pop = 100L, pso_n_gen = 1000L, pso_check_freq = -1L, pso_inertia_method = c("ld", "dampening"), pso_par_w_min = 0.1, pso_par_w_max = 0.99, pso_par_w_damp = 0.99,
    pso_velocity_method = c("fw", "ld"), pso_par_c_cog = 2.0, pso_par_c_soc = 2.0, pso_par_initial_c_cog = 2.5, pso_par_final_c_cog = 0.5, pso_par_initial_c_soc = 0.5, pso_par_final_c_soc = 2.5)

  ### CG method (starting from 1)
  algo_settings$cg_method <- as.integer(factor(match.arg(algo_settings$cg_method,
    c("FR", "PR", "FR-PR", "HS", "DY", "HZ")),
    levels = c("FR", "PR", "FR-PR", "HS", "DY", "HZ")))

  ### GD method (starting from 0)
  algo_settings$gd_method <- as.integer(factor(match.arg(algo_settings$gd_method,
    c("Basic", "Momentum", "NAG", "AdaGrad", "RMSprop", "AdaDelta", "Adam/AdaMax", "Nadam/NadaMax")),
    levels = c("Basic", "Momentum", "NAG", "AdaGrad", "RMSprop", "AdaDelta", "Adam/AdaMax", "Nadam/NadaMax"))) - 1L

  ### DE mutation method (starting from 1)
  algo_settings$de_mutation_method <- as.integer(factor(match.arg(algo_settings$de_mutation_method,
    c("rand", "best")),
    levels = c("rand", "best")))

  ### PSO inertia method (starting from 1)
  algo_settings$pso_inertia_method <- as.integer(factor(match.arg(algo_settings$pso_inertia_method,
    c("ld", "dampening")),
    levels = c("ld", "dampening")))

  ### PSO velocity method (starting from 1)
  algo_settings$pso_velocity_method <- as.integer(factor(match.arg(algo_settings$pso_velocity_method,
    c("fw", "ld")),
    levels = c("fw", "ld")))

  Ralgo_settings[names(algo_settings)] <- algo_settings

  ### fit
  fit <- .Call("fit", y, weights, as.integer(impact) - 1, start, model, Rcontrol, Ralgo_settings)

  return(fit)
}
