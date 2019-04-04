% FIXME
\name{fast2PL}
\alias{fast2PL}

%\alias{}

\title{Fit a single-group or special multiple-group 2PL}

\description{
  FIXME
}

\usage{
fast2PL(y, weights = NULL, impact = NULL, start = NULL,
  control = list(optimizer = c("Newton", "BFGS", "L-BFGS-B", "CG", "GD", "NM", "DE", "PSO"), accelerator = c("none", "Ramsay", "SQUAREM", "Zhou", "Anderson"), maxit = 500L, reltol = 1e-4, Q = 61L, global = FALSE, criterium = c("ll", "l2", "l2_itemwise")),
  algo_settings = list(err_tol = 1e-8, iter_max = 2000L, lbfgs_par_M = 10L,
    cg_method = c("FR", "PR", "FR-PR", "HS", "DY", "HZ"), cg_restart_threshold = 0.1,
    gd_method = c("Basic", "Momentum", "NAG", "AdaGrad", "RMSprop", "AdaDelta", "Adam/AdaMax", "Nadam/NadaMax"),
    nm_par_alpha = 1.0, nm_par_beta  = 0.5, nm_par_gamma = 2.0, nm_par_delta = 0.5,
    gd_step_size = 0.1, gd_momentum_par = 0.9, gd_norm_term = 10e-08, gd_ada_rho = 0.9, gd_adam_beta_1 = 0.9, gd_adam_beta_2 = 0.999, gd_ada_max = FALSE,
    de_n_pop = 200L, de_n_gen = 1000L, de_check_freq = -1L, de_mutation_method = c("rand", "best"), de_par_F = 0.8, de_par_CR = 0.9,
    pso_n_pop = 100L, pso_n_gen = 1000L, pso_check_freq = -1L, pso_inertia_method = c("ld", "dampening"), pso_par_w_min = 0.1, pso_par_w_max = 0.99, pso_par_w_damp = 0.99,
    pso_velocity_method = c("fw", "ld"), pso_par_c_cog = 2.0, pso_par_c_soc = 2.0, pso_par_initial_c_cog = 2.5, pso_par_final_c_cog = 0.5, pso_par_initial_c_soc = 0.5, pso_par_final_c_soc = 2.5))
}

\arguments{
  \item{y}{item response object that can be coerced
    (via \code{\link[base]{as.matrix}}) to a numeric matrix with scores 0, 1.
    Typically, either already a matrix or a data frame.}
  \item{weights}{an optional vector of weights (interpreted as case weights).}
  \item{impact}{an optional \code{factor} allowing for grouping the subjects
    (rows). If specified, a multiple-group model is fitted to account for
    impact (see details below). By default, no impact is modelled, i.e., a
    single-group model is fitted.}
  \item{start}{an optional list of starting values (see examples below).}
  \item{control}{an optional list of control parameters for the EM algorithm.}
  \item{algo_settings}{an optional list of control parameters for optim.}
}

\details{
  FIXME
}

\value{
  \code{fast2PL} returns an S3 object of class \code{"twopl"}, 
  i.e., a list of the following components:
  \item{ipars}{estimated item parameters in slope/intercept parametrization,}
  \item{mu}{(estimated) mean(s) of the normal distribution(s) for the person
    parameters (for each level of \code{impact}),}
  \item{sg}{(estimated) standard deviation(s) of the normal distribution(s) for
    the person parameters (for each level of \code{impact}),}
  \item{iter}{number of iterations until the EM algorithm converged,}
  \item{ll}{log-likelihood of the fitted model,}
  \item{converged}{convergence code (0 or 1),}
  \item{convergence}{vector of convergence with logicals indicating convergence
    in each EM iteration,}
  \item{critval}{critical value indicating termination of the EM,}
  \item{y_u}{processed item response matrix containing unique rows (possibly
    reordered),}
  \item{rgl}{number of occurences of each pattern (row) in y_ splitted for each
    level of \code{impact} (if supplied),}
  \item{time}{estimation time, broken into different sections.}
}

\references{
  FIXME
}

%\seealso{}

%\examples{ }

