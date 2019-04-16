\name{fast3PLu}
\alias{fast3PLu}

\alias{coef.threeplu}
\alias{logLik.threeplu}
\alias{nobs.threeplu}
\alias{print.threeplu}
\alias{print.summary.threeplu}
\alias{summary.threeplu}
\alias{vcov.threeplu}
\alias{weights.threeplu}
\alias{confint.threeplu}
\alias{estfun.threeplu}

\title{Fit a Single or Multiple Group 3PLu}

\description{
  FIXME
}

\usage{
fast3PLu(y, weights = NULL, impact = NULL, start = NULL,
  control = NULL, algo_settings = NULL)
}

\arguments{
  \item{y}{item response object that can be coerced (via
    \code{\link[base]{as.matrix}}) to a numeric matrix with scores 0, 1.
    Typically, either already a matrix or a data frame.}
  \item{weights}{an optional vector of weights (interpreted as case weights).}
  \item{impact}{an optional \code{factor} allowing for grouping the subjects
    (rows). If specified, a multiple group model is fitted to account for impact
    (see details below). By default, no impact is modelled, i.e., a single group
    model is fitted.}
  \item{start}{an optional list of starting values (see examples below).}
  \item{control}{an optional list of control parameters for the EM algorithm.}
  \item{algo_settings}{an optional list of control parameters for the numerical
    optimization method.}
}

\details{
  \code{control} allows for the following settings (their defaults being shown
    in the first place):
  \itemize{
  \item{optimizer = "BFGS", "Newton", "L-BFGS", "CG", "GD", "NM", "DE", "PSO"}
  \item{accelerator = "Zhou", "none", "Ramsay", "SQUAREM"}
  \item{maxit = 500}
  \item{reltol = 1e-4}
  \item{Q = 61}
  }
  \code{algo_settings} allows for the following settings (their defaults being
    shown in the first place):
  \itemize{
  \item{err_tol = 1e-8}
  \item{iter_max = 25}
  \item{lbfgs_par_M = 10}
  \item{cg_method = "FR", "PR", "FR-PR", "HS", "DY", "HZ"}
  \item{cg_restart_threshold = 0.1}
  \item{gd_method = "Basic", "Momentum", "NAG", "AdaGrad", "RMSprop",
    "AdaDelta", "Adam/AdaMax", "Nadam/NadaMax"}
  \item{nm_par_alpha = 1.0}
  \item{nm_par_beta = 0.5}
  \item{nm_par_gamma = 2.0}
  \item{nm_par_delta = 0.5}
  \item{gd_step_size = 0.1}
  \item{gd_momentum_par = 0.9}
  \item{gd_norm_term = 10e-08}
  \item{gd_ada_rho = 0.9}
  \item{gd_adam_beta_1 = 0.9}
  \item{gd_adam_beta_2 = 0.999}
  \item{gd_ada_max = FALSE}
  \item{de_n_pop = 200}
  \item{de_n_gen = 1000}
  \item{de_check_freq = -1}
  \item{de_mutation_method = "rand", "best"}
  \item{de_par_F = 0.8}
  \item{de_par_CR = 0.9}
  \item{pso_n_pop = 100}
  \item{pso_n_gen = 1000}
  \item{pso_check_freq = -1}
  \item{pso_inertia_method = "ld", "dampening"}
  \item{pso_par_w_min = 0.1}
  \item{pso_par_w_max = 0.99}
  \item{pso_par_w_damp = 0.99}
  \item{pso_velocity_method = "fw", "ld"}
  \item{pso_par_c_cog = 2.0}
  \item{pso_par_c_soc = 2.0}
  \item{pso_par_initial_c_cog = 2.5}
  \item{pso_par_final_c_cog = 0.5}
  \item{pso_par_initial_c_soc = 0.5}
  \item{pso_par_final_c_soc = 2.5}
  }
}

\value{
  \code{fast3PLu} returns an S3 object of class \code{"threeplu"},  i.e., a list of
    the following components:
  \item{ipars}{estimated item parameters in slope/intercept parametrization,}
  \item{mu}{(estimated) mean(s) of the normal distribution(s) for the person
    parameters (for each level of \code{impact}),}
  \item{sg}{(estimated) standard deviation(s) of the normal distribution(s) for
    the person parameters (for each level of \code{impact}),}
  \item{ll}{marginal log-likelihood of the fitted model,}
  \item{iter}{number of EM cycles,}
  \item{converged}{convergence code (0 or 1),}
  \item{convergence}{logical vector indicating convergence in each M-step,}
  \item{critval}{critical value indicating termination of the EM,}
  \item{y_u}{processed item response matrix containing unique rows (possibly
    reordered),}
  \item{rgl}{number of occurences of each pattern (row) in y_ splitted for each
    level of \code{impact} (if supplied),}
  \item{X}{quadrature nodes,}
  \item{AX}{quadrature weights splitted for each level of \code{impact} (if
    supplied),}
  \item{time}{estimation time in nanoseconds, broken into different sections.}
}

\references{
  FIXME
}

%\seealso{}

%\examples{}
