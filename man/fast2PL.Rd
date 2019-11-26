\name{fast2PL}
\alias{fast2PL}

\alias{coef.twopl}
\alias{confint.twopl}
\alias{estfun.twopl}
\alias{logLik.twopl}
\alias{nobs.twopl}
\alias{print.summary.twopl}
\alias{print.twopl}
\alias{summary.twopl}
\alias{vcov.twopl}
\alias{weights.twopl}

\title{Single or Multiple Group 2PL Fitting Function}

\description{
  \code{fast2PL} is a basic fitting function for single and multiple group 2PL
  (in slope/intercept parametrization) under the MML framework relying on the
  application of an EM algorithm.
}

\usage{
fast2PL(y, weights = NULL, impact = NULL, start = NULL,
  control = NULL, algo_settings = NULL)
}

\arguments{
  \item{y}{item response object that can be coerced (via
    \code{\link[base]{as.matrix}}) to a numeric matrix with scores 0, or 1.
    Typically, either already a matrix or a data frame.}
  \item{weights}{an optional vector of weights (interpreted as case weights).}
  \item{impact}{an optional \code{factor} allowing for grouping the subjects
  (rows). If specified, a multiple group model is fitted to account for impact.
  By default, no impact is modelled, i.e., a single group model is fitted.}
  \item{start}{an optional list of starting values.}
  \item{control}{an optional list of control parameters for the EM algorithm.}
  \item{algo_settings}{an optional list of control parameters for the numerical
    optimization method.}
}

\details{
  \code{fast2PL} returns an object of class \code{"twopl"} for which several
  basic methods are available, including \code{print}, \code{summary},
  \code{coef}, \code{vcov} (OPG estimate), \code{logLik}, and \code{estfun}.

  Optimization is performed relying on \code{OptimLib}.

  \code{control} allows for the following settings (their defaults being shown
    in the first place):
  \itemize{
  \item{optimizer = "BFGS", "Newton", "L-BFGS", "CG", "GD", "NM", "DE", "PSO"}
  \item{accelerator = "Zhou", "Ramsay", "SQUAREM", "none"}
  \item{maxit = 500}
  \item{reltol = 1e-4}
  \item{Q = 61}
  \item{global = FALSE}
  \item{criterion = "ll", "l2", "l2_itemwise"}
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
  \code{fast2PL} returns an S3 object of class \code{"twopl"}, i.e., a list of
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
  \item{rgl}{number of occurences of each pattern (row) in y_u splitted for each
    level of \code{impact} (if supplied),}
  \item{X}{quadrature nodes,}
  \item{AX}{quadrature weights splitted for each level of \code{impact} (if
    supplied),}
  \item{time}{estimation time in nanoseconds, broken into different sections.}
}

\references{
  Bock RD, Aitkin M (1981).
  Marginal Maximum Likelihood Estimation of Item Parameters: Application of
  an EM Algorithm.
  \emph{Psychometrika}, \bold{46}(4), 443--459.

  O'Hara K (2019).
  OptimLib: A lightweight C++ library of numerical optimization methods for
  nonlinear functions.
  \url{https://www.kthohr.com/optimlib.html}.

  Ramsay JO (1975).
  Solving Implicit Equations in Psychometric Data Analysis.
  \emph{Psychometrika}, \bold{40}(3), 337--360.

  Varadhan R, Roland C (2008).
  Simple and Globally Convergent Methods for Accelerating the Convergence of
  any EM Algorithm.
  \emph{Scandinavian Journal of Statistics}, \bold{35}(2), 335--353.

  Zhou H, Alexander D, Lange K (2011).
  A Quasi-Newton Acceleration for High-Dimensional Optimization Algorithms.
  \emph{Statistics and Computing}, \bold{21}(2), 261--273.
}

\examples{
  data(LSAT7, package = "sonic")
  m1 <- fast2PL(LSAT7)
}
