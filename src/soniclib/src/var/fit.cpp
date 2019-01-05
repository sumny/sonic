//#include <omp.h>
#include "optim/optim.hpp"
#include "var/var.hpp"
#include <math.h>
#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>

// fit function
// FIXME only use pointers where needed
// FIXME use structs instead of single local variables
// FIXME check order of arguments in all functions etc. and rearrange if needed
// FIXME use std::move and arma::transform where possible
// FIXME implement adaptive Q, err_tol, iter_max, etc.
// FIXME arma::uword instead of int
// FIXME faster way for G = 1
// FIXME faster way if only unique patterns present
// FIXME implement pure Newton method
// FIXME reduce maxit = 0 time
// FIXME get OPENMP right; currently disabled
RcppExport SEXP fit(SEXP Ry, SEXP Rweights, SEXP Rimpact, SEXP Rstart, SEXP Rcontrol, SEXP Ralgo_settings, SEXP Rgd_settings)
{
  BEGIN_RCPP

  Rcpp::Timer timer;
  timer.step("start");

  // initial declarations
  arma::mat y = Rcpp::as<arma::mat>(Ry);
  arma::vec weights = Rcpp::as<arma::vec>(Rweights);
  // FIXME impact as uvec
  arma::vec impact = Rcpp::as<arma::vec>(Rimpact);
  arma::vec groups = arma::unique(impact);
  Rcpp::List start = Rcpp::as<Rcpp::List>(Rstart);
  Rcpp::List control = Rcpp::as<Rcpp::List>(Rcontrol);
  Rcpp::List algo_settings = Rcpp::as<Rcpp::List>(Ralgo_settings);
  Rcpp::List gd_settings = Rcpp::as<Rcpp::List>(Rgd_settings);
  int M = y.n_rows, N = y.n_cols, G = groups.n_elem;
  int optimizer = Rcpp::as<int>(control[0]), accelerator = Rcpp::as<int>(control[1]), maxit = Rcpp::as<int>(control[2]), Q = Rcpp::as<int>(control[4]), criterium = Rcpp::as<int>(control[6]);
  double reltol = Rcpp::as<double>(control[3]);
  bool global = Rcpp::as<bool>(control[5]);
  timer.step("declarations");

  // get the first N primes;
  arma::vec primes(N);
  sonic::generatePrimes(&N, &primes);
  timer.step("primes");

  // unique patterns
  int P;
  arma::mat rgl;
  arma::mat y_;
  sonic::preprocess(&y, &weights, &impact, &G, &primes, &P, &rgl, &y_);
  timer.step("unique_data");

  // start values (slopes and intercepts, group means and variances)
  arma::vec ipars = Rcpp::as<arma::vec>(start[0]);
  arma::vec ipars_i(2);
  arma::vec mu = Rcpp::as<arma::vec>(start[1]);
  arma::vec sg = arma::sqrt(Rcpp::as<arma::vec>(start[2]));

  // quadrature stuff
  // FIXME implement other types
  arma::vec X = arma::linspace<arma::vec>(-6, 6, Q);
  arma::mat AX(Q, G);
  for(int g = 0; g < G; ++g) {
    AX.col(g) = arma::normpdf(X, mu(g), sg(g));
    AX.col(g) = AX.col(g) / arma::accu(AX.col(g));
  }
  
  // some more declarations
  arma::vec ll_g(G);
  arma::cube rj_g(Q, 2 * N, G);
  arma::mat Prob(Q, N);
  arma::mat Probs(Q, 2 * N);
  arma::uvec p_vec = arma::linspace<arma::uvec>(0, P - 1, P);
  arma::uvec n_vec = arma::linspace<arma::uvec>(0, N - 1, N);
  arma::umat cat_ind(N, 2);
  cat_ind.each_col() = 2 * n_vec;
  cat_ind.col(1) += 1;
  arma::mat pul(Q, P);
  arma::mat pgul(P, G);
  double ll_old, ll_new, critval;
  int iter = 0;
  bool run = true;
  arma::uvec a_ind = arma::linspace<arma::uvec>(0, N - 1, N);
  arma::uvec d_ind = arma::linspace<arma::uvec>(N, 2 * N - 1, N);

  // Mstep stuff
  sonic::Mstep_data_g opt_data_g;
  sonic::Mstep_data_i opt_data_i;
  if(global) {
    opt_data_g.G = G;
    opt_data_g.ll_g = arma::vec(G);
    opt_data_g.a_ind = a_ind;
    opt_data_g.d_ind = d_ind;
    opt_data_g.gr_out = arma::vec(2 * N);
    opt_data_g.X = X;
    opt_data_g.gr_a = arma::mat(G, N);
    opt_data_g.gr_d = arma::mat(G, N);
    opt_data_g.gr_tmp = arma::mat(Q, N);
    opt_data_g.e = arma::mat(Q, N);
    opt_data_g.Prob = arma::mat(Q, N);
    opt_data_g.Probs = arma::mat(Q, 2 * N);
    opt_data_g.cat_ind = cat_ind;
  } else {
    opt_data_i.j = 0;
    opt_data_i.G = G;
    opt_data_i.gr_out = arma::vec(2);
    opt_data_i.gr_a = arma::vec(G);
    opt_data_i.gr_d = arma::vec(G);
    opt_data_i.gr_tmp = arma::vec(Q);
    opt_data_i.ll_g = ll_g;
    opt_data_i.X = X;
    opt_data_i.e = arma::vec(Q);
    opt_data_i.p = arma::mat(Q, 2);
  }

  // optim algo settings
  // FIXME constructor or if statements
  // FIXME Rcpp::as
  optim::algo_settings_t settings;
  settings.err_tol = algo_settings[0];
  settings.iter_max = algo_settings[1];
  settings.lbfgs_par_M = algo_settings[2];
  settings.cg_method = algo_settings[3];
  settings.cg_restart_threshold = algo_settings[4];
  settings.gd_method = algo_settings[5];
  settings.nm_par_alpha = algo_settings[6];
  settings.nm_par_beta = algo_settings[7];
  settings.nm_par_gamma = algo_settings[8];
  settings.nm_par_delta = algo_settings[9];
  settings.de_n_pop = algo_settings[10];
  settings.de_n_gen = algo_settings[11];
  settings.de_check_freq = algo_settings[12];
  settings.de_mutation_method = algo_settings[13];
  settings.de_par_F = algo_settings[14];
  settings.de_par_CR = algo_settings[15];
  settings.pso_n_pop = algo_settings[16];
  settings.pso_n_gen = algo_settings[17];

  // FIXME check this
  //settings.de_check_freq = algo_settings[18];

  settings.pso_inertia_method = algo_settings[19];
  settings.pso_par_w_min = algo_settings[20];
  settings.pso_par_w_max = algo_settings[21];
  settings.pso_par_w_damp = algo_settings[22];
  settings.pso_velocity_method = algo_settings[23];
  settings.pso_par_c_cog = algo_settings[24];
  settings.pso_par_c_soc = algo_settings[25];
  settings.pso_par_initial_c_cog = algo_settings[26];
  settings.pso_par_final_c_cog = algo_settings[27];
  settings.pso_par_initial_c_soc = algo_settings[28];
  settings.pso_par_final_c_soc = algo_settings[29];
  settings.gd_settings.step_size = gd_settings[0];
  settings.gd_settings.momentum_par = gd_settings[1];
  settings.gd_settings.norm_term = gd_settings[2];
  settings.gd_settings.ada_rho = gd_settings[3];
  settings.gd_settings.adam_beta_1 = gd_settings[4];
  settings.gd_settings.adam_beta_2 = gd_settings[5];
  settings.gd_settings.ada_max = gd_settings[6];

  timer.step("preE1");
  sonic::Estep(&y_, &rgl, &G, &ipars, &X, &AX, &ll_g, &rj_g, &Prob, &Probs, &pul, &pgul, &p_vec, &n_vec, &a_ind, &d_ind, &cat_ind);
  timer.step("postE1");

  // debug stuff
  //arma::vec gradient_g = arma::vec(2 * N);
  //arma::vec gradient_i = arma::vec(2);

  // group parameter stuff
  arma::mat grp_tmp(Q, P);
  arma::mat grp_X(Q, P);
  grp_X.each_col() = X;
  arma::mat grp_rgl(Q, P);
  grp_X.each_col() = X;

  // accelerator stuff
  arma::vec preM2 = ipars;
  arma::vec preM1 = preM2;
  arma::vec accels = preM2;
  double val;
  int sqr_iter = 0;
  arma::vec sqr_ipars = ipars;
  arma::vec sqr_ll_g = ll_g;
  arma::cube sqr_rj_g = rj_g;

  // EM
  if(maxit > 0) {
    while(run) {
      timer.step("preEM");
      Rcpp::checkUserInterrupt();
      iter += 1;
      preM2 = preM1;
      preM1 = ipars;
      ll_old = arma::accu(ll_g);

      // debug stuff
      //double ll = sonic::llfun_g(ipars, &gradient_g, &opt_data_g);
      //double ll = sonic::llfun_i(ipars, &gradient_i, &opt_data_i);
      //run = false;

      // item parameters
      // optimizer, 0 == "BFGS", 1 == "L-BFGS-B", 2 == "CG", 3 == "GD", 4 == "NM", 5 == "DE", 6 == "PSO"
      if(global) {
        opt_data_g.rj_g = rj_g;
        if(optimizer == 0) {
          bool success = optim::bfgs(ipars, sonic::llfun_g, &opt_data_g, settings);
        } else if(optimizer == 1) {
          bool success = optim::lbfgs(ipars, sonic::llfun_g, &opt_data_g, settings);
        } else if(optimizer == 2) {
          bool success = optim::cg(ipars, sonic::llfun_g, &opt_data_g, settings);
        } else if(optimizer == 3) {
          bool success = optim::gd(ipars, sonic::llfun_g, &opt_data_g, settings);
        } else if(optimizer == 4) {
          bool success = optim::nm(ipars, sonic::llfun_g, &opt_data_g, settings);
        } else if(optimizer == 5) {
          bool success = optim::de(ipars, sonic::llfun_g, &opt_data_g, settings);
        } else if(optimizer == 6) {
          bool success = optim::pso(ipars, sonic::llfun_g, &opt_data_g, settings);
        }
      } else {
        for(int j = 0; j < N; ++j) {
          opt_data_i.j = j;
          opt_data_i.rj_g = rj_g.cols(2 * j, 2 * j + 1);
          ipars_i(0) = ipars(j);
          ipars_i(1) = ipars(N + j);
          if(optimizer == 0) {
            bool success = optim::bfgs(ipars_i, sonic::llfun_i, &opt_data_i, settings);
          } else if(optimizer == 1) {
            bool success = optim::lbfgs(ipars_i, sonic::llfun_i, &opt_data_i, settings);
          } else if(optimizer == 2) {
            bool success = optim::cg(ipars_i, sonic::llfun_i, &opt_data_i, settings);
          } else if(optimizer == 3) {
            bool success = optim::gd(ipars_i, sonic::llfun_i, &opt_data_i, settings);
          } else if(optimizer == 4) {
            bool success = optim::nm(ipars_i, sonic::llfun_i, &opt_data_i, settings);
          } else if(optimizer == 5) {
            bool success = optim::de(ipars_i, sonic::llfun_i, &opt_data_i, settings);
          } else if(optimizer == 6) {
            bool success = optim::pso(ipars_i, sonic::llfun_i, &opt_data_i, settings);
          }
          ipars(j) = ipars_i(0);
          ipars(N + j) = ipars_i(1);
        }
      }

      // group parameters
      // FIXME flexible mu = 0 and sd = 1 for g == 0 (sum restr or fixed)
      // FIXME find mu and sg numerically optimized? Probably not needed
      if(G > 1) {
        for(int g = 1; g < G; ++g) {
          grp_tmp = pul;
          grp_tmp.each_col() %= AX.col(g);
          grp_tmp.each_row() /= arma::sum(grp_tmp, 0);
          grp_rgl.each_row() = rgl.col(g);
          mu(g) = (1 / arma::accu(rgl.col(g))) * arma::accu(grp_X % grp_tmp % grp_rgl);
          sg(g) = std::sqrt((1 / arma::accu(rgl.col(g))) * arma::accu(arma::square(grp_X - mu(g)) % grp_tmp % grp_rgl));
          AX.col(g) = arma::normpdf(X, mu(g), sg(g));
          AX.col(g) = AX.col(g) / arma::accu(AX.col(g));
        }
      }

      // acceleration
      // EM acceleration, 0 == "none", 1 == "Ramsay", 2 == "SQUAREM", FIXME 3 = "VA-Steffensen", ...
      // FIXME itemwise / parameterwise acceleration
      if((accelerator != 0) & std::fmod(iter, 3) == 0) {
        sonic::accelerate(&y_, &rgl, &G, &ipars, &X, &AX, &sqr_ll_g, &sqr_rj_g, &Prob, &Probs, &pul, &pgul, &p_vec, &n_vec, &a_ind, &d_ind, &cat_ind, &preM1, &preM2, &accels, &val, &sqr_iter, &accelerator, &ll_g);
      }

      sonic::Estep(&y_, &rgl, &G, &ipars, &X, &AX, &ll_g, &rj_g, &Prob, &Probs, &pul, &pgul, &p_vec, &n_vec, &a_ind, &d_ind, &cat_ind);
      ll_new = arma::accu(ll_g);

      // termination criterium, 0 == "logL", 1 == "L2"
      if(criterium == 0) {
        critval = ll_new - ll_old;
      } else if(criterium == 1) {
        // FIXME apply this itemwise
        critval = std::sqrt(arma::accu(arma::square(ipars - preM1)));
      }

      // FIXME tune and implement this
      //if(settings.err_tol > reltol) {
      //  settings.err_tol /= 2;
      //}
      //if(settings.iter_max < 100) {
      //  settings.iter_max += 5;
      //}

      if(critval < reltol * 10) {
        accelerator = 0;
      }

      if((critval < reltol) | iter == maxit) {
        run = false;
      }
      timer.step("postEM");
    }
  } else {
    ll_new = arma::accu(ll_g);
  }

  timer.step("stop");
  Rcpp::NumericVector time(timer);
  Rcpp::List ret;
  ret["ipars"] = Rcpp::wrap(ipars);
  ret["mu"] = Rcpp::wrap(mu);
  ret["sg"] = Rcpp::wrap(sg);
  ret["iter"] = Rcpp::wrap(iter);
  ret["ll"] = Rcpp::wrap(ll_new);
  ret["y_"] = Rcpp::wrap(y_);
  ret["rgl"] = Rcpp::wrap(rgl);
  ret["time"] = Rcpp::wrap(time);

  return ret;

  END_RCPP
}

