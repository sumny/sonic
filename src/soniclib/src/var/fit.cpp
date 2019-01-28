//#include <omp.h>
#include "optim/optim.hpp"
#include "var/var.hpp"
#include <math.h>
#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>
#include <utility>

// fit function
// FIXME replace arma::accu and arma::sum with matrix multiplication
// FIXME check order of arguments in all functions etc. and rearrange if needed
// FIXME faster way for G = 1
// FIXME reduce maxit = 0 time
// FIXME get OPENMP right; currently disabled
// FIXME check iteration sequences
// FIXME check which par vectors needed
RcppExport SEXP fit(SEXP Ry, SEXP Rweights, SEXP Rimpact, SEXP Rstart, SEXP Rcontrol, SEXP Ralgo_settings)
{
  BEGIN_RCPP

  Rcpp::Timer timer;
  timer.step("start");

  // initial const declarations
  // FIXME try to use oneliners
  const arma::mat y = Rcpp::as<arma::mat>(Ry);
  const arma::vec weights = Rcpp::as<arma::vec>(Rweights);
  const arma::uvec impact = Rcpp::as<arma::uvec>(Rimpact);
  const arma::uvec groups = arma::unique(impact);
  const Rcpp::List start = Rcpp::as<Rcpp::List>(Rstart);
  const Rcpp::List control = Rcpp::as<Rcpp::List>(Rcontrol);
  const Rcpp::List algo_settings = Rcpp::as<Rcpp::List>(Ralgo_settings);
  const arma::uword M = y.n_rows, N = y.n_cols, G = groups.n_elem;
  const arma::uword optimizer = Rcpp::as<arma::uword>(control[0]), accelerator = Rcpp::as<arma::uword>(control[1]), maxit = Rcpp::as<arma::uword>(control[2]), Q = Rcpp::as<arma::uword>(control[4]), criterium = Rcpp::as<arma::uword>(control[6]);
  const double reltol = Rcpp::as<double>(control[3]);
  const bool global = Rcpp::as<bool>(control[5]);
  timer.step("declarations");

  // unique patterns
  arma::uword P;
  arma::mat y_u;
  arma::mat rgl = sonic::preprocess(y, weights, impact, N, G, P, y_u);
  arma::mat y_u_ = 1 - y_u;
  timer.step("unique_data");

  // start values (slopes and intercepts, group means and variances)
  arma::vec ipars = Rcpp::as<arma::vec>(start[0]);
  arma::vec mu = Rcpp::as<arma::vec>(start[1]);
  arma::vec sg = arma::sqrt(Rcpp::as<arma::vec>(start[2]));

  // quadrature stuff
  // FIXME implement other types
  const arma::vec X = arma::linspace<arma::vec>(-6, 6, Q);
  arma::mat AX(Q, G, arma::fill::none);
  for(arma::uword g = 0; g < G; ++g) {
    AX.col(g) = arma::normpdf(X, mu(g), sg(g));
    AX.col(g) /= arma::accu(AX.col(g));
  }
  
  // some more declarations
  double ll = 0;
  const arma::uvec p_vec = arma::linspace<arma::uvec>(0, P - 1, P);
  const arma::uvec n_vec = arma::linspace<arma::uvec>(0, N - 1, N);
  const arma::uvec z_ind = 2 * n_vec;
  const arma::uvec o_ind = z_ind + 1;
  arma::uvec item_ind(2, arma::fill::none);
  arma::mat pul(Q, P, arma::fill::none);
  arma::mat pgul(P, G, arma::fill::none);
  arma::cube rj_g(Q, 2 * N, G, arma::fill::none);
  double critval = reltol + 1, ll_new = 0, ll_old = 0;
  arma::vec itemnrm(N, arma::fill::none);
  arma::uvec itemopt = arma::ones<arma::uvec>(N);
  arma::uword iter = 0;
  bool run = true;
  const arma::uvec a_ind = arma::linspace<arma::uvec>(0, N - 1, N);
  const arma::uvec d_ind = arma::linspace<arma::uvec>(N, 2 * N - 1, N);
  Rcpp::LogicalVector convergence(maxit, false);

  // Mstep stuff
  sonic::Mstep_data_g opt_data_g;
  sonic::Mstep_data_g_wh opt_data_g_wh;
  sonic::Mstep_data_i opt_data_i;
  sonic::Mstep_data_i_wh opt_data_i_wh;
  if(global) {
    if(optimizer == 0) {
      opt_data_g_wh.G = G;
      opt_data_g_wh.ll = 0;
      opt_data_g_wh.n_vec = n_vec;
      opt_data_g_wh.a_ind = a_ind;
      opt_data_g_wh.d_ind = d_ind;
      opt_data_g_wh.gr_out = arma::zeros<arma::vec>(2 * N);
      opt_data_g_wh.X = X;
      opt_data_g_wh.gr_tmp = arma::mat(Q, N, arma::fill::none);
      opt_data_g_wh.W_tmp = arma::mat(Q, N, arma::fill::none);
      opt_data_g_wh.Prob = arma::mat(Q, N, arma::fill::none);
      opt_data_g_wh.Probs = arma::mat(Q, 2 * N, arma::fill::none);
      opt_data_g_wh.z_ind = z_ind;
      opt_data_g_wh.o_ind = o_ind;
      opt_data_g_wh.hess_out = arma::zeros<arma::mat>(2 * N, 2 * N);
    } else {
      opt_data_g.G = G;
      opt_data_g.ll = 0;
      opt_data_g.a_ind = a_ind;
      opt_data_g.d_ind = d_ind;
      opt_data_g.gr_out = arma::zeros<arma::vec>(2 * N);
      opt_data_g.X = X;
      opt_data_g.gr_tmp = arma::mat(Q, N, arma::fill::none);
      opt_data_g.Prob = arma::mat(Q, N, arma::fill::none);
      opt_data_g.Probs = arma::mat(Q, 2 * N, arma::fill::none);
      opt_data_g.z_ind = z_ind;
      opt_data_g.o_ind = o_ind;
    }
  } else {
    if(optimizer == 0) {
      opt_data_i_wh.G = G;
      opt_data_i_wh.ll = 0;
      opt_data_i_wh.gr_out = arma::zeros<arma::vec>(2);
      opt_data_i_wh.gr_tmp = arma::vec(Q, arma::fill::none);
      opt_data_i_wh.X = X;
      opt_data_i_wh.probs = arma::mat(Q, 2, arma::fill::none);
      opt_data_i_wh.W_tmp = arma::vec(Q, arma::fill::none);
      opt_data_i_wh.hess_out = arma::zeros<arma::mat>(2, 2);
    } else {
      opt_data_i.G = G;
      opt_data_i.ll = 0;
      opt_data_i.gr_out = arma::vec(2, arma::fill::none);
      opt_data_i.gr_tmp = arma::vec(Q, arma::fill::none);
      opt_data_i.X = X;
      opt_data_i.probs = arma::mat(Q, 2, arma::fill::none);
    }
  }

  // optim algo settings
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
  settings.gd_settings.step_size = algo_settings[10];
  settings.gd_settings.momentum_par = algo_settings[11];
  settings.gd_settings.norm_term = algo_settings[12];
  settings.gd_settings.ada_rho = algo_settings[13];
  settings.gd_settings.adam_beta_1 = algo_settings[14];
  settings.gd_settings.adam_beta_2 = algo_settings[15];
  settings.gd_settings.ada_max = algo_settings[16];
  settings.de_n_pop = algo_settings[17];
  settings.de_n_gen = algo_settings[18];
  settings.de_check_freq = algo_settings[19];
  settings.de_mutation_method = algo_settings[20];
  settings.de_par_F = algo_settings[21];
  settings.de_par_CR = algo_settings[22];
  settings.pso_n_pop = algo_settings[23];
  settings.pso_n_gen = algo_settings[24];
  // FIXME check this
  //settings.de_check_freq = algo_settings[25];
  settings.pso_inertia_method = algo_settings[26];
  settings.pso_par_w_min = algo_settings[27];
  settings.pso_par_w_max = algo_settings[28];
  settings.pso_par_w_damp = algo_settings[29];
  settings.pso_velocity_method = algo_settings[30];
  settings.pso_par_c_cog = algo_settings[31];
  settings.pso_par_c_soc = algo_settings[32];
  settings.pso_par_initial_c_cog = algo_settings[33];
  settings.pso_par_final_c_cog = algo_settings[34];
  settings.pso_par_initial_c_soc = algo_settings[35];
  settings.pso_par_final_c_soc = algo_settings[36];

  timer.step("preE1");
  sonic::Estep(y_u, y_u_, rgl, G, ipars, X, AX, ll, pul, pgul, rj_g, p_vec, n_vec, a_ind, d_ind);
  timer.step("postE1");

  // accelerator stuff
  arma::vec preM2 = ipars;
  arma::vec preM1 = preM2;
  arma::vec preM2_va = ipars;
  arma::vec preM1_va = preM2_va;
  arma::vec ipars_va = preM1_va;

  // debug stuff
  //opt_data_g.rj_g = rj_g;
  //arma::vec gradient_g = arma::vec(2 * N);
  //double ll = sonic::llfun_g(ipars, &gradient_g, &opt_data_g);
  //run = false;
  
  // EM
  if(maxit > 0) {
    while(run) {
      timer.step("preM");
      Rcpp::checkUserInterrupt();
      iter += 1;
      ll_old = ll;
      convergence[iter - 1] = sonic::Mstep_items(ipars, global, optimizer, N, rj_g, itemopt, &opt_data_g, &opt_data_g_wh, &opt_data_i, &opt_data_i_wh, settings);

      // group parameters
      if(G > 1) {
        sonic::Mstep_groups(G, P, Q, mu, sg, AX, X, pul, rgl);
      }
      timer.step("postM");

      // acceleration
      // EM acceleration, 0 == "none", 1 == "Ramsay", 2 == "SQUAREM", FIXME 3 = "VA-Steffensen", ...
      // FIXME itemwise / parameterwise acceleration? Probably not needed
      timer.step("preA");
      if((accelerator != 0) && (std::fmod(iter, 2) == 0)) {
        sonic::accelerate(y_u, y_u_, rgl, N, G, P, Q, ipars, X, AX, p_vec, n_vec, a_ind, d_ind, preM1, preM2, accelerator, ll, preM2_va, preM1_va, ipars_va, global, optimizer, rj_g, &opt_data_g, &opt_data_g_wh, &opt_data_i, &opt_data_i_wh, settings);
      }
      timer.step("postA");

      timer.step("preE");
      sonic::Estep(y_u, y_u_, rgl, G, ipars, X, AX, ll, pul, pgul, rj_g, p_vec, n_vec, a_ind, d_ind);
      ll_new = ll;
      timer.step("postE");

      // termination criterium
      // termination criterium, 0 == "ll", 1 == "l2", 2 == "l2_itemwise", FIXME 3 = "ll_itemwise"
      // FIXME seperate function for termination
      if(accelerator == 3) {
        if(std::fmod(iter, 2) == 0) {
          critval = arma::accu(arma::square(ipars - preM2_va));
        }
      } else {
        if(criterium == 0) {
          critval = ll_new - ll_old;
        } else if(criterium == 1) {
          critval = std::sqrt(arma::accu(arma::square(ipars_va - preM1_va)));
        } else if(criterium == 2) {
          // FIXME check this
          n_vec.for_each( [&itemnrm, &ipars, &preM1, &item_ind, &itemopt, &reltol](const arma::uword &j) {
            item_ind(0) = j;
            item_ind(1) = j + 1;
            itemnrm(j) = std::sqrt(arma::accu(arma::square(ipars.elem(item_ind) - preM1.elem(item_ind))));
            if(itemnrm(j) < reltol) {
              itemopt(j) = 0;
            }
          });
          if(arma::accu(itemopt) == 0) {
            critval = 0;
          } else {
            critval = reltol + 1;
          }
        }
      }

      // check for termination
      if((critval < reltol) || (iter == maxit)) {
        run = false;
      } else {
        preM2 = preM1;
        preM1 = ipars;
      }

    }
  } else {
    ll_new = ll;
  }

  timer.step("stop");
  Rcpp::NumericVector time(timer);
  Rcpp::List ret;
  ret["ipars"] = Rcpp::wrap(ipars);
  ret["mu"] = Rcpp::wrap(mu);
  ret["sg"] = Rcpp::wrap(sg);
  ret["iter"] = Rcpp::wrap(iter);
  ret["ll"] = Rcpp::wrap(ll_new);
  ret["critval"] = Rcpp::wrap(critval);
  ret["y_u"] = Rcpp::wrap(y_u);
  ret["rgl"] = Rcpp::wrap(rgl);
  ret["time"] = Rcpp::wrap(time);
  // FIXME iter - 1 and what if maxit == 0
  //ret["convergence"] = Rcpp::wrap(convergence[Rcpp::Range(0, iter)]);
  //ret["itemnrm"] = Rcpp::wrap(itemnrm);
  //ret["preM2_va"] = Rcpp::wrap(preM2_va);
  //ret["preM1_va"] = Rcpp::wrap(preM1_va);
  //ret["ipars_va"] = Rcpp::wrap(ipars_va);
  // debug stuff
  //ret["debug_ll"] = Rcpp::wrap(ll);
  //ret["debug_gr"] = Rcpp::wrap(gradient_g);
  //ret["debug_gr_a"] = Rcpp::wrap(opt_data_g.gr_a);
  //ret["debug_gr_d"] = Rcpp::wrap(opt_data_g.gr_d);
  //ret["debug_gr_tmp"] = Rcpp::wrap(opt_data_g.gr_tmp);

  return(ret);

  END_RCPP
}



// Register native routines
// tools::package_native_routine_registration_skeleton("sonic")
static const R_CallMethodDef CallEntries[] = {
  {"fit", (DL_FUNC) &fit, 6},
  {NULL, NULL, 0}
};



void R_init_sonic(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

