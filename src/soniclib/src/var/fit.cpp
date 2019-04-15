#include <cmath>
#include <math.h>
//#include <omp.h>
#include "optim/optim.hpp"
#include <Rcpp/Benchmark/Timer.h>
#include <RcppArmadillo.h>
//#include <time.h>
#include <utility>
#include "var/var.hpp"

// fit function
// FIXME check order of arguments in all functions etc. and rearrange if needed
// FIXME faster way for G = 1
// FIXME reduce maxit = 0 time
// FIXME fully implement OPENMP; currently disabled
// FIXME check which par vectors are needed
// FIXME infinite values in optimization
// FIXME NA values
// FIXME probs functions in misc.cpp
// FIXME always else if instead else
RcppExport SEXP fit(SEXP Ry, SEXP Rweights, SEXP Rimpact, SEXP Rstart, SEXP Rmodel, SEXP Rcontrol, SEXP Ralgo_settings)
{
  BEGIN_RCPP

  Rcpp::Timer timer;
  timer.step("start");

  // initial (const) declarations
  // FIXME try to use one-liners
  const arma::mat y = Rcpp::as<arma::mat>(Ry);
  const arma::vec weights = Rcpp::as<arma::vec>(Rweights);
  const arma::uvec impact = Rcpp::as<arma::uvec>(Rimpact);
  const arma::uvec groups = arma::unique(impact);
  const Rcpp::List start = Rcpp::as<Rcpp::List>(Rstart);
  const arma::uword model = Rcpp::as<arma::uword>(Rmodel); // 0 == "2PL", 1 == "RM", 2 == "3PL", 3 == "3PLu", 4 == "4PL"
  const Rcpp::List control = Rcpp::as<Rcpp::List>(Rcontrol);
  const Rcpp::List algo_settings = Rcpp::as<Rcpp::List>(Ralgo_settings);
  const arma::uword M = y.n_rows, N = y.n_cols, G = groups.n_elem;
  const arma::uword optimizer = Rcpp::as<arma::uword>(control[0]), accelerator = Rcpp::as<arma::uword>(control[1]), maxit = Rcpp::as<arma::uword>(control[2]), Q = Rcpp::as<arma::uword>(control[4]), criterion = Rcpp::as<arma::uword>(control[6]);
  const bool global = Rcpp::as<bool>(control[5]);
  const double reltol = Rcpp::as<double>(control[3]);
  const bool compute_ll = (criterion == 0);
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
  arma::uword npars;
  arma::uvec a_ind(N, arma::fill::none);
  arma::uvec d_ind(N, arma::fill::none);
  arma::uvec g_ind(N, arma::fill::none);
  arma::uvec u_ind(N, arma::fill::none);
  if(model == 0) {
    npars = 2;
    a_ind = arma::linspace<arma::uvec>(0, N - 1, N);
    d_ind = arma::linspace<arma::uvec>(N, (2 * N) - 1, N);
  } else if(model == 1) {
    npars = 1;
    d_ind = arma::linspace<arma::uvec>(0, N - 1, N);
  } else if(model == 2) {
    npars = 3;
    a_ind = arma::linspace<arma::uvec>(0, N - 1, N);
    d_ind = arma::linspace<arma::uvec>(N, (2 * N) - 1, N);
    g_ind = arma::linspace<arma::uvec>(2 * N, (3 * N) - 1, N);
  } else if(model == 3) {
    npars = 3;
    a_ind = arma::linspace<arma::uvec>(0, N - 1, N);
    d_ind = arma::linspace<arma::uvec>(N, (2 * N) - 1, N);
    u_ind = arma::linspace<arma::uvec>(2 * N, (3 * N) - 1, N);
  } else if(model == 4) {
    npars = 4;
    a_ind = arma::linspace<arma::uvec>(0, N - 1, N);
    d_ind = arma::linspace<arma::uvec>(N, (2 * N) - 1, N);
    g_ind = arma::linspace<arma::uvec>(2 * N, (3 * N) - 1, N);
    u_ind = arma::linspace<arma::uvec>(3 * N, (4 * N) - 1, N);
  }
  arma::uvec item_ind(npars, arma::fill::none);
  arma::mat pul(Q, P, arma::fill::none);
  arma::mat pgul(P, G, arma::fill::none);
  arma::cube rj_g(Q, 2 * N, G, arma::fill::none);
  double critval = reltol + 1, ll_new = 0, ll_old = 0;
  arma::vec itemnrm(N, arma::fill::none);
  arma::uvec itemopt = arma::ones<arma::uvec>(N);
  arma::uword iter = 0;
  bool run = true;
  Rcpp::LogicalVector convergence(maxit + 1, false);

  // Mstep stuff
  sonic::Mstep_data_g opt_data_g;
  sonic::Mstep_data_g_wh opt_data_g_wh;
  sonic::Mstep_data_i opt_data_i;
  sonic::Mstep_data_i_wh opt_data_i_wh;
  if(global) {
    if(optimizer == 0) {
      opt_data_g_wh.model = model;
      opt_data_g_wh.G = G;
      opt_data_g_wh.ll = 0;
      opt_data_g_wh.n_vec = n_vec;
      opt_data_g_wh.a_ind = a_ind;
      opt_data_g_wh.d_ind = d_ind;
      opt_data_g_wh.g_ind = g_ind;
      opt_data_g_wh.u_ind = u_ind;
      opt_data_g_wh.gr_out = arma::zeros<arma::vec>(npars * N);
      opt_data_g_wh.X = X;
      opt_data_g_wh.gr_tmp = arma::mat(Q, N, arma::fill::none);
      opt_data_g_wh.W_tmp = arma::mat(Q, N, arma::fill::none);
      opt_data_g_wh.Prob = arma::mat(Q, N, arma::fill::none);
      opt_data_g_wh.Probs = arma::mat(Q, 2 * N, arma::fill::none);
      opt_data_g_wh.z_ind = z_ind;
      opt_data_g_wh.o_ind = o_ind;
      opt_data_g_wh.hess_out = arma::zeros<arma::mat>(npars * N, npars * N);
    } else {
      opt_data_g.model = model;
      opt_data_g.G = G;
      opt_data_g.ll = 0;
      opt_data_g.a_ind = a_ind;
      opt_data_g.d_ind = d_ind;
      opt_data_g.g_ind = g_ind;
      opt_data_g.u_ind = u_ind;
      opt_data_g.gr_out = arma::zeros<arma::vec>(npars * N);
      opt_data_g.X = X;
      opt_data_g.gr_tmp = arma::mat(Q, N, arma::fill::none);
      opt_data_g.Prob = arma::mat(Q, N, arma::fill::none);
      opt_data_g.Probs = arma::mat(Q, 2 * N, arma::fill::none);
      opt_data_g.z_ind = z_ind;
      opt_data_g.o_ind = o_ind;
    }
  } else {
    if(optimizer == 0) {
      opt_data_i_wh.model = model;
      opt_data_i_wh.G = G;
      opt_data_i_wh.ll = 0;
      opt_data_i_wh.gr_out = arma::zeros<arma::vec>(npars);
      opt_data_i_wh.gr_tmp = arma::vec(Q, arma::fill::none);
      opt_data_i_wh.X = X;
      opt_data_i_wh.probs = arma::mat(Q, 2, arma::fill::none);
      opt_data_i_wh.W_tmp = arma::vec(Q, arma::fill::none);
      opt_data_i_wh.hess_out = arma::zeros<arma::mat>(npars, npars);
    } else {
      opt_data_i.model = model;
      opt_data_i.G = G;
      opt_data_i.ll = 0;
      opt_data_i.gr_out = arma::vec(npars, arma::fill::none);
      opt_data_i.gr_tmp = arma::vec(Q, arma::fill::none);
      opt_data_i.X = X;
      opt_data_i.probs = arma::mat(Q, 2, arma::fill::none);
    }
  }

  // optim algo settings
  optim::algo_settings_t settings;
  settings.err_tol = Rcpp::as<double>(algo_settings[0]);
  settings.iter_max = Rcpp::as<int>(algo_settings[1]);
  settings.lbfgs_par_M = Rcpp::as<int>(algo_settings[2]);
  settings.cg_method = Rcpp::as<int>(algo_settings[3]);
  settings.cg_restart_threshold = Rcpp::as<double>(algo_settings[4]);
  settings.gd_method = Rcpp::as<int>(algo_settings[5]);
  settings.nm_par_alpha = Rcpp::as<double>(algo_settings[6]);
  settings.nm_par_beta = Rcpp::as<double>(algo_settings[7]);
  settings.nm_par_gamma = Rcpp::as<double>(algo_settings[8]);
  settings.nm_par_delta = Rcpp::as<double>(algo_settings[9]);
  settings.gd_settings.step_size = Rcpp::as<double>(algo_settings[10]);
  settings.gd_settings.momentum_par = Rcpp::as<double>(algo_settings[11]);
  settings.gd_settings.norm_term = Rcpp::as<double>(algo_settings[12]);
  settings.gd_settings.ada_rho = Rcpp::as<double>(algo_settings[13]);
  settings.gd_settings.adam_beta_1 = Rcpp::as<double>(algo_settings[14]);
  settings.gd_settings.adam_beta_2 = Rcpp::as<double>(algo_settings[15]);
  settings.gd_settings.ada_max = Rcpp::as<bool>(algo_settings[16]);
  settings.de_n_pop = Rcpp::as<int>(algo_settings[17]);
  settings.de_n_gen = Rcpp::as<int>(algo_settings[18]);
  settings.de_check_freq = Rcpp::as<int>(algo_settings[19]);
  settings.de_mutation_method = Rcpp::as<int>(algo_settings[20]);
  settings.de_par_F = Rcpp::as<double>(algo_settings[21]);
  settings.de_par_CR = Rcpp::as<double>(algo_settings[22]);
  settings.pso_n_pop = Rcpp::as<int>( algo_settings[23]);
  settings.pso_n_gen = Rcpp::as<int>(algo_settings[24]);
  // FIXME check this
  //settings.pso_check_freq = algo_settings[25];
  settings.pso_inertia_method = Rcpp::as<int>(algo_settings[26]);
  settings.pso_par_w_min = Rcpp::as<double>(algo_settings[27]);
  settings.pso_par_w_max = Rcpp::as<double>(algo_settings[28]);
  settings.pso_par_w_damp = Rcpp::as<double>(algo_settings[29]);
  settings.pso_velocity_method = Rcpp::as<int>(algo_settings[30]);
  settings.pso_par_c_cog = Rcpp::as<double>(algo_settings[31]);
  settings.pso_par_c_soc = Rcpp::as<double>(algo_settings[32]);
  settings.pso_par_initial_c_cog = Rcpp::as<double>(algo_settings[33]);
  settings.pso_par_final_c_cog = Rcpp::as<double>(algo_settings[34]);
  settings.pso_par_initial_c_soc = Rcpp::as<double>(algo_settings[35]);
  settings.pso_par_final_c_soc = Rcpp::as<double>(algo_settings[36]);

  // accelerator stuff
  arma::vec preM2 = ipars;
  arma::vec preM1 = preM2;
  arma::mat U(npars * N, 3, arma::fill::zeros);
  arma::mat V(npars * N, 3, arma::fill::zeros);

  // debug stuff
  //opt_data_g.rj_g = rj_g;
  //arma::vec gradient_g = arma::vec(npars * N);
  //double ll = sonic::llfun_g(ipars, &gradient_g, &opt_data_g);
  //run = false;

  timer.step("preE1");
  sonic::Estep(y_u, y_u_, rgl, G, ipars, X, AX, ll, pul, pgul, rj_g, p_vec, n_vec, a_ind, d_ind, g_ind, u_ind, true, model);
  timer.step("postE1");

  // run simulations only 60 seconds CPU time (see bottom of the loop)
  //std::clock_t extra_time = std::clock();

  // EM
  if(maxit > 0) {
    while(run) {
      timer.step("preM");
      Rcpp::checkUserInterrupt();
      iter += 1;
      ll_old = ll;
      convergence[iter - 1] = sonic::Mstep_items(model, npars, ipars, global, optimizer, N, rj_g, itemopt, &opt_data_g, &opt_data_g_wh, &opt_data_i, &opt_data_i_wh, settings);

      // group parameters
      if((model == 1) || (G > 1)) {
        sonic::Mstep_groups(model, G, P, Q, mu, sg, AX, X, pul, rgl);
      }
      timer.step("postM");

      // EM acceleration, 0 == "none", 1 == "Ramsay", 2 == "SQUAREM", 3 == "Zhou"
      timer.step("preA");
      sonic::accelerate(y_u, y_u_, rgl, N, G, P, Q, ipars, X, AX, p_vec, n_vec, a_ind, d_ind, g_ind, u_ind, preM1, preM2, accelerator, ll, U, V, iter, model);
      timer.step("postA");

      timer.step("preE");
      sonic::Estep(y_u, y_u_, rgl, G, ipars, X, AX, ll, pul, pgul, rj_g, p_vec, n_vec, a_ind, d_ind, g_ind, u_ind, compute_ll, model);
      ll_new = ll;
      timer.step("postE");

      // termination criterion, 0 == "ll", 1 == "l2", 2 == "l2_itemwise"
      // FIXME seperate function for termination
      if(criterion == 0) {
        critval = ll_new - ll_old;
      } else if(criterion == 1) {
        critval = arma::norm(ipars - preM1, 2);
      } else if(criterion == 2) {
        // FIXME check this
        n_vec.for_each( [&itemnrm, &ipars, &preM1, &item_ind, &itemopt, &reltol, &N, &model](const arma::uword &j) {
          item_ind(0) = j;
          if(model != 1) {
            item_ind(1) = j + 1;
          } else if((model == 2) || (model == 3) || (model == 4)) {
            item_ind(2) = j + 2;
          } else if(model == 4) {
            item_ind(3) = j + 3;
          }
          itemnrm(j) = arma::norm(ipars.elem(item_ind) - preM1.elem(item_ind), 2);
          if(itemnrm(j) <= (reltol / N)) {
            itemopt(j) = 0;
          }
        });
        if(arma::accu(itemopt) == 0) {
          critval = arma::accu(itemnrm);
        } else {
          critval = reltol + 1;
        }
      }

      // check for termination
      //if((std::abs(critval) <= reltol) || (iter == maxit) || (((std::clock() - extra_time) / (double) CLOCKS_PER_SEC) >= 60)) {
      if((std::abs(critval) <= reltol) || (iter == maxit)) {
        run = false;
        convergence = convergence[Rcpp::Range(0, iter - 1)];
        if(criterion != 0) {
          sonic::Estep(y_u, y_u_, rgl, G, ipars, X, AX, ll, pul, pgul, rj_g, p_vec, n_vec, a_ind, d_ind, g_ind, u_ind, true, model);
          ll_new = ll;
        }
      } else {
        preM2 = preM1;
        preM1 = ipars;
      }
    }
  } else {
    ll_new = ll;
    convergence = convergence[0];
  }

  double converged = ((std::abs(critval) <= reltol) && iter < maxit);

  timer.step("stop");
  Rcpp::NumericVector time(timer);

  Rcpp::List ret;
  ret["ipars"] = Rcpp::wrap(ipars);
  ret["mu"] = Rcpp::wrap(mu);
  ret["sg"] = Rcpp::wrap(sg);
  ret["ll"] = Rcpp::wrap(ll_new);
  ret["iter"] = Rcpp::wrap(iter);
  ret["converged"] = Rcpp::wrap(converged);
  ret["convergence"] = Rcpp::wrap(convergence);
  ret["critval"] = Rcpp::wrap(critval);
  ret["y_u"] = Rcpp::wrap(y_u);
  ret["rgl"] = Rcpp::wrap(rgl);
  ret["X"] = Rcpp::wrap(X);
  ret["AX"] = Rcpp::wrap(AX);
  ret["time"] = Rcpp::wrap(time);
  // debug stuff
  //ret["debug_ll"] = Rcpp::wrap(ll);
  //ret["debug_gr"] = Rcpp::wrap(gradient_g);
  //ret["debug_gr_tmp"] = Rcpp::wrap(opt_data_g.gr_tmp);

  return(ret);

  END_RCPP
}



// Register native routines
// tools::package_native_routine_registration_skeleton("sonic")
static const R_CallMethodDef CallEntries[] = {
  {"fit", (DL_FUNC) &fit, 7},
  {NULL, NULL, 0}
};



void R_init_sonic(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
