#ifndef _Mstep_HPP
#define _Mstep_HPP

#include "optim/optim.hpp"

void Mstep_groups(const arma::uword G, const arma::uword P, const arma::uword Q, arma::vec& mu, arma::vec& sg, arma::mat& AX, const arma::vec& X, const arma::mat& pul, const arma::mat& rgl);

// log-likelihood function data struct global
struct Mstep_data_g
{
  arma::uword G;
  double ll;
  arma::uvec a_ind;
  arma::uvec d_ind;
  arma::vec gr_out;
  arma::vec X;
  arma::mat gr_tmp;
  arma::mat Prob;
  arma::mat Probs;
  arma::uvec z_ind;
  arma::uvec o_ind;
  arma::cube rj_g;
};

// log-likelihood function global
double llfun_g(const arma::vec& pars_in, arma::vec *grad_out, void *opt_data);



// log-likelihood function data struct global with hessian
struct Mstep_data_g_wh
{
  arma::uword G;
  double ll;
  arma::uvec n_vec;
  arma::uvec a_ind;
  arma::uvec d_ind;
  arma::vec gr_out;
  arma::vec X;
  arma::mat gr_tmp;
  arma::mat W_tmp;
  arma::mat Prob;
  arma::mat Probs;
  arma::uvec z_ind;
  arma::uvec o_ind;
  arma::mat hess_out;
  arma::cube rj_g;
};

// log-likelihood function global with hessian
double llfun_g_wh(const arma::vec& pars_in, arma::vec *grad_out, arma::mat *hess_out, void *opt_data);



// log-likelihood function data struct itemwise
struct Mstep_data_i
{
  arma::uword G;
  double ll;
  arma::vec gr_out;
  arma::vec gr_tmp;
  arma::vec X;
  arma::mat probs;
  arma::cube rj_g;
};

// log-likelihood function itemwise
double llfun_i(const arma::vec& pars_in, arma::vec *grad_out, void *opt_data);



// log-likelihood function data struct itemwise with hessian
struct Mstep_data_i_wh
{
  arma::uword G;
  double ll;
  arma::vec gr_out;
  arma::vec gr_tmp;
  arma::vec X;
  arma::vec W_tmp;
  arma::mat probs;
  arma::mat hess_out;
  arma::cube rj_g;
};

// log-likelihood function itemwise with hessian
double llfun_i_wh(const arma::vec& pars_in, arma::vec *grad_out, arma::mat *hess_out, void *opt_data);



bool Mstep_items(arma::vec& ipars, const bool global, const arma::uword optimizer, const arma::uword N, const arma::cube& rj_g, const arma::uvec& itemopt, Mstep_data_g *opt_data_g, Mstep_data_g_wh *opt_data_g_wh, Mstep_data_i *opt_data_i, Mstep_data_i_wh *opt_data_i_wh, optim::algo_settings_t& settings);

#endif

