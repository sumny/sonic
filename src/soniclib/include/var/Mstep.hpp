#ifndef _Mstep_HPP
#define _Mstep_HPP

// log-likelihood function data struct global
struct Mstep_data_g
{
  int G;
  arma::vec ll_g;
  arma::uvec a_ind;
  arma::uvec d_ind;
  arma::vec gr_out;
  arma::vec X;
  arma::mat gr_a;
  arma::mat gr_d;
  arma::mat gr_tmp;
  arma::mat e;
  arma::mat Prob;
  arma::mat Probs;
  arma::umat cat_ind;
  arma::cube rj_g;
};



// log-likelihood function global
double llfun_g(const arma::vec &pars_in, arma::vec *grad_out, void *opt_data);



// log-likelihood function data struct itemwise
struct Mstep_data_i
{
  int j;
  int G;
  arma::vec gr_out;
  arma::vec gr_a;
  arma::vec gr_d;
  arma::vec gr_tmp;
  arma::vec ll_g;
  arma::vec X;
  arma::vec e;
  arma::mat p;
  arma::cube rj_g;
};



// log-likelihood function itemwise
double llfun_i(const arma::vec &pars_in, arma::vec *grad_out, void *opt_data);

#endif

