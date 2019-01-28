#ifndef _accelerate_HPP
#define _accelerate_HPP

#include "optim/optim.hpp"
#include "var.hpp"

// accelerators
void accelerate(const arma::mat& y_u, const arma::mat& y_u_, const arma::mat& rgl, const arma::uword N, const arma::uword G, const arma::uword P, const arma::uword Q, arma::vec& ipars, const arma::vec& X, const arma::mat& AX, const arma::uvec& p_vec, const arma::uvec& n_vec, const arma::uvec& a_ind, const arma::uvec& d_ind, const arma::vec& preM1, const arma::vec& preM2, const arma::uword accelerator, const double ll, arma::vec& preM2_va, arma::vec& preM1_va, arma::vec& ipars_va, const bool global, const arma::uword optimizer, const arma::cube& rj_g, sonic::Mstep_data_g *opt_data_g, sonic::Mstep_data_g_wh *opt_data_g_wh, sonic::Mstep_data_i *opt_data_i, sonic::Mstep_data_i_wh *opt_data_i_wh, optim::algo_settings_t& settings);

#endif

