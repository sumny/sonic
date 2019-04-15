#ifndef _accelerate_HPP
#define _accelerate_HPP

#include "optim/optim.hpp"
#include "var.hpp"

// accelerators
void accelerate(const arma::mat& y_u, const arma::mat& y_u_, const arma::mat& rgl, const arma::uword N, const arma::uword G, const arma::uword P, const arma::uword Q, arma::vec& ipars, const arma::vec& X, const arma::mat& AX, const arma::uvec& p_vec, const arma::uvec& n_vec, const arma::uvec& a_ind, const arma::uvec& d_ind, const arma::uvec& g_ind, const arma::uvec& u_ind, const arma::vec& preM1, const arma::vec& preM2, const arma::uword accelerator, const double ll, arma::mat& U, arma::mat& V, const arma::uword iter, const arma::uword model);

#endif
