#ifndef _Estep_HPP
#define _Estep_HPP

// Estep
void Estep(const arma::mat& y_u, const arma::mat& y_u_, const arma::mat& rgl, const arma::uword G, const arma::vec& ipars, const arma::vec& X, const arma::mat& AX, double& ll, arma::mat& pul, arma::mat& pgul, arma::cube& rj_g, const arma::uvec& p_vec, const arma::uvec& n_vec, const arma::uvec& a_ind, const arma::uvec& d_ind, const bool compute_ll);

#endif
