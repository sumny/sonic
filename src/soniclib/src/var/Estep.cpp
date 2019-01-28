#include "var.hpp"

// Estep
void sonic::Estep(const arma::mat& y_u, const arma::mat& y_u_, const arma::mat& rgl, const arma::uword G, const arma::vec& ipars, const arma::vec& X, const arma::mat& AX, double& ll, arma::mat& pul, arma::mat& pgul, arma::cube& rj_g, const arma::uvec& p_vec, const arma::uvec& n_vec, const arma::uvec& a_ind, const arma::uvec& d_ind)
{
  // Prob and Prob_ = log probability solving and not solving each item given quadratures (Q x N)
  arma::mat Prob = - X * (ipars(a_ind)).t();
  Prob.each_row() -= (ipars(d_ind)).t();
  Prob = 1 / (1 + arma::trunc_exp(Prob));
  arma::mat Prob_ = arma::trunc_log(1 - Prob);
  Prob = arma::trunc_log(Prob);

  // pul = probability of patterns given quadratures (Q x P)
  p_vec.for_each( [&pul, &Prob, &Prob_, &y_u, &y_u_](const arma::uword &l) {
    pul.col(l) = arma::trunc_exp((Prob * (y_u.row(l)).t()) + (Prob_ * (y_u_.row(l)).t()));
  });

  // pgul = pul weighted with groupwise quadrature weights (P x G)
  // FIXME drop this into the loop below
  pgul = pul.t() * AX;

  // rj_g = expected number of persons solving (1, 3, ...) and not solving (0, 2, ...) each item given quadratures sliced for groups (Q, 2 * N, G)
  for(arma::uword g = 0; g < G; ++g) {
    (rj_g.slice(g)).each_col(2 * n_vec) = (pul * ((1 / pgul.col(g)) % rgl.col(g))) % AX.col(g);
    n_vec.for_each( [&rj_g, &g, &pul, &pgul, &rgl, &y_u, &AX](const arma::uword &j) {
      (rj_g.slice(g)).col((2 * j) + 1) = (pul * ((1 / pgul.col(g)) % rgl.col(g) % y_u.col(j))) % AX.col(g);
      (rj_g.slice(g)).col(2 * j) -= (rj_g.slice(g)).col((2 * j) + 1);
    });
  }

  // ll
  ll = arma::accu(rgl.col(0) % arma::trunc_log(pgul.col(0)));
  for(arma::uword g = 1; g < G; ++g) {
    ll += arma::accu(rgl.col(g) % arma::trunc_log(pgul.col(g)));
  }
}

