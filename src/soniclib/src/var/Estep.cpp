#include "var.hpp"

// Estep
// FIXME lambda functions where possible
// FIXME unnecessary transposing
void sonic::Estep(arma::mat *y_, arma::mat *rgl, int *G, arma::vec *ipars,
  arma::vec *X, arma::mat *AX, arma::vec *ll_g, arma::cube *rj_g,
  arma::mat *Prob, arma::mat *Probs, arma::mat *pul, arma::mat *pgul,
  arma::uvec *p_vec, arma::uvec *n_vec, arma::uvec *a_ind, arma::uvec *d_ind,
  arma::umat *cat_ind)
{
  // Probs = log probability solving (1, 3, ...) and not solving (0, 2, ...) items given quadratures (Q x 2 * N)
  *Prob = *X * ((*ipars)(*a_ind)).t();
  (*Prob).each_row() += ((*ipars)(*d_ind)).t();
  *Prob = arma::pow(1 + arma::trunc_exp(- *Prob), -1);
  (*Probs).cols((*cat_ind).col(0)) = 1 - *Prob;
  (*Probs).cols((*cat_ind).col(1)) = *Prob;
  *Probs = arma::trunc_log(*Probs);
 
  // pul = probability of patterns given quadratures (Q x P)
  // FIXME easier
  (*p_vec).for_each( [&pul, &Probs, &cat_ind, &y_](arma::uword &l) {
    (*pul).col(l) = arma::trunc_exp(((*Probs).cols((*cat_ind).col(1)) * ((*y_).row(l)).t()) + ((*Probs).cols((*cat_ind).col(0)) * (1 - (*y_).row(l)).t()));
  });

  // pgul = pul weighted with groupwise quadrature weights (P x G)
  *pgul = (*pul).t() * *AX;

  // rj_g = expected number of persons solving (1, 3, ...) and not solving (0, 2, ...) items given quadratures sliced for groups (Q, 2 * N, G)
  for(int g = 0; g < *G; ++g) {
    ((*rj_g).slice(g)).each_col((*cat_ind).col(0))= (*pul * (arma::pow((*pgul).col(g), -1) % (*rgl).col(g))) % (*AX).col(g);

    (*n_vec).for_each( [&g, &rj_g, &pul, &pgul, &rgl, &y_, &AX](arma::uword &j) {
      ((*rj_g).slice(g)).col((2 * j) + 1) = (*pul * (arma::pow((*pgul).col(g), -1) % (*rgl).col(g) % (*y_).col(j))) % (*AX).col(g);
      ((*rj_g).slice(g)).col(2 * j) -= ((*rj_g).slice(g)).col((2 * j) + 1);

    });
  }

  // ll_g = ll but only for group g (Q)
  for(int g = 0; g < *G; ++g) {
    (*ll_g)(g) = arma::accu((*rgl).col(g) % arma::trunc_log((*pgul).col(g)));
  }
}

