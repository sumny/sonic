#ifndef _Estep_HPP
#define _Estep_HPP

// Estep
void Estep(arma::mat *y_, arma::mat *rgl, int *G, arma::vec *ipars,
  arma::vec *X, arma::mat *AX, arma::vec *ll_g, arma::cube *rj_g,
  arma::mat *Prob, arma::mat *Probs, arma::mat *pul, arma::mat *pgul,
  arma::uvec *p_vec, arma::uvec *n_vec, arma::uvec *a_ind, arma::uvec *d_ind,
  arma::umat *cat_ind);

#endif

