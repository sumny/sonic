#ifndef _accelerate_HPP
#define _accelerate_HPP

// accelerators
void accelerate(arma::mat *y_, arma::mat *rgl, int *G, arma::vec *ipars,
  arma::vec *X, arma::mat *AX, arma::vec *sqr_ll_g, arma::cube *sqr_rj_g,
  arma::mat *Prob, arma:: mat *Probs, arma::mat *pul, arma::mat *pgul,
  arma::uvec *p_vec, arma::uvec *n_vec, arma::uvec *a_ind, arma::uvec *d_ind,
  arma::umat *cat_ind,
  arma::vec *preM1, arma::vec *preM2, arma::vec *accels, double *val,
  int *sqr_iter, int *accelerator, arma::vec *ll_g);

#endif

