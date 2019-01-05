#include "var.hpp"

// preprocess data
void sonic::preprocess(arma::mat *y, arma::vec *weights,
  arma::vec *impact, int *G, arma::vec *primes, int *P, arma::mat *rgl,
  arma::mat *y_)
{
  // arma::unique does not preserve the original ordering
  // instead, it sorts in ascending order
  // therefore, reordered data structures are being used internally

  // unique patterns
  arma::vec pats = (*y + 1) * arma::log(*primes);
  arma::uvec pid_unique = arma::find_unique(pats);
  arma::vec pats_unique = pats(pid_unique);
  *P = pid_unique.size();

  // unique patterns groupwise and rgl
  // FIXME redundancy
  arma::uvec p_vec = arma::linspace<arma::uvec>(0, *P - 1, *P);
  arma::uvec g_id;
  arma::mat y_g;
  arma::vec pats_g;
  arma::vec weights_g;
  arma::mat rgl_(*P, *G);
  for(int g = 0; g < *G; ++g) {
    g_id = arma::find(*impact == g);
    pats_g = pats(g_id);
    weights_g = (*weights)(g_id);
    p_vec.for_each( [&g, &weights_g, &pats_g, &pats_unique, &rgl_](arma::uword &l) {
      rgl_(l, g) = arma::accu(weights_g(arma::find(pats_g == pats_unique(l))));
    });
  }
  *rgl = rgl_;

  // unique data
  *y_ = (*y).rows(pid_unique);
}

