#include "var.hpp"

// accelerators
// FIXME look for redundant calculations
void sonic::accelerate(arma::mat *y_, arma::mat *rgl, int *G, arma::vec *ipars,
  arma::vec *X, arma::mat *AX, arma::vec *sqr_ll_g, arma::cube *sqr_rj_g,
  arma::mat *Prob, arma:: mat *Probs, arma::mat *pul, arma::mat *pgul,
  arma::uvec *p_vec, arma::uvec *n_vec, arma::uvec *a_ind, arma::uvec *d_ind,
  arma::umat *cat_ind,
  arma::vec *preM1, arma::vec *preM2, arma::vec *accels, double *val,
  int *sqr_iter, int *accelerator, arma::vec *ll_g)
{
  if(*accelerator == 1) {
    // Ramsay
    *val = 1 - (std::sqrt(arma::accu(arma::square(*ipars - *preM1)) / arma::accu(arma::square(*ipars - (2 * *preM1) + *preM2))));
    if(*val < -5) {
      *val = -5;
    }
    *accels = ((1 - *val) * *ipars) + (*val * *preM1);
  } else if(*accelerator == 2) {
    // SQUAREM
    *val = - std::sqrt(arma::accu(arma::square(*preM1 - *preM2)) / arma::accu(arma::square(*ipars - (2 * *preM1) + *preM2)));
    if(*val > -1) {
      *val = -1;
    } else {
      *sqr_iter = 0;
      while(true) {
        *accels = *preM2 - ((2 * *val) * (*preM1 - *preM2)) + ((*val * *val) * (*ipars - (2 * *preM1) + *preM2));
        sonic::Estep(y_, rgl, G, accels, X, AX, sqr_ll_g, sqr_rj_g, Prob, Probs, pul, pgul, p_vec, n_vec, a_ind, d_ind, cat_ind);
        if(arma::accu(*sqr_ll_g) <= arma::accu(*ll_g)) {
          *val <- (*val - 1) / 2;
          *sqr_iter += 1;
        } else {
          break;
        }
        if(*sqr_iter == 4) {
          *val = -1;
          break;
        }
      }
    }
    *accels = *preM2 - ((2 * *val) * (*preM1 - *preM2)) + ((*val * *val) * (*ipars - (2 * *preM1) + *preM2));
  }
  *ipars = *accels;
}

