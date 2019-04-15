#include "optim/optim.hpp"
#include "var.hpp"

// accelerator
void sonic::accelerate(const arma::mat& y_u, const arma::mat& y_u_, const arma::mat& rgl, const arma::uword N, const arma::uword G, const arma::uword P, const arma::uword Q, arma::vec& ipars, const arma::vec& X, const arma::mat& AX, const arma::uvec& p_vec, const arma::uvec& n_vec, const arma::uvec& a_ind, const arma::uvec& d_ind, const arma::uvec& g_ind, const arma::uvec& u_ind, const arma::vec& preM1, const arma::vec& preM2, const arma::uword accelerator, const double ll, arma::mat& U, arma::mat& V, const arma::uword iter, const arma::uword model)
{
  arma::vec accels = ipars;
  double val = 0, tmp_ll = ll;
  arma::uword tmp_iter = 0; // start counting at 0
  arma::mat tmp_pul(Q, P, arma::fill::none);
  arma::mat tmp_pgul(P, G, arma::fill::none);
  arma::cube tmp_rj_g(Q, 2 * N, G, arma::fill::none);
  
  if(accelerator == 1) {
    // Ramsay
    if(std::fmod(iter, 3) == 0) {
      //val = 1 - (std::sqrt(arma::accu(arma::square(ipars - preM1)) / arma::accu(arma::square(ipars - (2 * preM1) + preM2))));
      val = 1 - (arma::norm(ipars - preM1, 2) / arma::norm(ipars - (2 * preM1) + preM2, 2));
      if(val < -5) {
        val = -5;
      }
      accels = ((1 - val) * ipars) + (val * preM1);
   }
  } else if(accelerator == 2) {
    // SQUAREM SqS3
    if(std::fmod(iter, 3) == 0) {
      //val = - std::sqrt(arma::accu(arma::square(preM1 - preM2)) / arma::accu(arma::square(ipars - (2 * preM1) + preM2)));
      val = - (arma::norm(preM1 - preM2, 2) / arma::norm(ipars - (2 * preM1) + preM2, 2));
      if(val > -1) {
        val = -1;
      } else {
        tmp_iter = 0;
        while((tmp_ll <= ll) && (tmp_iter < 5)) {
          accels = preM2 - ((2 * val) * (preM1 - preM2)) + ((val * val) * (ipars - (2 * preM1) + preM2));
          sonic::Estep(y_u, y_u_, rgl, G, accels, X, AX, tmp_ll, tmp_pul, tmp_pgul, tmp_rj_g, p_vec, n_vec, a_ind, d_ind, g_ind, u_ind, true, model);
          if(tmp_ll <= ll) {
            val = (val - 1) / 2;
            tmp_iter += 1;
          }
        }
      }
      accels = preM2 - ((2 * val) * (preM1 - preM2)) + ((val * val) * (ipars - (2 * preM1) + preM2));
    }
  } else if(accelerator == 3) {
    // Zhou
    if(std::fmod(iter, 5) == 0) {
      accels = preM1 - V * arma::inv((U.t() * U - U.t() * V)) * U.t() * (preM2 - preM1);
      sonic::Estep(y_u, y_u_, rgl, G, accels, X, AX, tmp_ll, tmp_pul, tmp_pgul, tmp_rj_g, p_vec, n_vec, a_ind, d_ind, g_ind, u_ind, true, model);
      if(tmp_ll <= ll) {
        accels = ipars;
      }
    }
    U.col(0) = U.col(1);
    U.col(1) = U.col(2);
    U.col(2) = preM1 - preM2;
    V.col(0) = V.col(1);
    V.col(1) = V.col(2);
    V.col(2) = ipars - preM1;
  }

  ipars = accels;
}
