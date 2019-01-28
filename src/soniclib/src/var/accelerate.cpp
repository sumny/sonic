#include "optim/optim.hpp"
#include "var.hpp"

// accelerator
// FIXME check for redundancy
void sonic::accelerate(const arma::mat& y_u, const arma::mat& y_u_, const arma::mat& rgl, const arma::uword N, const arma::uword G, const arma::uword P, const arma::uword Q, arma::vec& ipars, const arma::vec& X, const arma::mat& AX, const arma::uvec& p_vec, const arma::uvec& n_vec, const arma::uvec& a_ind, const arma::uvec& d_ind, const arma::vec& preM1, const arma::vec& preM2, const arma::uword accelerator, const double ll, arma::vec& preM2_va, arma::vec& preM1_va, arma::vec& ipars_va, const bool global, const arma::uword optimizer, const arma::cube& rj_g, sonic::Mstep_data_g *opt_data_g, sonic::Mstep_data_g_wh *opt_data_g_wh, sonic::Mstep_data_i *opt_data_i, sonic::Mstep_data_i_wh *opt_data_i_wh, optim::algo_settings_t& settings)

{
  arma::vec accels(2 * N, arma::fill::none);
  double val = 0, tmp_ll = 0;
  arma::uword sqr_iter = 0; // start counting at 0
  arma::mat tmp_pul(Q, P, arma::fill::none);
  arma::mat tmp_pgul(P, G, arma::fill::none);
  arma::cube tmp_rj_g(Q, 2 * N, G, arma::fill::none);
  const arma::uvec itemopt_va = arma::ones<arma::uvec>(N);

  if(accelerator == 1) {
    // Ramsay
    val = 1 - (std::sqrt(arma::accu(arma::square(ipars - preM1)) / arma::accu(arma::square(ipars - (2 * preM1) + preM2))));
    if(val < -5) {
      val = -5;
    }
    accels = ((1 - val) * ipars) + (val * preM1);
  } else if(accelerator == 2) {
    // SQUAREM SqS3
    val = - std::sqrt(arma::accu(arma::square(preM1 - preM2)) / arma::accu(arma::square(ipars - (2 * preM1) + preM2)));
    if(val > -1) {
      val = -1;
    } else {
      sqr_iter = 0;
      while((tmp_ll < ll) && (sqr_iter < 5)) {
        accels = preM2 - ((2 * val) * (preM1 - preM2)) + ((val * val) * (ipars - (2 * preM1) + preM2));
        sonic::Estep(y_u, y_u_, rgl, G, accels, X, AX, tmp_ll, tmp_pul, tmp_pgul, tmp_rj_g, p_vec, n_vec, a_ind, d_ind);
        if(tmp_ll < ll) {
          val = (val - 1) / 2;
          sqr_iter += 1;
        }
      }
    }
    accels = preM2 - ((2 * val) * (preM1 - preM2)) + ((val * val) * (ipars - (2 * preM1) + preM2));
  } else if(accelerator == 3) {
    // VA-Steffensen (Guo et al. 2017)
    // FIXME
    preM2_va = preM2 - ((std::sqrt(arma::accu(arma::square(preM1 - preM2))) * (preM2 - 2 * preM1 + ipars)) / std::sqrt(arma::accu(arma::square(preM2 - 2 * preM1 + ipars))));
    sonic::Estep(y_u, y_u_, rgl, G, preM2_va, X, AX, tmp_ll, tmp_pul, tmp_pgul, tmp_rj_g, p_vec, n_vec, a_ind, d_ind);
    preM1_va = preM2_va;
    bool converged = sonic::Mstep_items(preM1_va, global, optimizer, N, tmp_rj_g, itemopt_va, opt_data_g, opt_data_g_wh, opt_data_i, opt_data_i_wh, settings);
    sonic::Estep(y_u, y_u_, rgl, G, preM1_va, X, AX, tmp_ll, tmp_pul, tmp_pgul, tmp_rj_g, p_vec, n_vec, a_ind, d_ind);
    ipars_va = preM1_va;
    converged = sonic::Mstep_items(ipars_va, global, optimizer, N, tmp_rj_g, itemopt_va, opt_data_g, opt_data_g_wh, opt_data_i, opt_data_i_wh, settings);
    accels  = preM2_va - ((std::sqrt(arma::accu(arma::square(preM1_va - preM2_va))) * (preM2_va - 2 * preM1_va + ipars_va)) / std::sqrt(arma::accu(arma::square(preM2_va - 2 * preM1_va + ipars_va))));
  }
  ipars = accels;
}

