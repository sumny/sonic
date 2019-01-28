#include "optim/optim.hpp"
#include "var.hpp"

// FIXME fix pointer syntax stuff
// FIXME const

// Mstep items
bool sonic::Mstep_items(arma::vec& ipars, const bool global, const arma::uword optimizer, const arma::uword N, const arma::cube& rj_g, const arma::uvec& itemopt, sonic::Mstep_data_g *opt_data_g, sonic::Mstep_data_g_wh *opt_data_g_wh, sonic::Mstep_data_i *opt_data_i, sonic::Mstep_data_i_wh *opt_data_i_wh, optim::algo_settings_t& settings)
{
  bool converged;
  if(global) {
    if(optimizer == 0) {
      (*opt_data_g_wh).rj_g = std::move(rj_g);
      converged = optim::newton(ipars, sonic::llfun_g_wh, opt_data_g_wh, settings);
    } else if(optimizer == 1) {
      (*opt_data_g).rj_g = std::move(rj_g);
      converged = optim::bfgs(ipars, sonic::llfun_g, opt_data_g, settings);
    } else if(optimizer == 2) {
      (*opt_data_g).rj_g = std::move(rj_g);
      converged = optim::lbfgs(ipars, sonic::llfun_g, opt_data_g, settings);
    } else if(optimizer == 3) {
      (*opt_data_g).rj_g = std::move(rj_g);
      converged = optim::cg(ipars, sonic::llfun_g, opt_data_g, settings);
    } else if(optimizer == 4) {
      (*opt_data_g).rj_g = std::move(rj_g);
      converged = optim::gd(ipars, sonic::llfun_g, opt_data_g, settings);
    } else if(optimizer == 5) {
      (*opt_data_g).rj_g = std::move(rj_g);
      converged = optim::nm(ipars, sonic::llfun_g, opt_data_g, settings);
    } else if(optimizer == 6) {
      (*opt_data_g).rj_g = std::move(rj_g);
      converged = optim::de(ipars, sonic::llfun_g, opt_data_g, settings);
    } else if(optimizer == 7) {
      (*opt_data_g).rj_g = std::move(rj_g);
      converged = optim::pso(ipars, sonic::llfun_g, opt_data_g, settings);
    }
  } else {
    arma::vec ipars_i(2, arma::fill::none);
    for(arma::uword j = 0; j < N; ++j) {
      if(itemopt(j) == 1) {
        ipars_i(0) = ipars(j);
        ipars_i(1) = ipars(N + j);
        if(optimizer == 0) {
          (*opt_data_i_wh).rj_g = std::move(rj_g.cols(2 * j, 2 * j + 1));
          converged = optim::newton(ipars_i, sonic::llfun_i_wh, opt_data_i_wh, settings);
        } else if(optimizer == 1) {
          (*opt_data_i).rj_g = std::move(rj_g.cols(2 * j, 2 * j + 1));
          converged = optim::bfgs(ipars_i, sonic::llfun_i, opt_data_i, settings);
        } else if(optimizer == 2) {
          (*opt_data_i).rj_g = std::move(rj_g.cols(2 * j, 2 * j + 1));
          converged = optim::lbfgs(ipars_i, sonic::llfun_i, opt_data_i, settings);
        } else if(optimizer == 3) {
          (*opt_data_i).rj_g = std::move(rj_g.cols(2 * j, 2 * j + 1));
          converged = optim::cg(ipars_i, sonic::llfun_i, opt_data_i, settings);
        } else if(optimizer == 4) {
          (*opt_data_i).rj_g = std::move(rj_g.cols(2 * j, 2 * j + 1));
          converged = optim::gd(ipars_i, sonic::llfun_i, opt_data_i, settings);
        } else if(optimizer == 5) {
          (*opt_data_i).rj_g = std::move(rj_g.cols(2 * j, 2 * j + 1));
          converged = optim::nm(ipars_i, sonic::llfun_i, opt_data_i, settings);
        } else if(optimizer == 6) {
          (*opt_data_i).rj_g = std::move(rj_g.cols(2 * j, 2 * j + 1));
          converged = optim::de(ipars_i, sonic::llfun_i, opt_data_i, settings);
        } else if(optimizer == 7) {
          (*opt_data_i).rj_g = std::move(rj_g.cols(2 * j, 2 * j + 1));
          converged = optim::pso(ipars_i, sonic::llfun_i, opt_data_i, settings);
        }
        ipars(j) = ipars_i(0);
        ipars(N + j) = ipars_i(1);
      }
    }
  }
  return(converged);
}



// Mstep groups
void sonic::Mstep_groups(const arma::uword G, const arma::uword P, const arma::uword Q, arma::vec& mu, arma::vec& sg, arma::mat& AX, const arma::vec& X, const arma::mat& pul, const arma::mat& rgl)
{
  arma::mat grp_tmp = pul;
  arma::mat grp_X(Q, P, arma::fill::none);
  // FIXME easier
  grp_X.each_col() = X;
  arma::mat grp_rgl(Q, P, arma::fill::none);

  for(arma::uword g = 1; g < G; ++g) {
    grp_tmp = pul;
    grp_tmp.each_col() %= AX.col(g);
    grp_tmp.each_row() /= arma::sum(grp_tmp, 0);
    grp_rgl.each_row() = rgl.col(g);
    mu(g) = (1 / arma::accu(rgl.col(g))) * arma::accu(grp_X % grp_tmp % grp_rgl);
    sg(g) = std::sqrt((1 / arma::accu(rgl.col(g))) * arma::accu(arma::square(grp_X - mu(g)) % grp_tmp % grp_rgl));
    AX.col(g) = arma::normpdf(X, mu(g), sg(g));
    AX.col(g) /= arma::accu(AX.col(g));
  }
}



// log-likelihood function itemwise with hessian
double sonic::llfun_i_wh(const arma::vec &pars_in, arma::vec *grad_out, arma::mat *hess_out, void *opt_data)
{
  // structure Mstep_data_i_wh
  sonic::Mstep_data_i_wh* objfn_data = reinterpret_cast<sonic::Mstep_data_i_wh*>(opt_data);

  // probs = (log) probability solving (1) and not solving (0) item j given quadratures (Q x 2)
  (objfn_data->probs).col(1) = arma::pow(1 + arma::trunc_exp(- (objfn_data->X * pars_in(0) + pars_in(1))), -1);
  (objfn_data->probs).col(0) = 1 - (objfn_data->probs).col(1);

  if(grad_out) {
    (objfn_data->gr_out).zeros();
    for(arma::uword g = 0; g < objfn_data->G; ++g) {
      objfn_data->gr_tmp = ((objfn_data->rj_g).slice(g)).col(1) - (arma::sum((objfn_data->rj_g).slice(g), 1) % (objfn_data->probs).col(1));
      (objfn_data->gr_out)(0) += arma::accu(objfn_data->gr_tmp % objfn_data->X);
      (objfn_data->gr_out)(1) += arma::accu(objfn_data->gr_tmp);
    }
  *grad_out = - objfn_data->gr_out;
  }

  if(hess_out) {
    (objfn_data->hess_out).zeros();
    objfn_data->W_tmp = (objfn_data->probs).col(1) % (objfn_data->probs).col(0);
    for(arma::uword g = 0; g < objfn_data->G; ++g) {
      (objfn_data->hess_out)(0, 0) += arma::accu(arma::sum((objfn_data->rj_g).slice(g), 1) % objfn_data->W_tmp % arma::square(objfn_data->X));
      (objfn_data->hess_out)(0, 1) += arma::accu(arma::sum((objfn_data->rj_g).slice(g), 1) % objfn_data->W_tmp % objfn_data->X);
      (objfn_data->hess_out)(1, 1) += arma::accu(arma::sum((objfn_data->rj_g).slice(g), 1) % objfn_data->W_tmp);
    }
  (objfn_data->hess_out)(1, 0) += (objfn_data->hess_out)(0, 1);
  *hess_out = objfn_data->hess_out;
  }

  objfn_data->probs = arma::trunc_log(objfn_data->probs);
  objfn_data->ll = arma::accu((objfn_data->rj_g).slice(0) % objfn_data->probs);
  for(arma::uword g = 1; g < objfn_data->G; ++g) {
    objfn_data->ll += arma::accu((objfn_data->rj_g).slice(g) % objfn_data->probs);
  }

  return(- objfn_data->ll);
}



// log-likelihood function itemwise
double sonic::llfun_i(const arma::vec &pars_in, arma::vec *grad_out, void *opt_data)
{
  // structure Mstep_data_i
  sonic::Mstep_data_i* objfn_data = reinterpret_cast<sonic::Mstep_data_i*>(opt_data);

  // probs = (log) probability solving (1) and not solving (0) item j given quadratures (Q x 2)
  (objfn_data->probs).col(1) = arma::pow(1 + arma::trunc_exp(- (objfn_data->X * pars_in(0) + pars_in(1))), -1);
  (objfn_data->probs).col(0) = 1 - (objfn_data->probs).col(1);

  if(grad_out) {
    (objfn_data->gr_out).zeros();
    for(arma::uword g = 0; g < objfn_data->G; ++g) {
      objfn_data->gr_tmp = ((objfn_data->rj_g).slice(g)).col(1) - (arma::sum((objfn_data->rj_g).slice(g), 1) % (objfn_data->probs).col(1));
      (objfn_data->gr_out)(0) += arma::accu(objfn_data->gr_tmp % objfn_data->X);
      (objfn_data->gr_out)(1) += arma::accu(objfn_data->gr_tmp);
    }
  *grad_out = - objfn_data->gr_out;
  }

  objfn_data->probs = arma::trunc_log(objfn_data->probs);
  objfn_data->ll = arma::accu((objfn_data->rj_g).slice(0) % objfn_data->probs);
  for(arma::uword g = 1; g < objfn_data->G; ++g) {
    objfn_data->ll += arma::accu((objfn_data->rj_g).slice(g) % objfn_data->probs);
  }

  return(- objfn_data->ll);
}



// log-likelihood function global with hessian
double sonic::llfun_g_wh(const arma::vec &pars_in, arma::vec *grad_out, arma::mat *hess_out, void *opt_data)
{
  // structure Mstep_data_g
  sonic::Mstep_data_g_wh* objfn_data = reinterpret_cast<sonic::Mstep_data_g_wh*>(opt_data);

  // FIXME easier
  objfn_data->Prob = objfn_data->X * pars_in(objfn_data->a_ind).t();
  (objfn_data->Prob).each_row() += pars_in(objfn_data->d_ind).t();
  objfn_data->Prob = arma::pow(1 + arma::trunc_exp(- (objfn_data->Prob)), -1);
  (objfn_data->Probs).cols(objfn_data->z_ind) = 1 - objfn_data->Prob;
  (objfn_data->Probs).cols(objfn_data->o_ind) = objfn_data->Prob;

  if(grad_out) {
    (objfn_data->gr_out).zeros();
    for(arma::uword g = 0; g < objfn_data->G; ++g) {
      // FIXME easier
      objfn_data->gr_tmp = ((objfn_data->rj_g).slice(g)).cols(objfn_data->o_ind) - ((((objfn_data->rj_g).slice(g)).cols(objfn_data->z_ind) + ((objfn_data->rj_g).slice(g)).cols(objfn_data->o_ind)) % (objfn_data->Probs).cols(objfn_data->o_ind));
      (objfn_data->gr_out)(objfn_data->d_ind) += arma::sum((objfn_data->gr_tmp), 0);
      (objfn_data->gr_tmp).each_col() %= objfn_data->X;
      (objfn_data->gr_out)(objfn_data->a_ind) += arma::sum((objfn_data->gr_tmp), 0);
    }
    *grad_out = - objfn_data->gr_out;
  }

  if(hess_out) {
    (objfn_data->hess_out).zeros();
    objfn_data->W_tmp = (objfn_data->Probs).cols(objfn_data->o_ind) % (objfn_data->Probs).cols(objfn_data->o_ind);
    for(arma::uword g = 0; g < objfn_data->G; ++g) {
      objfn_data->n_vec.for_each( [&objfn_data, &g](const arma::uword &j) {
        (objfn_data->hess_out)(2 * j, 2 * j) += arma::accu(arma::sum((objfn_data->rj_g).slice(g), 2 * j + 1) % objfn_data->W_tmp % arma::square(objfn_data->X));
        (objfn_data->hess_out)(2 * j, 2 * j + 1) += arma::accu(arma::sum((objfn_data->rj_g).slice(g), 2 * j + 1) % objfn_data->W_tmp % objfn_data->X);
        (objfn_data->hess_out)(2 * j + 1, 2 * j) += (objfn_data->hess_out)(2 * j, 2 * j + 1);
        (objfn_data->hess_out)(2 * j + 1, 2 * j + 1) += arma::accu(arma::sum((objfn_data->rj_g).slice(g), 2 * j + 1) % objfn_data->W_tmp);
      });
    }
  *hess_out = objfn_data->hess_out;
  }

  objfn_data->Probs = arma::trunc_log(objfn_data->Probs);
  objfn_data->ll = arma::accu((objfn_data->rj_g).slice(0) % objfn_data->Probs);
  for(arma::uword g = 1; g < objfn_data->G; ++g) {
     objfn_data->ll += arma::accu((objfn_data->rj_g).slice(g) % objfn_data->Probs);
  }
  
  return(- objfn_data->ll);
}



// log-likelihood function global
double sonic::llfun_g(const arma::vec &pars_in, arma::vec *grad_out, void *opt_data)
{
  // structure Mstep_data_g
  sonic::Mstep_data_g* objfn_data = reinterpret_cast<sonic::Mstep_data_g*>(opt_data);

  // FIXME easier
  objfn_data->Prob = objfn_data->X * pars_in(objfn_data->a_ind).t();
  (objfn_data->Prob).each_row() += pars_in(objfn_data->d_ind).t();
  objfn_data->Prob = arma::pow(1 + arma::trunc_exp(- (objfn_data->Prob)), -1);
  (objfn_data->Probs).cols(objfn_data->z_ind) = 1 - objfn_data->Prob;
  (objfn_data->Probs).cols(objfn_data->o_ind) = objfn_data->Prob;

  if(grad_out) {
    (objfn_data->gr_out).zeros();
    for(arma::uword g = 0; g < objfn_data->G; ++g) {
      // FIXME easier
      objfn_data->gr_tmp = ((objfn_data->rj_g).slice(g)).cols(objfn_data->o_ind) - ((((objfn_data->rj_g).slice(g)).cols(objfn_data->z_ind) + ((objfn_data->rj_g).slice(g)).cols(objfn_data->o_ind)) % (objfn_data->Probs).cols(objfn_data->o_ind));
      (objfn_data->gr_out)(objfn_data->d_ind) += arma::sum((objfn_data->gr_tmp), 0);
      (objfn_data->gr_tmp).each_col() %= objfn_data->X;
      (objfn_data->gr_out)(objfn_data->a_ind) += arma::sum((objfn_data->gr_tmp), 0);
    }
    *grad_out = - objfn_data->gr_out;
  }

  objfn_data->Probs = arma::trunc_log(objfn_data->Probs);
  objfn_data->ll = arma::accu((objfn_data->rj_g).slice(0) % objfn_data->Probs);
  for(arma::uword g = 1; g < objfn_data->G; ++g) {
     objfn_data->ll += arma::accu((objfn_data->rj_g).slice(g) % objfn_data->Probs);
  }
  
  return(- objfn_data->ll);
}

