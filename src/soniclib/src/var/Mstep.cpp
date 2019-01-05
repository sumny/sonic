#include "var.hpp"

// FIXME fix pointer syntax stuff
// FIXME comment

// log-likelihood function itemwise
double sonic::llfun_i(const arma::vec &pars_in, arma::vec *grad_out, void *opt_data)
{
  // structure Mstep_data_i
  sonic::Mstep_data_i* objfn_data = reinterpret_cast<sonic::Mstep_data_i*>(opt_data);

  // e = exp(a * X + d) (Q)
  // p = log probability solving (1) and not solving (0) item j given quadratures (Q x 2)
  objfn_data->e = objfn_data->X * pars_in(0) + pars_in(1);
  (objfn_data->p).col(1) = arma::pow(1 + arma::trunc_exp(- objfn_data->e), -1);
  (objfn_data->p).col(0) = 1 - (objfn_data->p).col(1);
  objfn_data->p = arma::trunc_log(objfn_data->p);
  objfn_data->e = arma::trunc_exp(objfn_data->e);

  for(int g = 0; g < objfn_data->G; ++g) {
    //(objfn_data->ll_g)(g) = arma::accu(((objfn_data->rj_g).slice(g)).cols((objfn_data->cat_ind).row(objfn_data->j)) % objfn_data->p);
    (objfn_data->ll_g)(g) = arma::accu((objfn_data->rj_g).slice(g) % objfn_data->p);
  }

  if(grad_out) {
    for(int g = 0; g < objfn_data->G; ++g) {
      //objfn_data->gr_tmp = (((objfn_data->rj_g).slice(g)).col(2 * (objfn_data->j) + 1) - (((objfn_data->rj_g).slice(g)).col(2 * (objfn_data->j)) % objfn_data->e)) / (objfn_data->e + 1);
      objfn_data->gr_tmp = (((objfn_data->rj_g).slice(g)).col(1) - (((objfn_data->rj_g).slice(g)).col(0) % objfn_data->e)) / (objfn_data->e + 1);
      (objfn_data->gr_a)(g) = arma::accu(objfn_data->X % objfn_data->gr_tmp);
      (objfn_data->gr_d)(g) = arma::accu(objfn_data->gr_tmp);
    }
  (objfn_data->gr_out)(0) = arma::accu(objfn_data->gr_a);
  (objfn_data->gr_out)(1) = arma::accu(objfn_data->gr_d);
  *grad_out = - objfn_data->gr_out;
  }

  return(- arma::accu(objfn_data->ll_g));
}



// log-likelihood function global
double sonic::llfun_g(const arma::vec &pars_in, arma::vec *grad_out, void *opt_data)
{
  // FIXME transposed stuff etc.
  // structure Mstep_data_g
  sonic::Mstep_data_g* objfn_data = reinterpret_cast<sonic::Mstep_data_g*>(opt_data);

  objfn_data->e = objfn_data->X * pars_in(objfn_data->a_ind).t();
  (objfn_data->e).each_row() += pars_in(objfn_data->d_ind).t();
  objfn_data->Prob = arma::pow(1 + arma::trunc_exp(- (objfn_data->e)), -1);
  (objfn_data->Probs).cols((objfn_data->cat_ind).col(0)) = 1 - objfn_data->Prob;
  (objfn_data->Probs).cols((objfn_data->cat_ind).col(1)) = objfn_data->Prob;
  objfn_data->Probs = arma::trunc_log(objfn_data->Probs);
  objfn_data->e = arma::trunc_exp((objfn_data->e));

  // log-likelihood function global is the sum over the groupwise ones:
  for(int g = 0; g < objfn_data->G; ++g) {
     (objfn_data->ll_g)(g) = arma::accu((objfn_data->rj_g).slice(g) % objfn_data->Probs);
  }

  // gradient of a and d global is the sum over the groupwise ones:
  if(grad_out) {
    for(int g = 0; g < objfn_data->G; ++g) {
      objfn_data->gr_tmp = (((objfn_data->rj_g).slice(g)).cols((objfn_data->cat_ind).col(1)) - (((objfn_data->rj_g).slice(g)).cols((objfn_data->cat_ind).col(0)) % objfn_data->e)) / (objfn_data->e + 1);
      (objfn_data->gr_d).row(g) = arma::sum((objfn_data->gr_tmp), 0);
      (objfn_data->gr_tmp).each_col() %= (objfn_data->X);
      (objfn_data->gr_a).row(g) = arma::sum((objfn_data->gr_tmp), 0);
    }
    (objfn_data->gr_out)(objfn_data->a_ind) = (arma::sum(objfn_data->gr_a, 0)).t();
    (objfn_data->gr_out)(objfn_data->d_ind) = (arma::sum(objfn_data->gr_d, 0)).t();
    *grad_out = - objfn_data->gr_out;
  }
  
  return(- arma::accu(objfn_data->ll_g));
}

