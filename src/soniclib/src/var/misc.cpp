#include "var.hpp"

// invert hessian
arma::mat sonic::solve(const arma::mat& H)
{
  arma::mat H_ = - H;
  if((H.n_rows == 2) && (H.n_cols == 2)) {
    H_(0, 0) = H(1, 1);
    H_(1, 1) = H(0, 0);
    H_ *= (1 / ((H(0, 0) * H(1, 1)) - (H(0, 1) * H(1, 0))));
  } else {
    //const arma::uword N = H.n_rows;
    //double det_minor = 0;
    //for(arma::uword j = 0; j < N; ++j) {
    //  det_minor = 1 / ((H(2 * j, 2 * j) * H((2 * j) + 1, (2 * j) + 1)) - (H(2 * j, (2 * j) + 1) * H((2 * j) + 1, 2 * j)));
    //  H_(2 * j, 2 * j) = det_minor * H((2 * j) + 1, (2 * j) + 1);
    //  H_(2 * j, (2 * j) + 1) = - det_minor * H(2 * j, (2 * j) + 1);
    //  H_((2 * j) + 1, 2 * j) = - det_minor * H((2 * j) + 1, 2 * j);
    //  H_((2 * j) + 1, (2 * j) + 1) = det_minor * H(2 * j, 2 * j);
    //}

    H_ = arma::inv_sympd(H);

    //arma::mat U;
    //arma::vec s;
    //arma::mat V;
    //arma::svd(U, s, V, H);
    //arma::mat S = arma::diagmat(1 / s);
    //H_ = V * (S * U.t());
  }

  return(H_);
}



// logit vectorwise
arma::vec sonic::logit(arma::vec x)
{
  return(arma::trunc_log(x / (1 - x)));
}

// logit
double sonic::logit_(double x)
{
  return(std::log(x / (1 - x)));
}

// antilogit vectorwise
arma::vec sonic::antilogit(arma::vec x)
{
  return(1 / (1 + arma::trunc_exp(- x)));
}

// antilogit
double sonic::antilogit_(double x)
{
  return(1 / (1 + std::exp(- x)));
}
