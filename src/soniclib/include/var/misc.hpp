#ifndef _misc_HPP
#define _misc_HPP

// invert hessian
arma::mat solve(const arma::mat& H);

// logit vectorwise
arma::vec logit(arma::vec x);

// logit
double logit_(double x);

// antilogit vectorwise
arma::vec antilogit(arma::vec x);

// antilogit
double antilogit_(double x);

#endif
