#ifndef _optim_backtracking_HPP
#define _optim_backtracking_HPP

double line_search_bt(double step, arma::vec& x, arma::vec& grad, const arma::vec& direc, const double* wolfe_cons_1_inp, std::function<double (const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data)> opt_objfn, void* opt_data);

#endif
