#include "optim.hpp"

double optim::line_search_bt(double step, arma::vec& x, arma::vec& grad, const arma::vec& direc, const double* wolfe_cons_1_inp, std::function<double (const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data)> opt_objfn, void* opt_data)
{
    const double rho = 0.5; // factor of changing the step length if the Armijo sufficient descrease condition is not met

    const double wolfe_cons_1 = (wolfe_cons_1_inp) ? *wolfe_cons_1_inp : 1E-03; // tolerence on the Armijo sufficient decrease condition

    arma::vec x_0 = x;

    // backtracking
    while(opt_objfn(x + step * direc, nullptr, opt_data) > (opt_objfn(x, nullptr, opt_data) + wolfe_cons_1 * step * arma::dot(grad, direc))) {
      step = rho * step;
    }

    // update x and grad
    x = x + step * direc;
    opt_objfn(x, &grad, opt_data);

    return(step);
}
