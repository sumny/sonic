#ifndef _preprocess_HPP
#define _preprocess_HPP

// generate the first N prime numbers
arma::vec generatePrimes(const arma::uword N);


// preprocess data
arma::mat preprocess(const arma::mat& y, const arma::vec& weights,
  const arma::uvec& impact, const arma::uword N, const arma::uword G,
  arma::uword& P, arma::mat& y_u);

#endif

