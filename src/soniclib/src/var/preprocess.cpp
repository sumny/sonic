#include "var.hpp"

// generate the first N prime numbers
// FIXME faster version / check for primesieve
arma::vec sonic::generatePrimes(const arma::uword N)
{
  arma::uword j = 0, num = 1, pc = 0;
  arma::vec primes(N, arma::fill::none);

  while(j < N) {
    pc = 0;
    ++num;
    for(arma::uword i = 1; i < num; ++i) {
      if(num % i == 0) {
        ++pc;
      }
    }
    if(pc == 1) {
      primes(j) = num;
      ++j;
    }
  }

  return(primes);
}



// preprocess data
arma::mat sonic::preprocess(const arma::mat& y, const arma::vec& weights,
  const arma::uvec& impact, const arma::uword N, const arma::uword G,
  arma::uword& P, arma::mat& y_u)
{
  // arma::unique does not preserve the original order
  // instead, it sorts in ascending order
  // therefore, reordered data structures are being used internally

  // first N prime numbers
  arma::vec primes = sonic::generatePrimes(N);

  // unique patterns
  arma::vec pats = (y + 1) * arma::log(primes);
  arma::uvec pid_unique = arma::find_unique(pats);
  arma::vec pats_unique = pats(pid_unique);
  P = pid_unique.size();

  // unique patterns groupwise and rgl
  arma::uvec g_id;
  arma::vec weights_g;
  arma::vec pats_g;
  arma::uvec p_vec = arma::linspace<arma::uvec>(0, P - 1, P);
  arma::mat rgl(P, G, arma::fill::none);
  for(arma::uword g = 0; g < G; ++g) {
    g_id = arma::find(impact == g);
    pats_g = pats(g_id);
    weights_g = (weights)(g_id);
    p_vec.for_each( [&g, &weights_g, &pats_g, &pats_unique, &rgl](const arma::uword &l) {
      rgl(l, g) = arma::accu(weights_g(arma::find(pats_g == pats_unique(l))));
    });
  }

  // unique data
  y_u = y.rows(pid_unique);

  return(rgl);
}
