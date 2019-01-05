#include "var.hpp"

// generate the first N prime numbers
// FIXME faster version / check for primesieve
void sonic::generatePrimes(int *N, arma::vec *primes)
{
  int j = 0, num = 1, pc = 0;
  
  while(j < *N) {
    pc = 0;
    ++num;
    for(int i = 1; i < num; ++i) {
      if(num % i == 0) {
        ++pc;
      }
    }
    if(pc == 1) {
      (*primes)(j) = num;
      ++j;
    }
  }
}

