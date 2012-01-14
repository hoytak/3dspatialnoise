#ifndef _NOISE_DISTRIBUTIONS_H_
#define _NOISE_DISTRIBUTIONS_H_

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cmath>

class BinaryNoiseDistribution {
 public:
 BinaryNoiseDistribution(gsl_rng *_r) 
   : r(_r)
    {}

  double operator()() const {
    return ((gsl_ran_bernoulli(r, 0.5) != 0) ? -1 : 1);
  }

 private:
  gsl_rng *r;
};


class NormalNoiseDistribution {
 public:
 NormalNoiseDistribution(gsl_rng *_r) 
   : r(_r)
    {}

  double operator()() const {
    return gsl_ran_gaussian(r, 1);
  }

 private:
  gsl_rng *r;
};

class UniformNoiseDistribution {
 public:
 UniformNoiseDistribution(csize_t num_cutting_planes, gsl_rng *_r)
   : r(_r)
    {}

  double operator()() const {
    return gsl_ran_flat(r, -1, 1);
  }

 private:
  gsl_rng *r;
};


#endif /* _NOISE_DISTRIBUTIONS_H_ */
