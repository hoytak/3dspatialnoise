#ifndef _NOISE_DISTRIBUTIONS_H_
#define _NOISE_DISTRIBUTIONS_H_

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cmath>

class BinaryNoiseDistribution {
 public:
 BinaryNoiseDistribution(csize_t num_cutting_planes, gsl_rng *_r) 
   : scale_factor(1.0 / sqrt(double(num_cutting_planes))), r(_r)
    {}

  double operator()() const {
    return ((gsl_ran_bernoulli(r, 0.5) == 0) ? -1 : 1)*scale_factor;
  }

 private:
  double scale_factor;
  gsl_rng *r;
};


class NormalNoiseDistribution {
 public:
 NormalNoiseDistribution(csize_t num_cutting_planes, gsl_rng *_r) 
   : scale_factor(1.0 / sqrt(double(num_cutting_planes))), r(_r)
    {}

  double operator()() const {
    return gsl_ran_gaussian(r, 1) * scale_factor;
  }

 private:
  double scale_factor;
  gsl_rng *r;
};

class UniformNoiseDistribution {
 public:
 UniformNoiseDistribution(csize_t num_cutting_planes, gsl_rng *_r)
   : scale_factor((4.0 / 12) / sqrt(double(num_cutting_planes))), r(_r)
    {}

  double operator()() const {
    return gsl_ran_flat(r, -1, 1) * scale_factor;
  }

 private:
  double scale_factor;
  gsl_rng *r;
};


#endif /* _NOISE_DISTRIBUTIONS_H_ */
