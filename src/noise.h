#ifndef _NOISE_H_
#define _NOISE_H_

/** The main class for inclusion by cython. */

#include "physics.h"
#include "ground_base.hpp"
#include "cutting_planes.h"
#include "boost_shared_ptr.h"

/* The main function for creating a new distribution. */
template<class NoiseGenerator> 
Ground* MakeNoisyGround(
    double xl, double xu, double yl, double yu, double zl, double zu,
    size_t noise_particle_gridsize_per_meter, 
    size_t noise_order, size_t random_seed)
{
    shared_ptr<NoiseGenerator> ng_ptr(new NoiseGenerator(noise_order, random_seed));

    return new NoisyGround<NoiseGenerator>(
	xl, xu, yl, yu, zl, zu, noise_particle_gridsize_per_meter, ng_ptr);
}

#endif /* _NOISE_H_ */
