#ifndef _GROUND_H_
#define _GROUND_H_

/** The main class for inclusion by cython. */

#include "ground.h"
#include "physics.h"
#include "ground_base.h"
#include "cutting_planes.h"
#include "boost_shared_ptr.h"

// Now wrap the above

NoisyGround<BinaryCuttingPlanes>* NoisyGround_BinaryCuttingPlanes(
    double xl, double xu, double yl, double yu, double zl, double zu,
    size_t noise_particle_gridsize_per_meter, 
    size_t noise_order, size_t random_seed);

NoisyGround<UniformCuttingPlanes>* NoisyGround_UniformCuttingPlanes(
    double xl, double xu, double yl, double yu, double zl, double zu,
    size_t noise_particle_gridsize_per_meter, 
    size_t noise_order, size_t random_seed);

NoisyGround<NormalCuttingPlanes>* NoisyGround_NormalCuttingPlanes(
    double xl, double xu, double yl, double yu, double zl, double zu,
    size_t noise_particle_gridsize_per_meter, 
    size_t noise_order, size_t random_seed);

#endif /* _NOISE_H_ */
