/** The main class for inclusion by cython. */
#include "ground.h"
#include "physics.h"
#include "ground_base.hpp"
#include "cutting_planes.h"
#include "boost_shared_ptr.h"

/* The main function for creating a new distribution. */
template<class NoiseGenerator> 
NoisyGround<NoiseGenerator>* MakeNoisyGround(
    double xl, double xu, double yl, double yu, double zl, double zu,
    size_t noise_particle_gridsize_per_meter, 
    size_t noise_order, size_t random_seed)
{
    shared_ptr<NoiseGenerator> ng_ptr(new NoiseGenerator(noise_order, random_seed));

    return new NoisyGround<NoiseGenerator>(
	xl, xu, yl, yu, zl, zu, noise_particle_gridsize_per_meter, ng_ptr);
}

// Now wrap the above

NoisyGround<BinaryCuttingPlanes>* NoisyGround_BinaryCuttingPlanes(
    double xl, double xu, double yl, double yu, double zl, double zu,
    size_t noise_particle_gridsize_per_meter, 
    size_t noise_order, size_t random_seed)
{
    return MakeNoisyGround<BinaryCuttingPlanes>(
	xl,xu,yl,yu,zl,zu, noise_particle_gridsize_per_meter, noise_order, random_seed);
}

NoisyGround<UniformCuttingPlanes>* NoisyGround_UniformCuttingPlanes(
    double xl, double xu, double yl, double yu, double zl, double zu,
    size_t noise_particle_gridsize_per_meter, 
    size_t noise_order, size_t random_seed)
{
    return MakeNoisyGround<UniformCuttingPlanes>(
	xl,xu,yl,yu,zl,zu, noise_particle_gridsize_per_meter, noise_order, random_seed);
}

NoisyGround<NormalCuttingPlanes>* NoisyGround_NormalCuttingPlanes(
    double xl, double xu, double yl, double yu, double zl, double zu,
    size_t noise_particle_gridsize_per_meter, 
    size_t noise_order, size_t random_seed)
{
    return MakeNoisyGround<NormalCuttingPlanes>(
	xl,xu,yl,yu,zl,zu, noise_particle_gridsize_per_meter, noise_order, random_seed);
}

