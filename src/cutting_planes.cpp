#include "cutting_planes.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cmath>

CuttingPlanes::CuttingPlanes(size_t _num_cutting_planes, size_t _random_seed)
    : num_cutting_planes(_num_cutting_planes), random_seed(_random_seed)
{}

void CuttingPlanes::setup(const Box& b)
{
    /* Init the random number generator. */
    gsl_rng_env_setup();

    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng *r = gsl_rng_alloc (T);
    gsl_rng_set(r, random_seed);
   
    /* Go through and create however many cutting planes we need. */
    cutting_planes.reserve(num_cutting_planes);

    for(size_t i = 0; i < num_cutting_planes; ++i)
    {
	double x, y, z;
	double vx, vy, vz;

	x = gsl_ran_flat(r, b.xl, b.xu);
	y = gsl_ran_flat(r, b.yl, b.yu);
	z = gsl_ran_flat(r, b.zl, b.zu);
	gsl_ran_dir_3d(r, &vx, &vy, &vz);

	cutting_planes.push_back(Plane(Ray(Point(x,y,z),vx,vy,vz)));
    }

    /* Now go set the weights.*/
    setCuttingPlaneWeights(num_cutting_planes, r);

    gsl_rng_free(r);
}

////////////////////////////////////////////////////////////////////////////////
// Now for specific weights used on the cutting planes

void BinaryCuttingPlanes::setCuttingPlaneWeights(csize_t num_cutting_planes, gsl_rng *)
{
    // Simply set the cutting planes with weights either -1 or 1.
    // This results in setting the weights to 1 / \sqrt(\sum_i Var(X_i)) = 1 / \sqrt(n)

    weights.resize(num_cutting_planes);
    fill(weights.begin(), weights.end(), 1.0 / sqrt(double(num_cutting_planes)));
}

void NormalCuttingPlanes::setCuttingPlaneWeights(csize_t num_cutting_planes, gsl_rng *r)
{
    weights.resize(num_cutting_planes);
    
    double scale_factor = 1.0 / sqrt(double(num_cutting_planes));

    for(size_t i = 0; i < num_cutting_planes; ++i)
	weights[i] = gsl_ran_gaussian(r, 1) * scale_factor;

}

void UniformCuttingPlanes::setCuttingPlaneWeights(csize_t num_cutting_planes, gsl_rng *r)
{
    weights.resize(num_cutting_planes);
    
    double scale_factor = (4.0 / 12) / sqrt(double(num_cutting_planes));

    for(size_t i = 0; i < num_cutting_planes; ++i)
	weights[i] = gsl_ran_flat(r, -1, 1) * scale_factor;

}
