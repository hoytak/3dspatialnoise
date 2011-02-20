#ifndef _CUTTING_PLANS_H_
#define _CUTTING_PLANS_H_

#include "ground_base.h"
#include "geometry.h"
#include "typedefs.h"

#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

class CuttingPlanes 
{
public:
    CuttingPlanes(size_t num_cutting_planes, size_t random_seed);
    void setup(const Box& world); 

    inline double getNoiseParticleWeight(cdouble x, cdouble y, cdouble z) const; 

protected:
    virtual void setCuttingPlaneWeights(csize_t num_cutting_planes, gsl_rng *r)=0;
    vector<double> weights;

private:
    vector<Plane> cutting_planes;
    csize_t num_cutting_planes;
    csize_t random_seed;
};

/********************************************************************************/
// Now specific types of cutting planes

class BinaryCuttingPlanes : public CuttingPlanes 
{
public:
BinaryCuttingPlanes(size_t num_cutting_planes, size_t random_seed)
    : CuttingPlanes(num_cutting_planes, random_seed) 
    {}
protected:
    void setCuttingPlaneWeights(csize_t num_cutting_planes, gsl_rng *r);
};

class NormalCuttingPlanes : public CuttingPlanes 
{
public:
NormalCuttingPlanes(size_t num_cutting_planes, size_t random_seed)
    : CuttingPlanes(num_cutting_planes, random_seed) 
    {}
    
protected:
    void setCuttingPlaneWeights(csize_t num_cutting_planes, gsl_rng *r);
};

class UniformCuttingPlanes : public CuttingPlanes 
{
public:
UniformCuttingPlanes(size_t num_cutting_planes, size_t random_seed)
    : CuttingPlanes(num_cutting_planes, random_seed) 
    {}

protected:
    void setCuttingPlaneWeights(csize_t num_cutting_planes, gsl_rng *r);
};

////////////////////////////////////////////////////////////////////////////////
// A few inline functions of the above.

inline double CuttingPlanes::getNoiseParticleWeight(cdouble x, cdouble y, cdouble z) const
{
    double mass = 0;

    for(size_t j = 0; j < cutting_planes.size(); ++j)
    {
	const int w = (cutting_planes[j].onFrontSide(Point(x, y, z))) ? 1 : -1;
	mass += w*weights[j];
    }

    return mass;
}

#endif /* _CUTTING_PLANS_H_ */
