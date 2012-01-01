#ifndef _GROUND_BASE_H_
#define _GROUND_BASE_H_

#include "physics.h"
#include "boost_shared_ptr.h"
#include "geometry.h"

#include <iostream>

#include <vector>
using namespace std;

class Ground {
public:

    virtual ~Ground(){}
    virtual void excludeRegionFromNoise(double xl, double xu, double yl, double yu, double zl, double zu)=0;
    virtual void addNoiseToMeasurements(double_ptr mdata, cdouble_ptr m_pos_vect, csize_t n_measurements, 
					cdouble noise_std_dev) const=0;

    virtual void getDensityAtPoints(double_ptr density_data, cdouble_ptr density_pos_vect, csize_t n, double noise_std_dev) const=0;
    virtual double getAverageDensityInRegion(double xl, double xu, double yl, double yu, double zl, double zu, double noise_std_dev) const=0;
};

template <class NoiseGenerator> class NoisyGround : public Ground
{
public:
    NoisyGround(double xl, double xu, double yl, double yu, double zl, double zu, 
		size_t noise_particle_gridsize_per_meter, shared_ptr<NoiseGenerator> ng_ptr);

    void excludeRegionFromNoise(double xl, double xu, double yl, double yu, double zl, double zu);

    void addNoiseToMeasurements(double_ptr mdata, cdouble_ptr m_pos_vect, csize_t n_measurements, 
				double noise_std_dev) const;

    void getDensityAtPoints(double_ptr density_data, cdouble_ptr density_pos_vect, csize_t n, double noise_std_dev) const;
    
    double getAverageDensityInRegion(double xl, double xu, double yl, double yu, double zl, double zu, double noise_std_dev) const;
    
private:
   
    const Box bounds;

    // Dealing with the indexing of the particles.
    const size_t n_part_x, n_part_y, n_part_z, n_part;
    const double part_vx, part_vy, part_vz;
    const double part_x_start, part_y_start, part_z_start;
    const double particle_density_weight;
    
    inline size_t getIndex(size_t xi, size_t yi, size_t zi) const;
    
    vector<bool>   noise_particle_active;
    vector<double> noise_particle_weight;

    double raw_noise_shift;
    double raw_noise_scale;

    void fillNoiseParticles();

    // The noise generator
    shared_ptr<const NoiseGenerator> noise_generator_ptr;

    vector<Box> excluded_regions;

};

template <class NoiseGenerator> inline size_t 
NoisyGround<NoiseGenerator>::getIndex(size_t xi, size_t yi, size_t zi) const
{
    return zi*n_part_x*n_part_y + yi*n_part_x + xi;
}

#endif /* _GROUND_H_ */
