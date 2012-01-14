#ifndef _GROUND_H_
#define _GROUND_H_

#include "physics.h"
#include "geometry.h"
#include "noise_distributions.h"
#include "typedefs.h"
#include "cutting_planes.hpp"

#include <Eigen/Dense>
#include <gsl/gsl_cdf.h>

template <class WeightGenerator> class NoisyGround 
{
private:
  const Box bounds;

  // Dealing with the indexing of the particles.
  const size_t n_part_x, n_part_y, n_part_z, n_part;
  const double part_vx, part_vy, part_vz;
  const double part_x_start, part_y_start, part_z_start;
    
  Eigen::ArrayXd noise_particle_weight;

  cdouble pm_reflection_boundary, pm_base_length, noise_scale;

  const CuttingPlanes<WeightGenerator> density;

public:
 NoisyGround(cdouble xl, cdouble xu, cdouble yl, cdouble yu, cdouble zl, cdouble zu, 
	     csize_t noise_particle_gridsize_per_meter, 
	     cdouble cutting_plane_linear_density,
	     csize_t random_seed, 
	     cdouble _pm_reflection_boundary, cdouble _pm_base_length)
   : bounds(xl, xu, yl, yu, zl, zu)
   , n_part_x(noise_particle_gridsize_per_meter*(xu - xl))
   , n_part_y(noise_particle_gridsize_per_meter*(yu - yl))
   , n_part_z(noise_particle_gridsize_per_meter*(zu - zl))
   , n_part(n_part_x * n_part_y * n_part_z)
   , part_vx( (xu - xl) / (n_part_x) )
   , part_vy( (yu - yl) / (n_part_y) )
   , part_vz( (zu - zl) / (n_part_z) )
   , part_x_start(xl + 0.5 * part_vx)
   , part_y_start(yl + 0.5 * part_vy)
   , part_z_start(zl + 0.5 * part_vz)
   , noise_particle_weight(n_part)
   , pm_reflection_boundary(_pm_reflection_boundary)
   , pm_base_length(_pm_base_length)
   , noise_scale(pm_reflection_boundary / pm_base_length)

   , density(cutting_plane_linear_density * ( (xu - xl) + (yu - yl) + (zu - zl)),
	     random_seed, Box(xl, xu, yl, yu, zl, zu))
  {

    cdouble particle_density_weight = (xu - xl)*(yu - yl)*(zu - zl) / n_part;

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(size_t i = 0; i < n_part; ++i) {
      size_t xi=0, yi=0, zi=0;

      getIndex(xi, yi, zi, i);

      cdouble x = part_x_start + xi * part_vx;
      cdouble y = part_y_start + yi * part_vy;
      cdouble z = part_z_start + zi * part_vz;

      noise_particle_weight[i] = particle_density_weight * getNoiseMass(x,y,z);
    }
  }

  void excludeRegionFromNoise(double xl, double xu, double yl, double yu, double zl, double zu) {
    // Get all the particles in this region and ensure that they are
    // all those particles are set to be inactive.

    if ( xl >= bounds.xu || xu <= bounds.xl
	 || yl >= bounds.yu || yu <= bounds.yl 
	 || zl >= bounds.zu || zu <= bounds.zl)
	return;

    const size_t x_start = size_t( ceil ( (xl - bounds.xl) / (bounds.xu - bounds.xl) * n_part_x) );
    const size_t x_end   = size_t( floor( (xu - bounds.xl) / (bounds.xu - bounds.xl) * n_part_x) );
    const size_t y_start = size_t( ceil ( (yl - bounds.yl) / (bounds.yu - bounds.yl) * n_part_y) );
    const size_t y_end   = size_t( floor( (yu - bounds.yl) / (bounds.yu - bounds.yl) * n_part_y) );
    const size_t z_start = size_t( ceil ( (zl - bounds.zl) / (bounds.zu - bounds.zl) * n_part_z) );
    const size_t z_end   = size_t( floor( (zu - bounds.zl) / (bounds.zu - bounds.zl) * n_part_z) );

    for(size_t zi = z_start; zi <= z_end; ++zi) {
      for(size_t yi = y_start; yi <= y_end; ++yi) {

	size_t idx_start = getIndex(x_start, yi, zi);
	size_t idx_end   = getIndex(x_end,   yi, zi);

	for(size_t i = idx_start; i <= idx_end; ++i)
	  noise_particle_weight[i] = 0;
      }
    }
  }
  
  void addNoiseToMeasurements(double_ptr mdata, cdouble_ptr m_pos_vect, csize_t n_measurements) const {

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(size_t idx = 0; idx < n_part; ++idx) {

      if(noise_particle_weight[idx] == 0)
	continue;

      size_t xi=0, yi=0, zi=0;

      getIndex(xi, yi, zi, idx);

      cdouble x = part_x_start + xi*part_vx;
      cdouble y = part_y_start + yi*part_vy;
      cdouble z = part_z_start + zi*part_vz;

      for(size_t mi = 0; mi < n_measurements; ++mi) {

	// Get the weight of the noise particle.
	cdouble mass = noise_particle_weight[idx];

	AddGF_Delta(mdata + 3*mi, mdata + 3*mi + 1, mdata + 3*mi + 2, 
		    mass,
		    x, y, z,
		    m_pos_vect[3*mi], m_pos_vect[3*mi+1], m_pos_vect[3*mi + 2]);
      }
    }
  }
    
  double getAverageDensityInRegion(double xl, double xu, double yl, double yu, double zl, double zu) const
  {
    if ( xl >= bounds.xu || xu <= bounds.xl
	 || yl >= bounds.yu || yu <= bounds.yl 
	 || zl >= bounds.zu || zu <= bounds.zl)
	return 0;

    csize_t x_start = size_t( ceil ( (xl - bounds.xl) / (bounds.xu - bounds.xl) * n_part_x) );
    csize_t x_end   = size_t( floor( (xu - bounds.xl) / (bounds.xu - bounds.xl) * n_part_x) );
    csize_t y_start = size_t( ceil ( (yl - bounds.yl) / (bounds.yu - bounds.yl) * n_part_y) );
    csize_t y_end   = size_t( floor( (yu - bounds.yl) / (bounds.yu - bounds.yl) * n_part_y) );
    csize_t z_start = size_t( ceil ( (zl - bounds.zl) / (bounds.zu - bounds.zl) * n_part_z) );
    csize_t z_end   = size_t( floor( (zu - bounds.zl) / (bounds.zu - bounds.zl) * n_part_z) );

    csize_t nx = x_end - x_start + 1;
    csize_t ny = y_end - y_start + 1;
    csize_t nz = z_end - z_start + 1;

    double total_density = 0;

#ifdef USE_OPENMP
#pragma omp parallel for reduction(+ : total_density)
#endif
    for(size_t index = 0; index < nx*ny*nz; ++index) {
      csize_t zi = index / (nx * ny); 
      csize_t yi = (index % (nx * ny)) / nx;
      csize_t xi = index % nx;

      total_density += noise_particle_weight[getIndex(xi, yi, zi)];
    }

    return total_density / ( (xu - xl) * (yu - yl) * (zu - zl));
  }

  void getDensityAtPoints(double_ptr d, double_cptr X, size_t n) const {
    
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(size_t i = 0; i <= n; ++i) {
      d[i] = getNoiseMass(X[3*i + 0], X[3*i + 1], X[3*i + 2]);
    }
  }


private:

  inline double getNoiseMass(double x, double y, double z) const {

    double mass = noise_scale * density(x, y, z);

    // if (mass < -pm_reflection_boundary)
    //   return -pm_reflection_boundary;
    // if (mass > pm_reflection_boundary)
    //   return pm_reflection_boundary;

    // return mass;

    cdouble& a = pm_reflection_boundary;
		
    cdouble mod_abs_mass = fmod(abs(mass + a), 4*a);
      
    double rval = ((mod_abs_mass > 2*a) ? (4*a - mod_abs_mass) : mod_abs_mass) - a;

    return rval;
  }

  inline size_t getIndex(size_t xi, size_t yi, size_t zi) const {
    return zi*n_part_x*n_part_y + yi*n_part_x + xi;
  }

  inline void getIndex(size_t &xi, size_t &yi, size_t &zi, size_t index) const {
    zi = index / (n_part_x * n_part_y); 
    yi = (index % (n_part_x * n_part_y)) / n_part_x;
    xi = index % n_part_x;
  }

};

#endif /* _GROUND_H_ */
