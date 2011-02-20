#ifndef _GROUND_BASE_HPP_
#define _GROUND_BASE_HPP_

#include "ground_base.h"
#include <algorithm>
#include <iostream>

template <class NoiseGenerator> 
NoisyGround<NoiseGenerator>::NoisyGround(
    double xl, double xu, double yl, double yu, double zl, double zu, 
    size_t noise_particle_gridsize_per_meter, shared_ptr<NoiseGenerator> ng_ptr)

    : bounds(xl, xu, yl, yu, zl, zu),
      n_part_x(noise_particle_gridsize_per_meter*(xu - xl)),
      n_part_y(noise_particle_gridsize_per_meter*(yu - yl)),
      n_part_z(noise_particle_gridsize_per_meter*(zu - zl)),
      n_part(n_part_x * n_part_y * n_part_z),
      part_vx( (xu - xl) / (n_part_x) ),
      part_vy( (yu - yl) / (n_part_y) ),
      part_vz( (zu - zl) / (n_part_z) ),
      part_x_start(xl + 0.5 * part_vx),
      part_y_start(yl + 0.5 * part_vy),
      part_z_start(zl + 0.5 * part_vz),
      particle_density_weight((xu - xl)*(yu - yl)*(zu - zl) / (n_part_x*n_part_y*n_part_z)),
      noise_particle_active(n_part_x*n_part_y*n_part_z, true),
      noise_generator_ptr( (ng_ptr->setup(bounds), ng_ptr) )
{
    // Now all the grid is essentially set up; just get the things in the buffers.
    fillNoiseParticles();
}

template <class NoiseGenerator> 
void NoisyGround<NoiseGenerator>::excludeRegionFromNoise(
    double xl, double xu, double yl, double yu, double zl, double zu)
{
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

    for(size_t zi = z_start; zi <= z_end; ++zi)
    {
	for(size_t yi = y_start; yi <= y_end; ++yi)
	{
	    size_t idx_start = getIndex(x_start, yi, zi);
	    size_t idx_end   = getIndex(x_end,   yi, zi);

	    fill(noise_particle_active.begin() + idx_start, 
		 noise_particle_active.begin() + idx_end + 1, false);
	}
    }

    excluded_regions.push_back(Box(xl, xu, yl, yu, zl, zu));
}

template <class NoiseGenerator> 
void NoisyGround<NoiseGenerator>::addNoiseToMeasurements(
    double_ptr mdata, cdouble_ptr m_pos_vect, csize_t n_measurements, 
    cdouble noise_std_dev) const 
{
    cdouble particle_scaling = particle_density_weight * noise_std_dev;

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(size_t mi = 0; mi < n_measurements; ++mi)
    {
	for(size_t zi = 0; zi < n_part_z; ++zi)
	{
	    for(size_t yi = 0; yi < n_part_y; ++yi)
	    {
		for(size_t xi = 0; xi < n_part_x; ++xi)
		{
		    size_t idx = getIndex(xi, yi, zi);

		    // skip if this noise particle is inactive 
		    if(!noise_particle_active[idx])
			continue;

		    cdouble x = part_x_start + xi*part_vx;
		    cdouble y = part_y_start + yi*part_vy;
		    cdouble z = part_z_start + zi*part_vz;
		
		    // Get the weight of the noise particle.
		    cdouble mass = noise_particle_weight[idx] * particle_scaling;

		    AddGF_Delta(mdata + 3*mi, mdata + 3*mi + 1, mdata + 3*mi + 2, 
				mass,
				x, y, z,
				m_pos_vect[3*mi], m_pos_vect[3*mi+1], m_pos_vect[3*mi + 2]);
		}
	    }
	}
    }
}

template <class NoiseGenerator> 
void NoisyGround<NoiseGenerator>::getDensityAtPoints(
    double_ptr density_data, cdouble_ptr density_pos_vect, csize_t n, double noise_std_dev) const
{   
    
    // first find out which boxes actually intersect the bounding box
    // of the particles here so we don't have to boxes check
    // unnecsarily
 
    const NoiseGenerator& ng = (*noise_generator_ptr);

    Box bb(density_pos_vect[0], density_pos_vect[0],
	   density_pos_vect[1], density_pos_vect[1],
	   density_pos_vect[2], density_pos_vect[2]);

    bb.xl = bb.xu = density_pos_vect[0];
    bb.yl = bb.yu = density_pos_vect[1];
    bb.zl = bb.zu = density_pos_vect[2];

    for(size_t i = 1; i < n; ++i)
    {
	cdouble x = density_pos_vect[3*i];
	cdouble y = density_pos_vect[3*i + 1];
	cdouble z = density_pos_vect[3*i + 2];

	bb.xl = min(x, bb.xl);
	bb.xu = max(x, bb.xu);
	bb.yl = min(y, bb.yl);
	bb.yu = max(y, bb.yu);
	bb.zl = min(z, bb.zl);
	bb.zu = max(z, bb.zu);
    }

    vector<Box> relevant_boxes;
    relevant_boxes.reserve(excluded_regions.size());

    for(size_t i = 0; i < excluded_regions.size(); ++i)
    {
	if(excluded_regions[i].intersects(bb))
	    relevant_boxes.push_back(excluded_regions[i]);
    }

#ifdef USE_OPENMP
#pragma omp parallel for shared(ng, noise_std_dev, relevant_boxes)
#endif

    for(size_t i = 0; i < n; ++i)
    {
	cdouble x = density_pos_vect[3*i];
	cdouble y = density_pos_vect[3*i + 1];
	cdouble z = density_pos_vect[3*i + 2];

	for(size_t j = 0; j < relevant_boxes.size(); ++j)
	{
	    if(relevant_boxes[j].inBox(x,y,z))
	    {
		density_data[i] = 0;
		goto outer_loop_continue_location;
	    }
	}

	density_data[i] = (ng.getNoiseParticleWeight(x, y, z) + raw_noise_shift) * raw_noise_scale;

    outer_loop_continue_location:;
    }
}

template <class NoiseGenerator> 
double NoisyGround<NoiseGenerator>::getAverageDensityInRegion(
    double xl, double xu, double yl, double yu, double zl, double zu, double noise_std_dev) const
{
    if ( xl >= bounds.xu || xu <= bounds.xl
	 || yl >= bounds.yu || yu <= bounds.yl 
	 || zl >= bounds.zu || zu <= bounds.zl)
	return 0;

    const NoiseGenerator& ng = (*noise_generator_ptr);

    const size_t x_start = size_t( ceil ( (xl - bounds.xl) / (bounds.xu - bounds.xl) * n_part_x) );
    const size_t x_end   = size_t( floor( (xu - bounds.xl) / (bounds.xu - bounds.xl) * n_part_x) );
    const size_t y_start = size_t( ceil ( (yl - bounds.yl) / (bounds.yu - bounds.yl) * n_part_y) );
    const size_t y_end   = size_t( floor( (yu - bounds.yl) / (bounds.yu - bounds.yl) * n_part_y) );
    const size_t z_start = size_t( ceil ( (zl - bounds.zl) / (bounds.zu - bounds.zl) * n_part_z) );
    const size_t z_end   = size_t( floor( (zu - bounds.zl) / (bounds.zu - bounds.zl) * n_part_z) );

    double total_noise = 0;

#ifdef USE_OPENMP
#pragma omp parallel for reduction(+:total_noise)
#endif

    for(size_t zi = z_start; zi <= z_end; ++zi)
    {
	for(size_t yi = y_start; yi <= y_end; ++yi)
	{
	    for(size_t xi = x_start; xi <= x_end; ++xi)
	    {
		size_t idx = getIndex(xi,yi,zi);

		if(noise_particle_active[idx])
		{
		    cdouble x = part_x_start + xi*part_vx;
		    cdouble y = part_y_start + yi*part_vy;
		    cdouble z = part_z_start + zi*part_vz;
		
		    if(bounds.inBox(x,y,z)) 
			total_noise += noise_particle_weight[idx];
		    else
			total_noise += (ng.getNoiseParticleWeight(x,y,z) + raw_noise_shift) 
			    * raw_noise_scale;
		}
	    }
	}
    }

    total_noise *= noise_std_dev;

    return total_noise;
}

template <class NoiseGenerator> 
void NoisyGround<NoiseGenerator>::fillNoiseParticles()
{
    const NoiseGenerator& ng = (*noise_generator_ptr);
    
    double total_zi_mean = 0, total_zi_M2 = 0;

    noise_particle_weight.resize(n_part);

    // Do a little extra work to make parallelization easier

#ifdef USE_OPENMP
#pragma omp parallel for shared(ng) reduction(+: total_zi_mean, total_zi_M2)
#endif

    for(size_t zi = 0; zi < n_part_z; zi += 1)
    {
	double local_mean = 0, local_M2 = 0, local_np_idx = 0;

	for(size_t yi = 0; yi < n_part_y; ++yi)
	{
	    for(size_t xi = 0; xi < n_part_x; ++xi)
	    {
		cdouble x = part_x_start + xi * part_vx;
		cdouble y = part_y_start + yi * part_vy;
		cdouble z = part_z_start + zi * part_vz;
		
		// Get the weight of the noise particle.
		cdouble mass = ng.getNoiseParticleWeight(x, y, z);

		noise_particle_weight[getIndex(xi,yi,zi)] = mass;

		// Record the running variance and mean

		++local_np_idx;

		double delta = mass - local_mean;

		local_mean  += delta / local_np_idx;
		local_M2    += delta * (mass - local_mean);
	    }
	}
	
	total_zi_mean += local_mean;
	total_zi_M2   += local_M2 / local_np_idx;
    }

    // Now normalize the noise and recenter it so it is exactly mean
    // zero variance 1.
    raw_noise_shift = -total_zi_mean / n_part_z;
    raw_noise_scale = sqrt(n_part_z / total_zi_M2);

    for(size_t i = 0; i < n_part; ++i)
    {
	noise_particle_weight[i] = (noise_particle_weight[i] + raw_noise_shift) * raw_noise_scale;
    }

    // double true_mean = 0;
    
    // // Now verify the mean and std.
    // for(size_t i = 0; i < n_part; ++i)
    // 	true_mean += noise_particle_weight[i];

    // true_mean /= n_part;
    
    // cout << "True mean = " << true_mean << endl;

    // double true_var = 0;

    // for(size_t i = 0; i < n_part; ++i)
    // 	true_var += noise_particle_weight[i]*noise_particle_weight[i] - true_mean*true_mean;

    // true_var /= n_part;

    // cout << "True Var = " << true_var << endl;

}

#endif
