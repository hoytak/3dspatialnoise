#include "physics.h"

#include <cstdlib>
#include <iostream>
#include <algorithm>

using namespace std;

void addDeltaToMeasurements(double_ptr mdata,
			    cdouble_ptr m_pos_vect,
			    csize_t n,
			    cdouble mass, 
			    cdouble x, cdouble y, cdouble z)
{
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(size_t i = 0; i < n; ++i)
    {
	AddGF_Delta(mdata + 3*i, mdata + 3*i + 1, mdata + 3*i + 2, 
		    mass, 
		    x, y, z,
		    m_pos_vect[3*i], m_pos_vect[3*i+1], m_pos_vect[3*i + 2]);
    }
}

void addDeltaVectToMeasurements(double_ptr mdata,
				cdouble_ptr m_pos_vect,
				csize_t n_measurements,
				cdouble_ptr mass_vect, 
				cdouble_ptr delta_pos_vect,
				csize_t n_locations)
{
    size_t i, j;

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(i = 0; i < n_measurements; ++i)
    {
	for(j = 0; j < n_locations; ++j)
	{
	    AddGF_Delta(mdata + 3*i, mdata + 3*i + 1, mdata + 3*i + 2, 
			mass_vect[j],
			delta_pos_vect[3*j], delta_pos_vect[3*j+1], delta_pos_vect[3*j+2],
			m_pos_vect[3*i], m_pos_vect[3*i+1], m_pos_vect[3*i + 2]);
	}
    }
}

void addBoxToMeasurements(double_ptr mdata, 
			  cdouble_ptr m_pos_vect,
			  csize_t n_measurements, 
			  csize_t particle_grid_size, 
			  cdouble box_mass, 
			  cdouble x_min, cdouble x_max,
			  cdouble y_min, cdouble y_max,
			  cdouble z_min, cdouble z_max)
{
    cdouble vx = (x_max - x_min) / (particle_grid_size);
    cdouble x = x_min + vx/2;

    cdouble vy = (y_max - y_min) / (particle_grid_size);
    cdouble y = y_min + vy/2;

    cdouble vz = (z_max - z_min) / (particle_grid_size);
    cdouble z = z_min + vz/2;

    size_t xi, yi, zi;

    cdouble scaled_mass = box_mass / (particle_grid_size*particle_grid_size*particle_grid_size);

    for(xi = 0; xi < particle_grid_size; ++xi)
    {
	for(yi = 0; yi < particle_grid_size; ++yi)
	{
	    for(zi = 0; zi < particle_grid_size; ++zi)
	    {
                addDeltaToMeasurements(mdata, m_pos_vect,
				       n_measurements, scaled_mass, 
				       x + vx*xi, y + vy*yi, z + vz*zi);
	    }
	}
    }
}

void setLatticeCoefficients(double_ptr coeff_vect,
			    cdouble_ptr m_pos_vect,
			    csize_t n_measurements,
			    csize_t particle_grid_size,
			    cdouble_ptr min_vect, 
			    cdouble_ptr max_vect,
			    csize_t n_voxels,
			    cdouble density)
{
    fill(coeff_vect, coeff_vect + 3*n_measurements*n_voxels, 0);

#ifdef USE_OPENMP
#pragma omp parallel for 
#endif
    for(size_t mi = 0; mi < n_measurements; ++mi)
    {
	for(size_t vi = 0; vi < n_voxels; ++vi)
	{
	    cdouble x_min = min_vect[3*vi+0];
	    cdouble x_max = max_vect[3*vi+0];
	    cdouble y_min = min_vect[3*vi+1];
	    cdouble y_max = max_vect[3*vi+1];
	    cdouble z_min = min_vect[3*vi+2];
	    cdouble z_max = max_vect[3*vi+2];

	    cdouble xw = x_max - x_min;
	    cdouble yw = y_max - y_min;
	    cdouble zw = z_max - z_min;

	    cdouble volume = xw*yw*zw;

	    cdouble vx = xw / (particle_grid_size);
	    cdouble x  = x_min + vx/2;

	    cdouble vy = yw / (particle_grid_size);
	    cdouble y  = y_min + vy/2;

	    cdouble vz = zw / (particle_grid_size);
	    cdouble z  = z_min + vz/2;
	    /*
	    cout << "volume*density = " << (volume*density) << endl;
	    cout << "particle_grid_size" << particle_grid_size << endl;
	    */
	    cdouble scaled_mass = volume*density / (particle_grid_size*particle_grid_size*particle_grid_size);
	    
	    /*
	    cout << "(x:" << x_min << "," << x_max << ")"
		 << "(y:" << y_min << "," << y_max << ")"
		 << "(z:" << z_min << "," << z_max << ") = "
		 << scaled_mass << endl;
	    */

	    for(size_t xi = 0; xi < particle_grid_size; ++xi)
	    {
		for(size_t yi = 0; yi < particle_grid_size; ++yi)
		{
		    for(size_t zi = 0; zi < particle_grid_size; ++zi)
		    {
			/*
			cout << "Sampling at " << x + vx*xi 
			     << ", " << y + vy*yi 
			     << ", " << z + vz*zi << endl;

			cout << "Measurement to ("  
			     << m_pos_vect[3*mi] << ", " 
			     << m_pos_vect[3*mi+1] << ", "
			     << m_pos_vect[3*mi + 2] << ") " << endl;
			*/
			AddGF_Delta(coeff_vect + 3*n_voxels*mi + 3*vi + 0,
				    coeff_vect + 3*n_voxels*mi + 3*vi + 1,
				    coeff_vect + 3*n_voxels*mi + 3*vi + 2,
				    scaled_mass,
				    x + vx*xi, y + vy*yi, z + vz*zi,
				    m_pos_vect[3*mi], m_pos_vect[3*mi+1], m_pos_vect[3*mi + 2]);
			/*
			cout << "Txx = " << *(coeff_vect + 3*n_voxels*mi + 3*vi + 0)
			     << ", Txy = " << *(coeff_vect + 3*n_voxels*mi + 3*vi + 1)
			     << ", Tyy = " << *(coeff_vect + 3*n_voxels*mi + 3*vi + 2)
			     << endl;
			*/
		    }
		}
	    }
	}
    }
}

void setPointfieldCoefficients(double_ptr coeff_vect,
			       cdouble_ptr m_pos_vect,
			       cdouble_ptr p_pos_vect,
			       cdouble_ptr p_mass,
			       csize_t n_measurements,
			       csize_t n_points)
{
    fill(coeff_vect, coeff_vect + 3*n_measurements*n_points, 0);

#ifdef USE_OPENMP
#pragma omp parallel for 
#endif
    for(size_t mi = 0; mi < n_measurements; ++mi)
    {
	double mx = m_pos_vect[3*mi + 0];
	double my = m_pos_vect[3*mi + 1];
	double mz = m_pos_vect[3*mi + 2];

	for(size_t pi = 0; pi < n_points; ++pi)
	{
	    double px = p_pos_vect[3*pi + 0];
	    double py = p_pos_vect[3*pi + 1];
	    double pz = p_pos_vect[3*pi + 2];

	    coeff_vect[3*n_points*mi + 3*pi + 0] = 0;
	    coeff_vect[3*n_points*mi + 3*pi + 1] = 0;
	    coeff_vect[3*n_points*mi + 3*pi + 2] = 0;
	    
	    AddGF_Delta(coeff_vect + 3*n_points*mi + 3*pi + 0,
			coeff_vect + 3*n_points*mi + 3*pi + 1,
			coeff_vect + 3*n_points*mi + 3*pi + 2,
			p_mass[pi],
			px, py, pz, mx, my, mz);
	    
	}
    }
}

