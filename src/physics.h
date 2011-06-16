
#ifndef PHYSICS_H
#define PHYSICS_H

#define G 1 //6.67300e-11

#include "typedefs.h"
#include <cmath>

#define AddGF_Delta(Txx_ptr, Txy_ptr, Tyy_ptr,				\
		    mass,						\
		    px, py, pz,						\
		    mx, my, mz)						\
    do{									\
	cdouble dx = (px) - (mx);					\
	cdouble dy = (py) - (my);					\
	cdouble dz = (pz) - (mz);					\
	cdouble dl = dx*dx + dy*dy + dz*dz;				\
	cdouble common_coeff = (mass)*G*(pow(dl, -2.5));		\
	cdouble Txx_diff = (3*dx*dx - dl)*common_coeff;			\
	cdouble Txy_diff = 3*dx*dy*common_coeff;			\
	cdouble Tyy_diff = (3*dy*dy - dl)*common_coeff;			\
	*(Txx_ptr) += Txx_diff;						\
	*(Txy_ptr) += Txy_diff;						\
	*(Tyy_ptr) += Tyy_diff;						\
									\
	/*cout << "px=" << px << "; py=" << py << "; pz=" << pz		\
	     << "; mx=" << mx << "; my=" << my << "; mz=" << mz		\
	     << "; Txx=" << Txx_diff << "; Txy=" << Txy_diff << "; Tyy=" << Tyy_diff << endl; \
	*/								\
    }while(0)

void addDeltaToMeasurements(double_ptr mdata,
			    cdouble_ptr m_pos_vect,
			    csize_t n,
			    cdouble mass, 
			    cdouble x, cdouble y, cdouble z);

void addDeltaVectToMeasurements(double_ptr mdata,
				cdouble_ptr m_pos_vect,
				csize_t n_measurements,
				cdouble_ptr mass_vect, 
				cdouble_ptr delta_pos_vect,
				csize_t n_locations);

void addBoxToMeasurements(double_ptr mdata, 
			  cdouble_ptr m_pos_vect,
			  csize_t n_measurements, 
			  csize_t particle_grid_size, 
			  cdouble box_mass, 
			  cdouble x_min, cdouble x_max,
			  cdouble y_min, cdouble y_max,
			  cdouble z_min, cdouble z_max);

void setLatticeCoefficients(double_ptr coeff_vect,
			    cdouble_ptr m_pos_vect,
			    csize_t n_measurements,
			    csize_t particle_grid_size,
			    cdouble_ptr min_vect, 
			    cdouble_ptr max_vect,
			    csize_t n_voxels,
			    cdouble density);

void setPointfieldCoefficients(double_ptr coeff_vect,
			       cdouble_ptr m_pos_vect,
			       cdouble_ptr p_pos_vect,
			       cdouble_ptr p_points,
			       csize_t n_measurements,
			       csize_t n_points);

#endif
