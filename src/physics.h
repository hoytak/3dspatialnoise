
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
    } while(0)

#endif
