#ifndef _CUTTING_PLANS_H_
#define _CUTTING_PLANS_H_

#include "geometry.h"
#include "typedefs.h"

#include <vector>
#include <utility>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>
#include <Eigen/Dense>

using Eigen::ArrayXd;
using std::vector;
using std::pair;
using std::make_pair;

template <class WeightGenerator> class CuttingPlanes 
{
public:
  CuttingPlanes(size_t _num_cutting_planes, size_t random_seed, const Box& world)
    : plx(_num_cutting_planes), ply(_num_cutting_planes), plz(_num_cutting_planes)
    , plvx(_num_cutting_planes), plvy(_num_cutting_planes), plvz(_num_cutting_planes)
    , weights(_num_cutting_planes)
    , num_cutting_planes(_num_cutting_planes)
  {
    gsl_rng_env_setup();

    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng *r = gsl_rng_alloc (T);
    gsl_rng_set(r, random_seed);

    // First, generate the weights
    WeightGenerator wg(num_cutting_planes, r);
    for(size_t i = 0; i < num_cutting_planes; ++i)
      weights[i] = wg();

    // Now all the points
    Point center = world.center();
    Point corner = world.centeredCorner();

    cdouble R = dist(center, corner);
    
    for(size_t i = 0; i < num_cutting_planes; ++i)
    {
    plane_invalid:

      // Choose a point randomly in the containing sphere
      cdouble radius = R*pow(gsl_ran_flat(r, 0, 1), 1.0 / 3);

      double vx, vy, vz;     
      gsl_ran_dir_3d(r, &vx, &vy, &vz);

      // does it intersect our world?
      double abs_r_vx = radius * abs(vx);
      double abs_r_vy = radius * abs(vy);
      double abs_r_vz = radius * abs(vz);

      Plane cp = Plane(Ray(Point(abs_r_vx,abs_r_vy,abs_r_vz),
			   abs_r_vx,abs_r_vy,abs_r_vz));
	
      if(!cp.onFrontSide(corner))
	goto plane_invalid;

      // Otherwise, 
      plx[i] = radius*vx;
      ply[i] = radius*vy;
      plz[i] = radius*vz;
      plvx[i] = vx;
      plvy[i] = vy;
      plvz[i] = vz;
    }

    // Now I just have to callibrate the weights; choose a bunch of
    // random points; get the standard deviation of difference, then scale so
    // it's 1 per meter

    const size_t N = 500;
    double values[N];

    vector<pair<Point, Point> > pos;
    pos.reserve(N);

    for(size_t i = 0; i < N; ++i) {
      pos.push_back(make_pair(Point(gsl_ran_flat(r, world.xl, world.xu),
				    gsl_ran_flat(r, world.yl, world.yu),
				    gsl_ran_flat(r, world.zl, world.zu)),
			      Point(gsl_ran_flat(r, world.xl, world.xu),
				    gsl_ran_flat(r, world.yl, world.yu),
				    gsl_ran_flat(r, world.zl, world.zu))));
    }

    gsl_rng_free(r);   

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(size_t i = 0; i < N; ++i) {
      const Point& p1 = pos[i].first;
      const Point& p2 = pos[i].second;

      double d = dist(p1, p2);

      values[i] = ( (*this)(p1.x, p1.y, p1.z) - (*this)(p2.x, p2.y, p2.z) ) / d;
    }

    // Rescale the weights
    weights /= gsl_stats_sd_with_fixed_mean(values, 1, N, 0);
  }

  inline double operator()(cdouble x, cdouble y, cdouble z) const {

    ArrayXd dist = (x - plx)*plvx + (y - ply)*plvy + (z - plz)*plvz;

    double density = 0;

    for(size_t i = 0; i < num_cutting_planes; ++i) 
      density += ( (dist[i] > 0) - (dist[i] < 0) ) * weights[i];

    return density;
  }

private:
  ArrayXd plx, ply, plz, plvx, plvy, plvz, weights;
  csize_t num_cutting_planes;
};

#endif /* _CUTTING_PLANS_H_ */
