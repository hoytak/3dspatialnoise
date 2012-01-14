from numpy cimport ndarray as ar
from numpy import empty, zeros, asarray


cdef extern from "ground.hpp":

    cdef cppclass Ground "NoisyGround<BinaryNoiseDistribution>":

        Ground(double xl, double xu, double yl, double yu, double zl, double zu, 
               size_t noise_particle_gridsize_per_meter, 
               double cutting_plane_linear_density,
               size_t random_seed, 
               double pm_reflection_boundary, double pm_crossing_rate)

        void excludeRegionFromNoise(double xl, double xu, double yl, double yu, double zl, double zu)

        void getDensityAtPoints(double *res, double *pos, size_t n)

        void addNoiseToMeasurements(double * mdata, double * m_pos_vect, size_t n_measurements)

        double getAverageDensityInRegion(double xl, double xu, double yl, double yu, double zl, double zu)

cdef class NoisyGround:

    cdef Ground *ground

    def __init__(self, npt):

        # get the bounding box

        cdef double xl, xu, yl, yu, zl, zu
        bb = npt.bounding_box
        xl, xu, yl, yu, zl, zu = asarray(bb).ravel()

        self.ground = new Ground(
            xl, xu, yl, yu, zl, zu,
            npt.gridsize_per_meter,
            npt.cutting_plane_linear_density,
            npt.random_seed,
            npt.pm_reflection_boundary,
            npt.pm_crossing_rate)

    def __dealloc__(self):
        del self.ground

    def excludeRegionFromNoise(self, bb):

        cdef ar[double] bounds = asarray(bb).ravel()

        cdef double xl = bounds[0]
        cdef double xu = bounds[1]
        cdef double yl = bounds[2]
        cdef double yu = bounds[3]
        cdef double zl = bounds[4]
        cdef double zu = bounds[5]

        self.ground.excludeRegionFromNoise(xl,xu,yl,yu,zl,zu)

    def noiseDensityInRegion(self, bb):

        cdef ar[double] bounds = asarray(bb).ravel()

        cdef double xl = bounds[0]
        cdef double xu = bounds[1]
        cdef double yl = bounds[2]
        cdef double yu = bounds[3]
        cdef double zl = bounds[4]
        cdef double zu = bounds[5]

        return self.ground.getAverageDensityInRegion(xl,xu,yl,yu,zl,zu)

    def getNoiseDensity(self, ar[double, ndim=2, mode = "c"] X):
        """
        Returns the un-shifted noise distribution.
        """

        assert X.shape[1] == 3

        cdef ar[double, mode="c"] density = empty( X.shape[0])

        self.ground.getDensityAtPoints(<double*>density.data, <double*>X.data, X.shape[0])

        return density

    def noiseAtPoints(self, ar[double, ndim=2, mode='c'] mpts):

        cdef ar[double, ndim=2, mode="c"] T = zeros( (mpts.shape[0], 3) )

        self.ground.addNoiseToMeasurements(<double*>(T.data), <double*>(mpts.data), mpts.shape[0])

        return T
