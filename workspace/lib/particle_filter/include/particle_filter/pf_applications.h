#pragma once

#include <particle_filter/detail/particle_filter.h>
#include <math/Types.h>
#include <math/Geometry.h>
#include <utility>

namespace pf {

/**
 * Utility functions for particle filter applications with field of view constraints
 */
class PFApplications {
public:
    /**
     * Process a particle filter update with field-of-view constraints and return state estimate
     * 
     * @param filter The particle filter to process
     * @param ego_pos Position of the ego robot (observer)
     * @param neighbor_pos Position of the neighbor robot (if visible)
     * @param fov_beta Field of view angle in radians
     * @param fov_Rs Field of view sensing range
     * @param weight_reduction_factor Factor to reduce weights for particles in FOV
     * @return A pair containing the state estimate vector and covariance matrix
     */
    template<typename T, unsigned int DIM>
    static std::pair<math::Vector<double>, math::Matrix<double>> 
    processFovUpdate(ParticleFilter& filter, 
                     const math::VectorDIM<T, DIM>& ego_pos,
                     const math::VectorDIM<T, DIM>& neighbor_pos,
                     double fov_beta, 
                     double fov_Rs,
                     double weight_reduction_factor = 10.0);

};

} // namespace pf
