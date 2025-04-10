#include <particle_filter/pf_applications.h>

namespace pf {

// Implementation of template method
template<typename T, unsigned int DIM>
std::pair<math::Vector<double>, math::Matrix<double>> 
PFApplications::processFovUpdate(ParticleFilter& filter, 
                                const math::VectorDIM<T, DIM>& ego_pos,
                                const math::VectorDIM<T, DIM>& neighbor_pos,
                                double fov_beta, 
                                double fov_Rs,
                                double weight_reduction_factor) {
    
    // Propagate particles forward according to motion model
    filter.predict();
    
    // Get and update particle weights based on field of view
    math::Vector<double> weights = filter.getWeights();
    math::Matrix<double> samples = filter.getParticles();
    
    for (int s = 0; s < samples.cols(); ++s) {
        // Reduce weight of particles that are in robot's field of view
        // but don't match observation (implicit negative information)
        if (math::insideFOV(ego_pos, samples.col(s), fov_beta, fov_Rs)) {
            weights[s] /= weight_reduction_factor;
        }
    }
    filter.setWeights(weights);

    // Update filter with new measurement if neighbor is visible
    if (math::insideFOV(ego_pos, neighbor_pos, fov_beta, fov_Rs)) {
        math::Vector<double> neighbor_xy(DIM-1);
        neighbor_xy << neighbor_pos(0), neighbor_pos(1);
        filter.update(neighbor_xy);
    }

    // Resample particles and compute state estimate
    filter.resample();
    filter.estimateState();

    // Extract state estimate and covariance from the filter
    math::Vector<double> estimate = filter.getState();
    math::Matrix<double> cov = filter.getDistribution();
    
    return std::make_pair(estimate, cov);
}

// Explicit instantiations for common dimensions
template std::pair<math::Vector<double>, math::Matrix<double>> 
PFApplications::processFovUpdate<double, 2>(
    ParticleFilter& filter,
    const math::VectorDIM<double, 2>& ego_pos,
    const math::VectorDIM<double, 2>& neighbor_pos,
    double fov_beta,
    double fov_Rs,
    double weight_reduction_factor);

template std::pair<math::Vector<double>, math::Matrix<double>> 
PFApplications::processFovUpdate<double, 3>(
    ParticleFilter& filter,
    const math::VectorDIM<double, 3>& ego_pos,
    const math::VectorDIM<double, 3>& neighbor_pos,
    double fov_beta,
    double fov_Rs,
    double weight_reduction_factor);

} // namespace pf
