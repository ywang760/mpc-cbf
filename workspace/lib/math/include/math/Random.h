#pragma once

#include <math/Types.h>
#include <random>

namespace math {

/**
 * Add random Gaussian noise to state vector
 * @param xt The state vector
 * @param pos_std Standard deviation for position noise
 * @param vel_std Standard deviation for velocity noise
 * @param dim The dimension of the state space (position components)
 * @return The state vector with added noise
 */
Vector<double> addRandomNoise(
    const Vector<double>& xt, 
    double pos_std, 
    double vel_std, 
    unsigned int dim = 3U);
    // Can potential use a template for dim

} // namespace math
