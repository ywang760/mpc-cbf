#pragma once

#include <math/Types.h>
#include <model/DoubleIntegrator.h>
#include <random>

namespace math {
    /**
     * Template function to add random Gaussian noise to state
     * @tparam T Scalar type (double or float)
     * @tparam DIM Dimension of the state space
     * @param state The state object
     * @param pos_std Standard deviation for position noise
     * @param vel_std Standard deviation for velocity noise
     * @return The state with added noise
     */
template <typename T, unsigned int DIM>
model::State<T, DIM> addRandomNoise(const model::State<T, DIM>& state, T pos_std, T vel_std);

} // namespace math
