#pragma once

#include <math/Types.h>
#include <model/DoubleIntegrator.h>
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
// TODO: deprecate this function while usage in mpc_cbf is updated
Vector<double> addRandomNoise(const Vector<double>& xt, double pos_std, double vel_std,
                              unsigned int dim = 3U);

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
