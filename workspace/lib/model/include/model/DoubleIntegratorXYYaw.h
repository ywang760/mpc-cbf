//
// Created by lishuo on 8/16/24.
//

#ifndef MODEL_DOUBLEINTEGRATORXYYAW_H
#define MODEL_DOUBLEINTEGRATORXYYAW_H

#include <model/DoubleIntegrator.h>

namespace model {
/**
     * @brief Double integrator model for x, y position and yaw angle
     *
     * State vector (6 dimensions):
     * - px:   x position
     * - py:   y position
     * - pyaw: yaw angle
     * - vx:   x velocity
     * - vy:   y velocity
     * - vyaw: yaw angular velocity
     *
     * Control input (3 dimensions):
     * - ax:   x acceleration
     * - ay:   y acceleration
     * - ayaw: yaw angular acceleration
     */
template <typename T>
class DoubleIntegratorXYYaw : public DoubleIntegrator<T, 3U> {
  public:
    using DoubleIntegrator<T, 3U>::A_;
    using DoubleIntegrator<T, 3U>::B_;
    using DoubleIntegrator<T, 3U>::dim_;
    using Matrix = math::Matrix<T>;

    /**
         * @brief Constructor
         * @param ts Time step (seconds)
         */
    DoubleIntegratorXYYaw(T ts);
    ~DoubleIntegratorXYYaw() = default;
};

} // namespace model

#endif //MODEL_DOUBLEINTEGRATORXYYAW_H
