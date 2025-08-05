#pragma once

#include <Eigen/Dense>
#include <cmath>
#include <math/Types.h>
#include <model/DoubleIntegrator.h>

namespace math {

/**
     * Critically damped spring control law for smooth trajectory following
     * @param current_state The current state of the system
     * @param target The target position
     * @param spring_constant The spring constant (stiffness)
     * @return Control input (acceleration)
     */
template <typename T = double, unsigned int DIM = 3>
VectorDIM<T, DIM> criticallyDampedSpringControl(const model::State<T, DIM>& current_state,
                                                const VectorDIM<T, DIM>& target,
                                                const double spring_constant) {

    VectorDIM<T, DIM> Fs;
    Fs = spring_constant * (target - current_state.pos_);
    VectorDIM<T, DIM> Fd;
    Fd = -current_state.vel_ * 2 * sqrt(spring_constant);
    return Fs + Fd;
}

template <typename T = double, unsigned int DIM = 3>
struct PIDParams {
    T kp, ki, kd, dt;
};
template <typename T = double, unsigned int DIM = 3>
class PID {
  public:
    using VectorDIM = math::VectorDIM<T, DIM>;
    using State = model::State<T, DIM>;
    PID(const PIDParams<T, DIM>& params);

    VectorDIM control(State& state, VectorDIM& ref_pos, VectorDIM& ref_vel, VectorDIM& ref_acc);

  private:
    T k_p_, k_i_, k_d_, dt_;
    VectorDIM integral_err_;
};

} // namespace math
