#pragma once

#include <math/Types.h>
#include <model/DoubleIntegrator.h>
#include <cmath>

namespace math {

/**
 * Critically damped spring control law for smooth trajectory following
 * @param current_state The current state of the system
 * @param target The target position
 * @param spring_constant The spring constant (stiffness)
 * @return Control input (acceleration)
 */
template<typename T = double, unsigned int DIM = 3>
VectorDIM<T, DIM> criticallyDampedSpringControl(
    const model::State<T, DIM>& current_state, 
    const VectorDIM<T, DIM>& target, 
    const double spring_constant) {
    
    VectorDIM<T, DIM> Fs;
    Fs = spring_constant * (target - current_state.pos_);
    VectorDIM<T, DIM> Fd;
    Fd = -current_state.vel_ * 2 * sqrt(spring_constant);
    return Fs + Fd;
}

} // namespace math
