#pragma once
#include "PID3D.h"

namespace pid {

template<typename T>
PID3D<T>::PID3D(const PIDParams<T>& params)
    : k_p_(params.kp), k_i_(params.ki), k_d_(params.kd), dt_(params.dt), integral_err_(VectorDIM::Zero()) {}

template<typename T>
typename PID3D<T>::VectorDIM
PID3D<T>::control(State& state,
                  VectorDIM& ref_pos,
                  VectorDIM& ref_vel,
                  VectorDIM& ref_acc) {
    VectorDIM& pos = state.pos_;
    VectorDIM& vel = state.vel_;
    VectorDIM pos_err = ref_pos - pos;
    VectorDIM vel_err = ref_vel - vel;
    integral_err_ += pos_err * dt_;
    VectorDIM control = ref_acc + k_p_ * pos_err + k_i_ * integral_err_ + k_d_ * vel_err;
    return control;
}

}  // namespace pid
