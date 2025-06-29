#pragma once

#include <math/Types.h>  // 包含 VectorDIM 定义
#include <model/DoubleIntegrator.h>
#include <Eigen/Dense>

namespace pid {

template<typename T>
struct PIDParams {
    T kp, ki, kd, dt;
};

template<typename T>
class PID3D {
public:
    using VectorDIM = math::VectorDIM<T, 3U>;
    using State = model::State<T, 3>;
    PID3D(const PIDParams<T>& params);

    VectorDIM control(State& state, VectorDIM& ref_pos, VectorDIM& ref_vel, VectorDIM& ref_acc);

private:
    T k_p_, k_i_, k_d_, dt_;
    VectorDIM integral_err_;
};

}  // namespace pid

#include "PID3D.cpp"
