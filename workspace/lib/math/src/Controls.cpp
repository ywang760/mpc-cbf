#include <math/Controls.h>

namespace math
{

    template <typename T, unsigned int DIM>
    PID<T, DIM>::PID(const PIDParams<T, DIM> &params)
        : k_p_(params.kp), k_i_(params.ki), k_d_(params.kd), dt_(params.dt), integral_err_(VectorDIM::Zero()) {}

    template <typename T, unsigned int DIM>
    typename PID<T, DIM>::VectorDIM
    PID<T, DIM>::control(State &state,
                    VectorDIM &ref_pos,
                    VectorDIM &ref_vel,
                    VectorDIM &ref_acc)
    {
        VectorDIM &pos = state.pos_;
        VectorDIM &vel = state.vel_;
        VectorDIM pos_err = ref_pos - pos;
        VectorDIM vel_err = ref_vel - vel;
        integral_err_ += pos_err * dt_;
        VectorDIM control = ref_acc + k_p_ * pos_err + k_i_ * integral_err_ + k_d_ * vel_err;
        return control;
    }

    // Explicit template instantiations
    template struct PIDParams<double, 2U>;
    template struct PIDParams<float, 2U>;
    template struct PIDParams<double, 3U>;
    template struct PIDParams<float, 3U>;
    
    template class PID<double, 2U>;
    template class PID<float, 2U>;
    template class PID<double, 3U>;
    template class PID<float, 3U>;

} // namespace math