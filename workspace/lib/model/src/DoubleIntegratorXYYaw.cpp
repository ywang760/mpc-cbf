//
// Created by lishuo on 8/16/24.
//

#include <model/DoubleIntegratorXYYaw.h>

namespace model {
    template <typename T>
    DoubleIntegratorXYYaw<T>::DoubleIntegratorXYYaw(T ts) {
        assert(dim_ == 3U);
        // Linear system directly depend on position, velocity and acceleration.
        A_ = Matrix::Zero(6, 6);
        B_ = Matrix::Zero(6, 3);
        A_ << 1, 0, 0, ts, 0, 0,
              0, 1, 0, 0, ts, 0,
              0, 0, 1, 0, 0, ts,
              0, 0, 0, 1, 0, 0,
              0, 0, 0, 0, 1, 0,
              0, 0, 0, 0, 0, 1;

        B_ << 0.5 * pow(ts, 2.0), 0, 0,
              0, 0.5 * pow(ts, 2.0), 0,
              0, 0, 0.5 * pow(ts, 2.0),
              ts, 0, 0,
              0, ts, 0,
              0, 0, ts;
    }

    template class DoubleIntegratorXYYaw<double>;
    template class DoubleIntegratorXYYaw<float>;
} // model