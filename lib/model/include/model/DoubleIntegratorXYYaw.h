//
// Created by lishuo on 8/16/24.
//

#ifndef MODEL_DOUBLEINTEGRATORXYYAW_H
#define MODEL_DOUBLEINTEGRATORXYYAW_H

#include <model/DoubleIntegrator.h>

namespace model {
    template <typename T>
    class DoubleIntegratorXYYaw : public DoubleIntegrator<T, 3U> {
    public:
        using DoubleIntegrator<T, 3U>::A_;
        using DoubleIntegrator<T, 3U>::B_;
        using DoubleIntegrator<T, 3U>::dim_;
        using Matrix = math::Matrix<T>;

        DoubleIntegratorXYYaw(T ts);
        ~DoubleIntegratorXYYaw()=default;

    };

} // model

#endif //MODEL_DOUBLEINTEGRATORXYYAW_H
