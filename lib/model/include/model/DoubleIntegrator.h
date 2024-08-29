//
// Created by lishuo on 8/15/24.
//

#ifndef MODEL_DOUBLEINTEGRATOR_H
#define MODEL_DOUBLEINTEGRATOR_H

#include <math/Types.h>

namespace model {
    template <typename T, unsigned int DIM>
    struct State {
        using VectorDIM = math::VectorDIM<T, DIM>;
        VectorDIM pos_, vel_;
    };

    template <typename T>
    struct StatePropagator {
        using Matrix = math::Matrix<T>;
        Matrix pos_, vel_;
    };

    template <typename T, unsigned int DIM>
    class DoubleIntegrator {
    public:
        using State = model::State<T, DIM>;
        using Matrix = math::Matrix<T>;
        using Vector = math::Vector<T>;

        DoubleIntegrator()= default;
        ~DoubleIntegrator()= default;

        StatePropagator<T> get_A0(int K);
        StatePropagator<T> get_lambda(int K);

        State applyInput(const State& state, const Vector& u);

    protected:
        Matrix A_;
        Matrix B_;
        int dim_ = DIM;
    };

} // model


#endif //MODEL_DOUBLEINTEGRATOR_H
