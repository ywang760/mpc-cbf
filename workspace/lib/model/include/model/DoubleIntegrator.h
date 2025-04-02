//
// Created by lishuo on 8/15/24.
//

#ifndef MODEL_DOUBLEINTEGRATOR_H
#define MODEL_DOUBLEINTEGRATOR_H

#include <math/Types.h>

namespace model {
    /**
     * @brief State struct for double integrator dynamics
     *
     * Stores position and velocity vectors for a system with DIM dimensions
     *
     * @tparam T Numeric type (float, double, etc.)
     * @tparam DIM Number of dimensions in the state space
     */
    template <typename T, unsigned int DIM>
    struct State {
        using VectorDIM = math::VectorDIM<T, DIM>;
        VectorDIM pos_, vel_; // Position and velocity vectors
    };

    /**
     * @brief Propagator for state transitions
     *
     * Contains matrices used for state propagation in discrete-time systems
     *
     * @tparam T Numeric type (float, double, etc.)
     */
    template <typename T>
    struct StatePropagator {
        using Matrix = math::Matrix<T>;
        Matrix pos_, vel_; // Matrices for position and velocity propagation
    };

    /**
     * @brief Double integrator dynamics model
     *
     * Implements a double integrator system with dynamics:
     * ẋ = Ax + Bu, where x = [pos, vel]^T
     *
     * @tparam T Numeric type (float, double, etc.)
     * @tparam DIM Number of dimensions in the state space
     */
    template <typename T, unsigned int DIM>
    class DoubleIntegrator {
    public:
        using State = model::State<T, DIM>;
        using Matrix = math::Matrix<T>;
        using Vector = math::Vector<T>;

        /** @brief Default constructor */
        DoubleIntegrator()= default;

        /** @brief Default destructor */
        ~DoubleIntegrator()= default;

        /**
         * @brief Get the A0 matrix for K steps ahead prediction
         * @param K Number of time steps
         * @return StatePropagator containing propagation matrices
         */
        StatePropagator<T> get_A0(int K);

        /**
         * @brief Get the lambda matrices for input effect over K steps
         * @param K Number of time steps
         * @return StatePropagator containing input effect matrices
         */
        StatePropagator<T> get_lambda(int K);

        /**
         * @brief Apply control input to current state
         * @param state Current state
         * @param u Control input vector
         * @return Next state after applying input
         */
        State applyInput(const State& state, const Vector& u);

    protected:
        // Initialize A matrix for double integrator system
        // [I  dt*I]
        // [0    I ]
        Matrix A_ = Matrix::Identity(2 * DIM, 2 * DIM);

        // Initialize B matrix for double integrator system
        // [0.5*dt^2*I]
        // [dt*I      ]
        Matrix B_ = Matrix::Zero(2 * DIM, DIM);

        int dim_ = DIM; // Dimension of the state space
    };

} // model


#endif //MODEL_DOUBLEINTEGRATOR_H
