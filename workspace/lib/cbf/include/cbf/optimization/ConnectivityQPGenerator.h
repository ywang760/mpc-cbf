#ifndef MPC_CONNECTIVITYQPGENERATOR_H
#define MPC_CONNECTIVITYQPGENERATOR_H

#include <cbf/optimization/CBFQPGeneratorBase.h>
#include <cbf/detail/ConnectivityCBF.h>

namespace cbf {
    /**
     * @brief Connectivity Quadratic Program Generator
     * 
     * This class implements the connectivity-specific constraint methods for CBF QP generation.
     * It handles connectivity constraints to maintain network connectivity between agents.
     * 
     * @tparam T The numeric type (typically double)
     * @tparam DIM Dimension of the control input
     */
    template <typename T, unsigned int DIM>
    class ConnectivityQPGenerator : public CBFQPGeneratorBase<T, DIM> {
    public:
        using typename CBFQPGeneratorBase<T, DIM>::VectorDIM;
        using typename CBFQPGeneratorBase<T, DIM>::Vector;
        using typename CBFQPGeneratorBase<T, DIM>::Row;
        using typename CBFQPGeneratorBase<T, DIM>::LinearConstraint;
        using typename CBFQPGeneratorBase<T, DIM>::Matrix;


        /**
         * @brief Constructor that takes a pointer to the ConnectivityCBF instance
         *
         * @param cbf Shared pointer to the Connectivity CBF implementation
         * @param num_robots Number of robots (used for slack variables)
         * @param slack_mode If true, add slack variables to allow constraint relaxation
         */
        explicit ConnectivityQPGenerator(std::shared_ptr<ConnectivityCBF> cbf, int num_robots = 0, bool slack_mode = false);

        /**
         * @brief Adds connectivity constraint between agents to maintain network connectivity
         *
         * @param state Current state of the robot
         * @param robot_states Matrix of states for all robots in the network
         * @param use_slack Whether to use slack variables for this constraint
         */
        void addConnConstraint(const Vector &state, const Eigen::MatrixXd &robot_states, size_t self_idx,
                               bool use_slack = false);

        /**
         * @brief Adds safety constraint between agents to prevent collisions
         *
         * @param state Current state of the robot
         * @param neighbor_state State of the target/neighboring robot
         * @param use_slack Whether to use slack variables for this constraint
         * @param slack_idx Index of slack variable to use (if use_slack is true)
         */
        void addSafetyConstraint(const Vector& state, const Vector& neighbor_state, 
                               bool use_slack = false, std::size_t slack_idx = 0) override;

        /**
         * @brief Adds minimum velocity constraints
         *
         * @param state Current state of the robot
         */
        void addMinVelConstraints(const Vector& state) override;

        /**
         * @brief Adds maximum velocity constraints
         *
         * @param state Current state of the robot
         */
        void addMaxVelConstraints(const Vector& state) override;

        /**
        * @brief Accessor for the internal ConnectivityCBF instance
        * @return Shared pointer to the ConnectivityCBF object
        */
        const std::shared_ptr<ConnectivityCBF>& getCBF() const {
            return cbf_;
        }


    private:
        std::shared_ptr<ConnectivityCBF> cbf_; ///< Pointer to the Connectivity CBF implementation

        // Using methods from the base class
        using CBFQPGeneratorBase<T, DIM>::addLinearConstraintForControlInput;
        using CBFQPGeneratorBase<T, DIM>::addLinearConstraintForControlInputWithSlackVariables;
    };

} // cbf

#endif //MPC_CONNECTIVITYQPGENERATOR_H