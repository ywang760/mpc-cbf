//
// Created by lishuo on 8/31/24.
//

#ifndef CBF_OPTIMIZATION_CBFQPOPERATIONS_H
#define CBF_OPTIMIZATION_CBFQPOPERATIONS_H

#include <qpcpp/QPOperations.h>
#include <cbf/detail/cbf.h>

namespace cbf {
    /**
     * @brief Operations for Control Barrier Function Quadratic Programming
     *
     * This class provides the operations needed to formulate a Quadratic Program (QP)
     * for maintaining safety and visibility constraints through Control Barrier Functions (CBFs).
     * It handles the formulation of cost functions and constraints based on the underlying CBF.
     *
     * @tparam T The numeric type (typically double)
     * @tparam DIM Dimension of the control input
     */
    template <typename T, unsigned int DIM>
    class CBFQPOperations {
    public:
        using CostAddition = typename qpcpp::QPOperations<T>::CostAddition;
        using LinearConstraint = typename qpcpp::QPOperations<T>::LinearConstraint;
        using DecisionVariableBounds = typename qpcpp::QPOperations<T>::DecisionVariableBounds;
        using VectorDIM = math::VectorDIM<T, DIM>;
        using Vector = math::Vector<T>;
        using Matrix = math::Matrix<T>;

        /**
         * @brief Constructor
         *
         * @param cbf Shared pointer to the Field of View CBF implementation
         */
        CBFQPOperations(std::shared_ptr<FovCBF> cbf);

        /**
         * @brief Returns the number of decision variables in the QP
         *
         * @return The dimension of the control input
         */
        size_t numDecisionVariables() const {return DIM;}

        /**
         * @brief Generates cost terms to track a desired control input
         *
         * Creates a quadratic cost function to minimize the squared distance between
         * the actual and desired control input: ||u - desired_u||^2
         *
         * @param desired_u The desired control input to track
         * @return Cost function components (quadratic, linear, and constant terms)
         */
        CostAddition desiredControlCost(const VectorDIM& desired_u);

        /**
         * @brief Generates cost terms for slack variables
         *
         * Creates a cost function to penalize constraint violations through slack variables
         *
         * @param slack_weights Weights for each slack variable
         * @return Cost function components for slack variables
         */
        CostAddition slackCost(const std::vector<T>& slack_weights);

        /**
         * @brief Generates safety constraint between agents
         *
         * Creates a constraint to prevent collisions between robots
         *
         * @param state Current state of the robot
         * @param target_state State of the target/neighboring robot
         * @return Linear constraint representation for safety
         */
        LinearConstraint safetyConstraint(const Vector& state, const Vector& target_state);

        /**
         * @brief Generates field-of-view left border constraint
         *
         * Creates a constraint to keep the target within the left border of the field of view
         *
         * @param state Current state of the robot
         * @param target_state State of the target/neighboring robot
         * @return Linear constraint representation for left border
         */
        LinearConstraint leftBorderConstraint(const Vector& state, const Vector& target_state);

        /**
         * @brief Generates field-of-view right border constraint
         *
         * Creates a constraint to keep the target within the right border of the field of view
         *
         * @param state Current state of the robot
         * @param target_state State of the target/neighboring robot
         * @return Linear constraint representation for right border
         */
        LinearConstraint rightBorderConstraint(const Vector& state, const Vector& target_state);

        /**
         * @brief Generates range constraint
         *
         * Creates a constraint to keep the target within the sensing range
         *
         * @param state Current state of the robot
         * @param target_state State of the target/neighboring robot
         * @return Linear constraint representation for range
         */
        LinearConstraint rangeConstraint(const Vector& state, const Vector& target_state);

        /**
         * @brief Generates minimum velocity constraints
         *
         * Creates constraints to enforce minimum velocity requirements
         *
         * @param state Current state of the robot
         * @return Vector of linear constraints for minimum velocities
         */
        std::vector<LinearConstraint> minVelConstraints(const Vector& state);

        /**
         * @brief Generates maximum velocity constraints
         *
         * Creates constraints to enforce maximum velocity limits
         *
         * @param state Current state of the robot
         * @return Vector of linear constraints for maximum velocities
         */
        std::vector<LinearConstraint> maxVelConstraints(const Vector& state);

        /**
         * @brief Generates left border constraint with slack variable
         *
         * Creates a relaxed constraint for the left field-of-view border
         *
         * @param state Current state of the robot
         * @param target_state State of the target/neighboring robot
         * @param slack Slack variable value for constraint relaxation
         * @return Relaxed linear constraint representation
         */
        LinearConstraint leftBorderConstraintWithSlackVar(const Vector& state, const Vector& target_state, const T &slack);

        /**
         * @brief Generates right border constraint with slack variable
         *
         * Creates a relaxed constraint for the right field-of-view border
         *
         * @param state Current state of the robot
         * @param target_state State of the target/neighboring robot
         * @param slack Slack variable value for constraint relaxation
         * @return Relaxed linear constraint representation
         */
        LinearConstraint rightBorderConstraintWithSlackVar(const Vector& state, const Vector& target_state, const T &slack);

        /**
         * @brief Generates range constraint with slack variable
         *
         * Creates a relaxed constraint for sensing range
         *
         * @param state Current state of the robot
         * @param target_state State of the target/neighboring robot
         * @param slack Slack variable value for constraint relaxation
         * @return Relaxed linear constraint representation
         */
        LinearConstraint rangeConstraintWithSlackVar(const Vector& state, const Vector& target_state, const T &slack);

        /**
         * @brief Generates bounds for control inputs
         *
         * Creates upper and lower bound constraints for control inputs
         *
         * @param u_min Minimum allowed values for control inputs
         * @param u_max Maximum allowed values for control inputs
         * @return Bound constraints for decision variables
         */
        DecisionVariableBounds controlBoundConstraint(const VectorDIM& u_min, const VectorDIM& u_max);

    private:
        std::shared_ptr<FovCBF> cbf_; ///< Pointer to the Field of View CBF implementation
    };

} // cbf

#endif //CBF_OPTIMIZATION_CBFQPOPERATIONS_H
