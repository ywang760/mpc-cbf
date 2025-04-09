//
// Created by lishuo on 8/31/24.
//

#ifndef MPC_CBFQPGENERATOR_H
#define MPC_CBFQPGENERATOR_H

#include <cbf/optimization/CBFQPOperations.h>
#include <math/Helpers.h>
#include <qpcpp/Problem.h>

namespace cbf {
    /**
     * @brief Quadratic Program Generator for Control Barrier Functions
     *
     * This class is responsible for generating a quadratic program (QP) for control barrier functions (CBF).
     * It handles cost additions, constraint formulations, and provides mechanisms to solve the QP to get
     * a control input that satisfies CBF constraints.
     *
     * @tparam T The numeric type (typically double)
     * @tparam DIM Dimension of the control input
     */
    template <typename T, unsigned int DIM>
    class CBFQPGenerator {
    public:
        using CBFQPOperations = cbf::CBFQPOperations<T, DIM>;
        using CostAddition = typename cbf::CBFQPOperations<T, DIM>::CostAddition;
        using LinearConstraint = typename cbf::CBFQPOperations<T, DIM>::LinearConstraint;
        using DecisionVariableBounds = typename cbf::CBFQPOperations<T, DIM>::DecisionVariableBounds;
        using VectorDIM = math::VectorDIM<T, DIM>;
        using Vector = math::Vector<T>;
        using Row = math::Row<T>;

        /**
         * @brief Initializes the QP with CBF operations and optionally adds slack variables
         *
         * @param cbf_operations_ptr Pointer to the CBF operations implementation
         * @param num_neighbors Number of neighboring agents (used for slack variables)
         * @param slack_mode If true, add slack variables to allow constraint relaxation
         */
        void addCBFOperations(std::unique_ptr<CBFQPOperations> cbf_operations_ptr, int num_neighbors=0, bool slack_mode=false);

        /**
         * @brief Returns a reference to the underlying QP problem
         */
        qpcpp::Problem<T>& problem();

        /**
         * @brief Adds cost term to minimize difference from desired control input
         *
         * @param desired_u The desired control input to track
         */
        void addDesiredControlCost(const VectorDIM& desired_u);

        /**
         * @brief Adds cost terms for slack variables to penalize constraint violations
         *
         * @param slack_weights Weights for each slack variable in the cost function
         */
        void addSlackCost(const std::vector<T>& slack_weights);

        /**
         * @brief Adds safety constraint between agents to prevent collisions
         *
         * @param state Current state of the robot
         * @param target_state State of the target/neighboring robot
         */
        void addSafetyConstraint(const Vector& state, const Vector& target_state);

        /**
         * @brief Adds field-of-view left border constraint
         *
         * @param state Current state of the robot
         * @param target_state State of the target/neighboring robot
         */
        void addLeftBorderConstraint(const Vector& state, const Vector& target_state);

        /**
         * @brief Adds field-of-view right border constraint
         *
         * @param state Current state of the robot
         * @param target_state State of the target/neighboring robot
         */
        void addRightBorderConstraint(const Vector& state, const Vector& target_state);

        /**
         * @brief Adds range constraint to maintain target within sensor range
         *
         * @param state Current state of the robot
         * @param target_state State of the target/neighboring robot
         */
        void addRangeConstraint(const Vector& state, const Vector& target_state);

        /**
         * @brief Adds minimum velocity constraints
         *
         * @param state Current state of the robot
         */
        void addMinVelConstraints(const Vector& state);

        /**
         * @brief Adds maximum velocity constraints
         *
         * @param state Current state of the robot
         */
        void addMaxVelConstraints(const Vector& state);

        /**
         * @brief Adds safety constraint with slack variable to allow relaxation
         *
         * @param state Current state of the robot
         * @param target_state State of the target/neighboring robot
         * @param neighbor_idx Index of the neighbor for slack variable assignment
         */
        void addSafetyCBFConstraintWithSlackVariables(const Vector& state, const Vector& target_state, std::size_t neighbor_idx);

        /**
         * @brief Adds left border constraint with slack variable
         *
         * @param state Current state of the robot
         * @param target_state State of the target/neighboring robot
         * @param neighbor_idx Index of the neighbor for slack variable assignment
         */
        void addLeftBorderConstraintWithSlackVariables(const Vector& state, const Vector& target_state, std::size_t neighbor_idx);

        /**
         * @brief Adds right border constraint with slack variable
         *
         * @param state Current state of the robot
         * @param target_state State of the target/neighboring robot
         * @param neighbor_idx Index of the neighbor for slack variable assignment
         */
        void addRightBorderConstraintWithSlackVariables(const Vector& state, const Vector& target_state, std::size_t neighbor_idx);

        /**
         * @brief Adds range constraint with slack variable
         *
         * @param state Current state of the robot
         * @param target_state State of the target/neighboring robot
         * @param neighbor_idx Index of the neighbor for slack variable assignment
         */
        void addRangeConstraintWithSlackVariables(const Vector& state, const Vector& target_state, std::size_t neighbor_idx);

        /**
         * @brief Adds bounds on control inputs
         *
         * @param u_min Minimum values for each control input dimension
         * @param u_max Maximum values for each control input dimension
         */
        void addControlBoundConstraint(const VectorDIM &u_min, const VectorDIM &u_max);

        /**
         * @brief Retrieves the optimal control input after solving the QP
         *
         * @return The optimal control input vector
         */
        VectorDIM generatorCBFControlInput() const;

    private:
        /**
         * @brief Adds cost terms for control input variables
         *
         * @param cost_addition The cost terms to add (quadratic, linear, constant)
         */
        void addCostAdditionForControlInput(const CostAddition& cost_addition);

        /**
         * @brief Adds cost terms for slack variables
         *
         * @param cost_addition The cost terms to add (quadratic, linear, constant)
         */
        void addCostAdditionForSlackVariables(const CostAddition& cost_addition);

        /**
         * @brief Adds linear constraint for control input variables
         *
         * @param linear_constraint The constraint to add
         */
        void addLinearConstraintForControlInput(const LinearConstraint& linear_constraint);

        /**
         * @brief Adds linear constraint involving both control input and slack variables
         *
         * @param linear_constraint The constraint for control input variables
         * @param slack_coefficients Coefficients for slack variables in the constraint
         */
        void addLinearConstraintForControlInputWithSlackVariables(const LinearConstraint& linear_constraint,
                                                                  const Row &slack_coefficients);

        /**
         * @brief Adds bounds on decision variables
         *
         * @param decision_variable_bounds Upper and lower bounds
         */
        void addDecisionVariableBoundsForControlInput(const DecisionVariableBounds& decision_variable_bounds);

        std::unique_ptr<CBFQPOperations> cbf_operations_ptr_; ///< Operations related to CBF constraints

        qpcpp::Problem<T> problem_;                         ///< The QP problem being constructed
        std::vector<qpcpp::Variable<T> *> variables_;       ///< Control input variables
        std::vector<qpcpp::Variable<T> *> slack_variables_; ///< Slack variables for constraint relaxation
    };

} // cbf

#endif //MPC_CBFQPGENERATOR_H
