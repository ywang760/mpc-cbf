#ifndef MPC_CBFQPGENERATORBASE_H
#define MPC_CBFQPGENERATORBASE_H

#include <qpcpp/QPOperations.h>
#include <math/Helpers.h>
#include <qpcpp/Problem.h>
#include <memory>

namespace cbf {
    /**
     * @brief Base Quadratic Program Generator for Control Barrier Functions
     *
     * This class provides the common infrastructure for QP generation in different CBF applications.
     * Derived classes will implement specific constraint types for different applications (FOV, connectivity, etc.).
     *
     * @tparam T The numeric type (typically double)
     * @tparam DIM Dimension of the control input
     */
    template <typename T, unsigned int DIM>
    class CBFQPGeneratorBase {
    public:
        using CostAddition = typename qpcpp::QPOperations<T>::CostAddition;
        using LinearConstraint = typename qpcpp::QPOperations<T>::LinearConstraint;
        using DecisionVariableBounds = typename qpcpp::QPOperations<T>::DecisionVariableBounds;
        using VectorDIM = math::VectorDIM<T, DIM>;
        using Vector = math::Vector<T>;
        using Row = math::Row<T>;
        using Matrix = math::Matrix<T>;

        /**
         * @brief Constructor with number of robots and slack mode
         *
         * @param num_robots Number of robots (used for slack variables)
         * @param slack_mode If true, add slack variables to allow constraint relaxation
         */
        explicit CBFQPGeneratorBase(int num_robots = 0, bool slack_mode = false);

        /**
         * @brief Virtual destructor
         */
        virtual ~CBFQPGeneratorBase() = default;

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
         * This is a generic interface that will be implemented by derived classes
         * based on their specific barrier function formulations.
         *
         * @param state Current state of the robot
         * @param neighbor_state State of the target/neighboring robot
         * @param use_slack Whether to use slack variables for this constraint
         * @param slack_idx Index of slack variable to use (if use_slack is true)
         */
        virtual void addSafetyConstraint(const Vector& state, const Vector& neighbor_state, 
                                bool use_slack = false, std::size_t slack_idx = 0) = 0;
                                
        /**
         * @brief Adds minimum velocity constraints
         *
         * This is a generic interface that will be implemented by derived classes
         * based on their specific barrier function formulations.
         *
         * @param state Current state of the robot
         */
        virtual void addMinVelConstraints(const Vector& state) = 0;

        /**
         * @brief Adds maximum velocity constraints
         *
         * This is a generic interface that will be implemented by derived classes
         * based on their specific barrier function formulations.
         *
         * @param state Current state of the robot
         */
        virtual void addMaxVelConstraints(const Vector& state) = 0;

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

        /**
         * @brief Returns the number of decision variables in the QP
         *
         * @return The dimension of the control input
         */
        size_t numDecisionVariables() const { return DIM; }

    protected:
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

        qpcpp::Problem<T> problem_;                         ///< The QP problem being constructed
        std::vector<qpcpp::Variable<T> *> variables_;       ///< Control input variables
        std::vector<qpcpp::Variable<T> *> slack_variables_; ///< Slack variables for constraint relaxation
    };

} // cbf

#endif //MPC_CBFQPGENERATORBASE_H