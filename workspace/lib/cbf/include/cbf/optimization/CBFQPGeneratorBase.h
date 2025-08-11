#ifndef MPC_CBFQPGENERATORBASE_H
#define MPC_CBFQPGENERATORBASE_H

#include <math/Helpers.h>
#include <memory>
#include <qpcpp/Problem.h>
#include <qpcpp/QPOperations.h>

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
         * @brief Constructor for base QP generator
         */
    explicit CBFQPGeneratorBase();

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
         * @param slack_variables Vector of slack variables to apply costs to (defaults to member slack_variables_)
         */
    void addSlackCost(const std::vector<T>& slack_weights,
                      const std::vector<qpcpp::Variable<T>*>& slack_variables = {});

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
    void addControlBoundConstraint(const VectorDIM& u_min, const VectorDIM& u_max);

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
    size_t numDecisionVariables() const {
        return DIM;
    }

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
         * @param slack_variables Vector of slack variables to apply costs to
         */
    void addCostAdditionForSlackVariables(const CostAddition& cost_addition,
                                          const std::vector<qpcpp::Variable<T>*>& slack_variables);

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
    void addLinearConstraintForControlInputWithSlackVariables(
        const LinearConstraint& linear_constraint, const Row& slack_coefficients,
        const std::vector<qpcpp::Variable<T>*>& slack_variable);

    /**
         * @brief Adds bounds on decision variables
         *
         * @param decision_variable_bounds Upper and lower bounds
         */
    void addDecisionVariableBoundsForControlInput(
        const DecisionVariableBounds& decision_variable_bounds);

    qpcpp::Problem<T> problem_;                  ///< The QP problem being constructed
    std::vector<qpcpp::Variable<T>*> variables_; ///< Control input variables
    bool slack_mode_; ///< Whether to use slack variables for constraint relaxation
    std::vector<qpcpp::Variable<T>*>
        slack_variables_; ///< Slack variables for constraint relaxation
};

} // namespace cbf

#endif //MPC_CBFQPGENERATORBASE_H