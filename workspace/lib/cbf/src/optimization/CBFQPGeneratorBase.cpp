//
// Created on 4/10/25.
//

#include <cbf/optimization/CBFQPGeneratorBase.h>

namespace cbf {
template <typename T, unsigned int DIM>
CBFQPGeneratorBase<T, DIM>::CBFQPGeneratorBase(int num_robots, bool slack_mode) {
    // Allocate the decision variables for the control input
    size_t num_decision_variables = DIM;
    for (size_t decision_variable_idx = 0; decision_variable_idx < num_decision_variables;
         ++decision_variable_idx) {
        qpcpp::Variable<T>* variable_ptr = problem_.addVariable();
        variables_.push_back(variable_ptr);
    }

    // Add slack variables to the problem if in slack mode
    // Slack variables allow constraints to be violated by paying a penalty in the cost function
    slack_mode_ = slack_mode;
    if (slack_mode) {
        for (size_t i = 0; i < num_robots; ++i) {
            qpcpp::Variable<T>* variable_ptr = problem().addVariable(
                0, std::numeric_limits<T>::max()); // slack variables are non-negative
            slack_variables_.push_back(variable_ptr);
        }
    }
}

template <typename T, unsigned int DIM>
qpcpp::Problem<T>& CBFQPGeneratorBase<T, DIM>::problem() {
    // Return a reference to the underlying QP problem
    return problem_;
}

template <typename T, unsigned int DIM>
void CBFQPGeneratorBase<T, DIM>::addDesiredControlCost(const VectorDIM& desired_u) {
    // Create quadratic cost for tracking desired control
    Matrix quadratic_term(DIM, DIM);
    quadratic_term.setZero();
    Vector linear_term(DIM);
    linear_term.setZero();

    // Set the quadratic term to identity matrix for a standard quadratic norm
    quadratic_term.setIdentity();

    // Set the linear term as -2 * desired_u to create the expanded form of ||u - desired_u||^2
    linear_term = -2.0 * desired_u;

    // Set the constant term as the squared norm of desired_u
    T constant_term = desired_u.transpose() * desired_u;

    // Create cost addition and add to the QP problem
    CostAddition cost_addition(quadratic_term, linear_term, constant_term);
    addCostAdditionForControlInput(cost_addition);
}

template <typename T, unsigned int DIM>
void CBFQPGeneratorBase<T, DIM>::addSlackCost(const std::vector<T>& slack_weights) {
    // Initialize cost terms for slack variables
    Matrix quadratic_term(slack_weights.size(), slack_weights.size());
    quadratic_term.setZero();
    Vector linear_term(slack_weights.size());
    linear_term.setZero();

    // For slack variables, we only use linear costs to penalize their usage
    for (std::size_t i = 0; i < slack_weights.size(); ++i) {
        linear_term(i) = slack_weights.at(i);
    }

    // Create cost addition and add to the QP problem
    CostAddition cost_addition(quadratic_term, linear_term, 0);
    addCostAdditionForSlackVariables(cost_addition);
}

template <typename T, unsigned int DIM>
void CBFQPGeneratorBase<T, DIM>::addControlBoundConstraint(const VectorDIM& u_min,
                                                           const VectorDIM& u_max) {
    // Initialize vectors for the lower and upper bounds on control inputs
    VectorDIM lower_bounds;
    VectorDIM upper_bounds;

    // Set the individual bounds for each control dimension
    for (size_t d = 0; d < DIM; ++d) {
        lower_bounds(d) = u_min(d);
        upper_bounds(d) = u_max(d);
    }

    // Create decision variable bounds and add to the QP problem
    DecisionVariableBounds decision_variable_bounds(lower_bounds, upper_bounds);
    addDecisionVariableBoundsForControlInput(decision_variable_bounds);
}

template <typename T, unsigned int DIM>
void CBFQPGeneratorBase<T, DIM>::addCostAdditionForControlInput(const CostAddition& cost_addition) {
    // Small value to check if a coefficient is approximately zero
    constexpr T epsilon = std::numeric_limits<T>::epsilon() * T(100.0);
    size_t num_decision_variables = variables_.size();

    // Validate that the cost addition structure matches our variable count
    if (num_decision_variables != cost_addition.linear_term().rows() ||
        num_decision_variables != cost_addition.quadratic_term().rows() ||
        num_decision_variables != cost_addition.quadratic_term().cols()) {
        throw std::runtime_error(
            "CBFQPGeneratorBase::addCostAdditionForControlInput:"
            " number of decision variables of the controlInput does not match the "
            "CostAddition structure");
    }

    // Add the constant term to the cost function
    if (cost_addition.constant() != 0) {
        problem_.cost_function()->add_constant(cost_addition.constant());
    }

    for (size_t i = 0; i < num_decision_variables; ++i) {
        const qpcpp::Variable<T>* var1_ptr = variables_.at(i);

        // Add the linear term if it's not approximately zero
        if (!math::isApproximatelyEqual<T>(cost_addition.linear_term()(i), T(0.0), epsilon)) {
            problem_.cost_function()->addLinearTerm(var1_ptr, cost_addition.linear_term()(i));
        }

        // Add quadratic terms if they're not approximately zero
        for (std::size_t j = 0; j < num_decision_variables; ++j) {
            const qpcpp::Variable<T>* var2_ptr = variables_.at(j);
            if (!math::isApproximatelyEqual<T>(cost_addition.quadratic_term()(i, j), T(0.0),
                                               epsilon)) {
                problem_.cost_function()->addQuadraticTerm(var1_ptr, var2_ptr,
                                                           cost_addition.quadratic_term()(i, j));
            }
        }
    }
}

template <typename T, unsigned int DIM>
void CBFQPGeneratorBase<T, DIM>::addCostAdditionForSlackVariables(
    const CostAddition& cost_addition) {
    // Small value to check if a coefficient is approximately zero
    constexpr T epsilon = std::numeric_limits<T>::epsilon() * T(100.0);
    size_t num_slack_variables = slack_variables_.size();

    // Validate that the cost addition structure matches our slack variable count
    if (num_slack_variables != cost_addition.linear_term().rows() ||
        num_slack_variables != cost_addition.quadratic_term().rows() ||
        num_slack_variables != cost_addition.quadratic_term().cols()) {
        throw std::runtime_error("CBFQPGeneratorBase::addCostAdditionForSlackVariables:"
                                 " number of slack variables does not match the "
                                 "CostAddition structure");
    }

    // Add the constant term to the cost function
    if (cost_addition.constant() != 0) {
        problem().cost_function()->add_constant(cost_addition.constant());
    }

    for (size_t i = 0; i < num_slack_variables; ++i) {
        const qpcpp::Variable<T>* var1_ptr = slack_variables_.at(i);

        // Add the linear term if it's not approximately zero
        if (!math::isApproximatelyEqual<T>(cost_addition.linear_term()(i), T(0.0), epsilon)) {
            problem().cost_function()->addLinearTerm(var1_ptr, cost_addition.linear_term()(i));
        }

        // Add quadratic terms if they're not approximately zero
        for (std::size_t j = 0; j < num_slack_variables; ++j) {
            const qpcpp::Variable<T>* var2_ptr = slack_variables_.at(j);
            if (!math::isApproximatelyEqual<T>(cost_addition.quadratic_term()(i, j), T(0.0),
                                               epsilon)) {
                problem().cost_function()->addQuadraticTerm(var1_ptr, var2_ptr,
                                                            cost_addition.quadratic_term()(i, j));
            }
        }
    }
}

template <typename T, unsigned int DIM>
void CBFQPGeneratorBase<T, DIM>::addLinearConstraintForControlInput(
    const LinearConstraint& linear_constraint) {
    size_t num_decision_variables = variables_.size();

    // Validate that the constraint structure matches our variable count
    if (num_decision_variables != linear_constraint.coefficients().cols()) {
        throw std::runtime_error(
            "CBFQPGeneratorBase::"
            "addLinearConstraintForControlInput:"
            " number of decision variables of the controlInput does not match the "
            "LinearConstraint structure");
    }

    // Initialize a new linear constraint with lower and upper bounds
    qpcpp::LinearConstraint<T>* qpcpp_linear_constraint = problem_.addLinearConstraint(
        linear_constraint.lower_bound(), linear_constraint.upper_bound());

    // Set coefficients for each decision variable in this constraint
    for (std::size_t decision_variable_idx = 0; decision_variable_idx < num_decision_variables;
         ++decision_variable_idx) {
        const qpcpp::Variable<T>* var_ptr = variables_.at(decision_variable_idx);
        qpcpp_linear_constraint->setCoefficient(
            var_ptr, linear_constraint.coefficients()(decision_variable_idx));
    }
}

template <typename T, unsigned int DIM>
void CBFQPGeneratorBase<T, DIM>::addLinearConstraintForControlInputWithSlackVariables(
    const LinearConstraint& linear_constraint, const Row& slack_coefficients) {
    size_t num_decision_variables = variables_.size();
    size_t num_slack_variables = slack_variables_.size();

    // Validate that the constraint and slack coefficient structures match our variable counts
    if (num_decision_variables != linear_constraint.coefficients().cols()) {
        throw std::runtime_error(
            "CBFQPGeneratorBase::"
            "addLinearConstraintForControlInputWithSlackVariables:"
            " number of decision variables of the controlInput does not match the "
            "LinearConstraint structure");
    }
    if (num_slack_variables != slack_coefficients.size()) {
        throw std::runtime_error("CBFQPGeneratorBase::"
                                 "addLinearConstraintForControlInputWithSlackVariables:"
                                 " number of slack variables does not match the "
                                 "slack_coefficients structure");
    }

    // Initialize a new linear constraint with lower and upper bounds
    qpcpp::LinearConstraint<T>* qpcpp_linear_constraint = problem_.addLinearConstraint(
        linear_constraint.lower_bound(), linear_constraint.upper_bound());

    // Set coefficients for each decision variable in this constraint
    for (std::size_t decision_variable_idx = 0; decision_variable_idx < num_decision_variables;
         ++decision_variable_idx) {
        const qpcpp::Variable<T>* var_ptr = variables_.at(decision_variable_idx);
        qpcpp_linear_constraint->setCoefficient(
            var_ptr, linear_constraint.coefficients()(decision_variable_idx));
    }

    // Set coefficients for each slack variable in this constraint
    for (std::size_t slack_variable_idx = 0; slack_variable_idx < num_slack_variables;
         ++slack_variable_idx) {
        const qpcpp::Variable<T>* var_ptr = slack_variables_.at(slack_variable_idx);
        qpcpp_linear_constraint->setCoefficient(var_ptr, slack_coefficients(slack_variable_idx));
    }
}

template <typename T, unsigned int DIM>
void CBFQPGeneratorBase<T, DIM>::addDecisionVariableBoundsForControlInput(
    const DecisionVariableBounds& decision_variable_bounds) {
    size_t num_decision_variables = variables_.size();

    // Validate that the bounds structure matches our variable count
    if (num_decision_variables != decision_variable_bounds.lower_bounds().rows() ||
        num_decision_variables != decision_variable_bounds.upper_bounds().rows()) {
        throw std::invalid_argument(
            "CBFQPGeneratorBase::"
            "addDecisionVariableBoundsForControlInput:"
            " number of decision variables of the controlInput does not match the "
            "DecisionVariablesBounds structure");
    }

    // Set the bounds for each decision variable
    for (int decision_variable_idx = 0; decision_variable_idx < DIM; ++decision_variable_idx) {
        qpcpp::Variable<T>* var_ptr = variables_.at(decision_variable_idx);

        // Take the more restrictive of the existing and new lower bounds
        var_ptr->set_min(std::max(var_ptr->min(),
                                  decision_variable_bounds.lower_bounds()(decision_variable_idx)));

        // Take the more restrictive of the existing and new upper bounds
        var_ptr->set_max(std::min(var_ptr->max(),
                                  decision_variable_bounds.upper_bounds()(decision_variable_idx)));
    }
}

template <typename T, unsigned int DIM>
typename CBFQPGeneratorBase<T, DIM>::VectorDIM
CBFQPGeneratorBase<T, DIM>::generatorCBFControlInput() const {
    assert(DIM == numDecisionVariables());

    // Extract the solution values from the QP solver into a vector
    VectorDIM solution;
    for (size_t d = 0; d < DIM; ++d) {
        solution(d) = variables_.at(d)->solution_value();
    }
    return solution;
}

// Explicit template instantiation
template class CBFQPGeneratorBase<double, 3U>;

} // namespace cbf