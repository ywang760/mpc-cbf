//
// Created by lishuo on 8/31/24.
//

#include <cbf/optimization/CBFQPGenerator.h>

namespace cbf {
    template <typename T, unsigned int DIM>
    void
    CBFQPGenerator<T, DIM>::addCBFOperations(std::unique_ptr<CBFQPOperations> cbf_operations_ptr, int num_neighbors, bool slack_mode) {
        // Store the operations pointer for later use
        cbf_operations_ptr_ = std::move(cbf_operations_ptr);

        // Allocate the decision variables for the control input
        size_t num_decision_variables_of_control_input = cbf_operations_ptr_->numDecisionVariables();
        for (size_t decision_variable_idx = 0;
             decision_variable_idx < num_decision_variables_of_control_input;
             ++decision_variable_idx) {
            qpcpp::Variable<T>* variable_ptr = problem_.addVariable();
            variables_.push_back(variable_ptr);
        }

        // Add slack variables to the problem if in slack mode
        // Slack variables allow constraints to be violated by paying a penalty in the cost function
        if (slack_mode) {
            for (size_t i = 0; i < num_neighbors; ++i) {
                qpcpp::Variable<T>* variable_ptr = problem().addVariable(0, std::numeric_limits<T>::max()); // slack variables are non-negative
                slack_variables_.push_back(variable_ptr);
            }
        }
    }

    template <typename T, unsigned int DIM>
    qpcpp::Problem<T>& CBFQPGenerator<T, DIM>::problem() {
        // Return a reference to the underlying QP problem
        return problem_;
    }

    template <typename T, unsigned int DIM>
    void CBFQPGenerator<T, DIM>::addDesiredControlCost(const VectorDIM& desired_u) {
        // Get the cost terms for tracking the desired control input
        CostAddition cost_addition = cbf_operations_ptr_->desiredControlCost(desired_u);
        // Add these cost terms to the QP problem
        addCostAdditionForControlInput(cost_addition);
    }

    template <typename T, unsigned int DIM>
    void CBFQPGenerator<T, DIM>::addSlackCost(const std::vector<T> &slack_weights) {
        // Get the cost terms for the slack variables based on the provided weights
        CostAddition cost_addition = cbf_operations_ptr_->slackCost(slack_weights);
        // Add these cost terms to the QP problem
        addCostAdditionForSlackVariables(cost_addition);
    }

    template <typename T, unsigned int DIM>
    void CBFQPGenerator<T, DIM>::addSafetyConstraint(const Vector &state,
                                                     const Vector &target_state) {
        // Get the safety constraint based on current states
        LinearConstraint linear_constraint = cbf_operations_ptr_->safetyConstraint(state, target_state);
        // Add this constraint to the QP problem
        addLinearConstraintForControlInput(linear_constraint);
    }

    template <typename T, unsigned int DIM>
    void CBFQPGenerator<T, DIM>::addLeftBorderConstraint(const Vector &state,
                                                         const Vector &target_state) {
        // Get the left field-of-view border constraint
        LinearConstraint linear_constraint = cbf_operations_ptr_->leftBorderConstraint(state, target_state);
        // Add this constraint to the QP problem
        addLinearConstraintForControlInput(linear_constraint);
    }

    template <typename T, unsigned int DIM>
    void CBFQPGenerator<T, DIM>::addRightBorderConstraint(const Vector &state,
                                                          const Vector &target_state) {
        // Get the right field-of-view border constraint
        LinearConstraint linear_constraint = cbf_operations_ptr_->rightBorderConstraint(state, target_state);
        // Add this constraint to the QP problem
        addLinearConstraintForControlInput(linear_constraint);
    }

    template <typename T, unsigned int DIM>
    void CBFQPGenerator<T, DIM>::addRangeConstraint(const Vector &state,
                                                    const Vector &target_state) {
        // Get the sensing range constraint
        LinearConstraint linear_constraint = cbf_operations_ptr_->rangeConstraint(state, target_state);
        // Add this constraint to the QP problem
        addLinearConstraintForControlInput(linear_constraint);
    }

    template <typename T, unsigned int DIM>
    void CBFQPGenerator<T, DIM>::addMinVelConstraints(const Vector &state) {
        // Get all minimum velocity constraints
        std::vector<LinearConstraint> linear_constraints = cbf_operations_ptr_->minVelConstraints(state);
        // Add each constraint to the QP problem
        for (size_t i = 0; i < linear_constraints.size(); ++i) {
            LinearConstraint &linear_constraint = linear_constraints.at(i);
            addLinearConstraintForControlInput(linear_constraint);
        }
    }

    template <typename T, unsigned int DIM>
    void CBFQPGenerator<T, DIM>::addMaxVelConstraints(const Vector &state) {
        // Get all maximum velocity constraints
        std::vector<LinearConstraint> linear_constraints = cbf_operations_ptr_->maxVelConstraints(state);
        // Add each constraint to the QP problem
        for (size_t i = 0; i < linear_constraints.size(); ++i) {
            LinearConstraint &linear_constraint = linear_constraints.at(i);
            addLinearConstraintForControlInput(linear_constraint);
        }
    }

    template <typename T, unsigned int DIM>
    void CBFQPGenerator<T, DIM>::addSafetyCBFConstraintWithSlackVariables(const Vector &state,
                                                                          const Vector &target_state,
                                                                          std::size_t neighbor_idx) {
        // Get the safety constraint
        LinearConstraint linear_constraint = cbf_operations_ptr_->safetyConstraint(state, target_state);

        // Create a row vector for slack variable coefficients (all zeros except at neighbor_idx)
        Row slack_coefficients = Row::Zero(slack_variables_.size());
        slack_coefficients(neighbor_idx) = -1; // Negative coefficient allows constraint relaxation

        // Add the constraint with slack variable to the QP problem
        addLinearConstraintForControlInputWithSlackVariables(linear_constraint, slack_coefficients);
    }

    template <typename T, unsigned int DIM>
    void CBFQPGenerator<T, DIM>::addLeftBorderConstraintWithSlackVariables(const Vector &state,
                                                                          const Vector &target_state,
                                                                          std::size_t neighbor_idx) {
        // Get the left field-of-view border constraint
        LinearConstraint linear_constraint = cbf_operations_ptr_->leftBorderConstraint(state, target_state);

        // Create a row vector for slack variable coefficients (all zeros except at neighbor_idx)
        Row slack_coefficients = Row::Zero(slack_variables_.size());
        slack_coefficients(neighbor_idx) = -1; // Negative coefficient allows constraint relaxation

        // Add the constraint with slack variable to the QP problem
        addLinearConstraintForControlInputWithSlackVariables(linear_constraint, slack_coefficients);
    }

    template <typename T, unsigned int DIM>
    void CBFQPGenerator<T, DIM>::addRightBorderConstraintWithSlackVariables(const Vector &state,
                                                                           const Vector &target_state,
                                                                           std::size_t neighbor_idx) {
        // Get the right field-of-view border constraint
        LinearConstraint linear_constraint = cbf_operations_ptr_->rightBorderConstraint(state, target_state);

        // Create a row vector for slack variable coefficients (all zeros except at neighbor_idx)
        Row slack_coefficients = Row::Zero(slack_variables_.size());
        slack_coefficients(neighbor_idx) = -1; // Negative coefficient allows constraint relaxation

        // Add the constraint with slack variable to the QP problem
        addLinearConstraintForControlInputWithSlackVariables(linear_constraint, slack_coefficients);
    }

    template <typename T, unsigned int DIM>
    void CBFQPGenerator<T, DIM>::addRangeConstraintWithSlackVariables(const Vector &state,
                                                                      const Vector &target_state,
                                                                      std::size_t neighbor_idx) {
        // Get the sensing range constraint
        LinearConstraint linear_constraint = cbf_operations_ptr_->rangeConstraint(state, target_state);

        // Create a row vector for slack variable coefficients (all zeros except at neighbor_idx)
        Row slack_coefficients = Row::Zero(slack_variables_.size());
        slack_coefficients(neighbor_idx) = -1; // Negative coefficient allows constraint relaxation

        // Add the constraint with slack variable to the QP problem
        addLinearConstraintForControlInputWithSlackVariables(linear_constraint, slack_coefficients);
    }

    template <typename T, unsigned int DIM>
    void CBFQPGenerator<T, DIM>::addControlBoundConstraint(const VectorDIM &u_min, const VectorDIM &u_max) {
        // Get bounds for control inputs
        DecisionVariableBounds decision_variable_bounds = cbf_operations_ptr_->controlBoundConstraint(u_min, u_max);
        // Add these bounds to the QP problem
        addDecisionVariableBoundsForControlInput(decision_variable_bounds);
    }

    template <typename T, unsigned int DIM>
    void
    CBFQPGenerator<T, DIM>::addCostAdditionForControlInput(const CostAddition &cost_addition) {
        // Small value to check if a coefficient is approximately zero
        constexpr T epsilon = std::numeric_limits<T>::epsilon() * T(100.0);
        size_t num_decision_variables = variables_.size();

        // Validate that the cost addition structure matches our variable count
        if (num_decision_variables != cost_addition.linear_term().rows() ||
            num_decision_variables != cost_addition.quadratic_term().rows() ||
            num_decision_variables != cost_addition.quadratic_term().cols()) {
            throw std::runtime_error("CBFQPGenerator::addCostAdditionForControlInput:"
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
                const qpcpp::Variable<T> *var2_ptr = variables_.at(j);
                if (!math::isApproximatelyEqual<T>(cost_addition.quadratic_term()(i,j), T(0.0), epsilon)) {
                    problem_.cost_function()->addQuadraticTerm(var1_ptr, var2_ptr, cost_addition.quadratic_term()(i,j));
                }
            }
        }
    }

    template <typename T, unsigned int DIM>
    void CBFQPGenerator<T, DIM>::addCostAdditionForSlackVariables(const CostAddition &cost_addition) {
        // Small value to check if a coefficient is approximately zero
        constexpr T epsilon = std::numeric_limits<T>::epsilon() * T(100.0);
        size_t num_slack_variables = slack_variables_.size();

        // Validate that the cost addition structure matches our slack variable count
        if (num_slack_variables != cost_addition.linear_term().rows() ||
            num_slack_variables != cost_addition.quadratic_term().rows() ||
            num_slack_variables != cost_addition.quadratic_term().cols()) {
            throw std::runtime_error("CBFQPGenerator::addCostAdditionForSlackVariables:"
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
                const qpcpp::Variable<T> *var2_ptr = slack_variables_.at(j);
                if (!math::isApproximatelyEqual<T>(cost_addition.quadratic_term()(i,j), T(0.0), epsilon)) {
                    problem().cost_function()->addQuadraticTerm(var1_ptr, var2_ptr, cost_addition.quadratic_term()(i,j));
                }
            }
        }
    }

    template <typename T, unsigned int DIM>
    void
    CBFQPGenerator<T, DIM>::addLinearConstraintForControlInput(const LinearConstraint &linear_constraint) {
        size_t num_decision_variables = variables_.size();

        // Validate that the constraint structure matches our variable count
        if (num_decision_variables != linear_constraint.coefficients().cols()) {
            throw std::runtime_error("CBFQPGenerator::"
                                     "addLinearConstraintForControlInput:"
                                     " number of decision variables of the controlInput does not match the "
                                     "LinearConstraint structure");
        }

        // Initialize a new linear constraint with lower and upper bounds
        qpcpp::LinearConstraint<T>* qpcpp_linear_constraint = problem_.addLinearConstraint(
                linear_constraint.lower_bound(),
                linear_constraint.upper_bound());

        // Set coefficients for each decision variable in this constraint
        for (std::size_t decision_variable_idx = 0;
             decision_variable_idx < num_decision_variables;
             ++decision_variable_idx) {
            const qpcpp::Variable<T>* var_ptr = variables_.at(decision_variable_idx);
            qpcpp_linear_constraint->setCoefficient(var_ptr, linear_constraint.coefficients()(decision_variable_idx));
        }
    }

    template <typename T, unsigned int DIM>
    void CBFQPGenerator<T, DIM>::addLinearConstraintForControlInputWithSlackVariables(
            const LinearConstraint &linear_constraint,
            const Row &slack_coefficients) {
        size_t num_decision_variables = variables_.size();
        size_t num_slack_variables = slack_variables_.size();

        // Validate that the constraint and slack coefficient structures match our variable counts
        if (num_decision_variables != linear_constraint.coefficients().cols()) {
            throw std::runtime_error("CBFQPGenerator::"
                                     "addLinearConstraintForControlInputWithSlackVariables:"
                                     " number of decision variables of the controlInput does not match the "
                                     "LinearConstraint structure");
        }
        if (num_slack_variables != slack_coefficients.size()) {
            throw std::runtime_error("CBFQPGenerator::"
                                     "addLinearConstraintForControlInputWithSlackVariables:"
                                     " number of slack variables does not match the "
                                     "slack_coefficients structure");
        }

        // Initialize a new linear constraint with lower and upper bounds
        qpcpp::LinearConstraint<T>* qpcpp_linear_constraint = problem_.addLinearConstraint(
                linear_constraint.lower_bound(),
                linear_constraint.upper_bound());

        // Set coefficients for each decision variable in this constraint
        for (std::size_t decision_variable_idx = 0;
             decision_variable_idx < num_decision_variables;
             ++decision_variable_idx) {
            const qpcpp::Variable<T>* var_ptr = variables_.at(decision_variable_idx);
            qpcpp_linear_constraint->setCoefficient(var_ptr, linear_constraint.coefficients()(decision_variable_idx));
        }

        // Set coefficients for each slack variable in this constraint
        for (std::size_t slack_variable_idx = 0;
             slack_variable_idx < num_slack_variables;
             ++slack_variable_idx) {
            const qpcpp::Variable<T>* var_ptr = slack_variables_.at(slack_variable_idx);
            qpcpp_linear_constraint->setCoefficient(var_ptr, slack_coefficients(slack_variable_idx));
        }
    }

    template <typename T, unsigned int DIM>
    void
    CBFQPGenerator<T, DIM>::addDecisionVariableBoundsForControlInput(const DecisionVariableBounds& decision_variable_bounds) {
        size_t num_decision_variables = variables_.size();

        // Validate that the bounds structure matches our variable count
        if (num_decision_variables != decision_variable_bounds.lower_bounds().rows() ||
            num_decision_variables != decision_variable_bounds.upper_bounds().rows()) {
            throw std::invalid_argument("CBFQPGenerator::"
                                        "addDecisionVariableBoundsForControlInput:"
                                        " number of decision variables of the controlInput does not match the "
                                        "DecisionVariablesBounds structure");
        }

        // Set the bounds for each decision variable
        for (int decision_variable_idx = 0;
             decision_variable_idx < DIM;
             ++decision_variable_idx) {
            qpcpp::Variable<T>* var_ptr =
                    variables_.at(decision_variable_idx);

            // Take the more restrictive of the existing and new lower bounds
            var_ptr->set_min(std::max(
                    var_ptr->min(),
                    decision_variable_bounds.lower_bounds()(decision_variable_idx)));

            // Take the more restrictive of the existing and new upper bounds
            var_ptr->set_max(std::min(
                    var_ptr->max(),
                    decision_variable_bounds.upper_bounds()(decision_variable_idx)));
        }
    }

    template <typename T, unsigned int DIM>
    typename CBFQPGenerator<T, DIM>::VectorDIM
    CBFQPGenerator<T, DIM>::generatorCBFControlInput() const {
        assert(DIM == cbf_operations_ptr_->numDecisionVariables());

        // Extract the solution values from the QP solver into a vector
        VectorDIM solution;
        for (size_t d = 0; d < DIM; ++d) {
            solution(d) = variables_.at(d)->solution_value();
        }
        return solution;
    }

    template class CBFQPGenerator<double, 3U>;
//    template class CBFQPGenerator<float, 3U>;

} // cbf