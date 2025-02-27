//
// Created by lishuo on 8/31/24.
//

#include <cbf/optimization/CBFQPGenerator.h>

namespace cbf {
    template <typename T, unsigned int DIM>
    void
    CBFQPGenerator<T, DIM>::addCBFOperations(std::unique_ptr<CBFQPOperations> cbf_operations_ptr, int num_neighbors, bool slack_mode) {
        cbf_operations_ptr_ = std::move(cbf_operations_ptr);
        // allocate the variable
        size_t num_decision_variables_of_control_input = cbf_operations_ptr_->numDecisionVariables();
        for (size_t decision_variable_idx = 0;
             decision_variable_idx < num_decision_variables_of_control_input;
             ++decision_variable_idx) {
            qpcpp::Variable<T>* variable_ptr = problem_.addVariable();
            variables_.push_back(variable_ptr);
        }
        // add slack variable to the problem
        if (slack_mode) {
            for (size_t i = 0; i < num_neighbors; ++i) {
                qpcpp::Variable<T>* variable_ptr = problem().addVariable(0, std::numeric_limits<T>::max()); // slack variables are non-negative
                slack_variables_.push_back(variable_ptr);
            }
        }
    }

    template <typename T, unsigned int DIM>
    qpcpp::Problem<T>& CBFQPGenerator<T, DIM>::problem() {
        return problem_;
    }

    template <typename T, unsigned int DIM>
    void CBFQPGenerator<T, DIM>::addDesiredControlCost(const VectorDIM& desired_u) {
        CostAddition cost_addition = cbf_operations_ptr_->desiredControlCost(desired_u);
        addCostAdditionForControlInput(cost_addition);
    }

    template <typename T, unsigned int DIM>
    void CBFQPGenerator<T, DIM>::addSlackCost(const std::vector<T> &slack_weights) {
        CostAddition cost_addition = cbf_operations_ptr_->slackCost(slack_weights);
        addCostAdditionForSlackVariables(cost_addition);
    }

    template <typename T, unsigned int DIM>
    void CBFQPGenerator<T, DIM>::addSafetyConstraint(const Vector &state,
                                                     const Vector &target_state) {
        LinearConstraint linear_constraint = cbf_operations_ptr_->safetyConstraint(state, target_state);
        addLinearConstraintForControlInput(linear_constraint);
    }

    template <typename T, unsigned int DIM>
    void CBFQPGenerator<T, DIM>::addLeftBorderConstraint(const Vector &state,
                                                         const Vector &target_state) {
        LinearConstraint linear_constraint = cbf_operations_ptr_->leftBorderConstraint(state, target_state);
        addLinearConstraintForControlInput(linear_constraint);
    }

    template <typename T, unsigned int DIM>
    void CBFQPGenerator<T, DIM>::addRightBorderConstraint(const Vector &state,
                                                          const Vector &target_state) {
        LinearConstraint linear_constraint = cbf_operations_ptr_->rightBorderConstraint(state, target_state);
        addLinearConstraintForControlInput(linear_constraint);
    }

    template <typename T, unsigned int DIM>
    void CBFQPGenerator<T, DIM>::addRangeConstraint(const Vector &state,
                                                    const Vector &target_state) {
        LinearConstraint linear_constraint = cbf_operations_ptr_->rangeConstraint(state, target_state);
        addLinearConstraintForControlInput(linear_constraint);
    }

    template <typename T, unsigned int DIM>
    void CBFQPGenerator<T, DIM>::addMinVelConstraints(const Vector &state) {
        std::vector<LinearConstraint> linear_constraints = cbf_operations_ptr_->minVelConstraints(state);
        for (size_t i = 0; i < linear_constraints.size(); ++i) {
            LinearConstraint &linear_constraint = linear_constraints.at(i);
            addLinearConstraintForControlInput(linear_constraint);
        }
    }

    template <typename T, unsigned int DIM>
    void CBFQPGenerator<T, DIM>::addMaxVelConstraints(const Vector &state) {
        std::vector<LinearConstraint> linear_constraints = cbf_operations_ptr_->maxVelConstraints(state);
        for (size_t i = 0; i < linear_constraints.size(); ++i) {
            LinearConstraint &linear_constraint = linear_constraints.at(i);
            addLinearConstraintForControlInput(linear_constraint);
        }
    }

    template <typename T, unsigned int DIM>
    void CBFQPGenerator<T, DIM>::addSafetyCBFConstraintWithSlackVariables(const Vector &state,
                                                                          const Vector &target_state,
                                                                          std::size_t neighbor_idx) {
        LinearConstraint linear_constraint = cbf_operations_ptr_->safetyConstraint(state, target_state);
        Row slack_coefficients = Row::Zero(slack_variables_.size());
        slack_coefficients(neighbor_idx) = -1;
        addLinearConstraintForControlInputWithSlackVariables(linear_constraint, slack_coefficients);
    }

    template <typename T, unsigned int DIM>
    void CBFQPGenerator<T, DIM>::addLeftBorderConstraintWithSlackVariables(const Vector &state,
                                                                          const Vector &target_state,
                                                                          std::size_t neighbor_idx) {
        LinearConstraint linear_constraint = cbf_operations_ptr_->leftBorderConstraint(state, target_state);
        Row slack_coefficients = Row::Zero(slack_variables_.size());
        slack_coefficients(neighbor_idx) = -1;
        addLinearConstraintForControlInputWithSlackVariables(linear_constraint, slack_coefficients);
    }

    template <typename T, unsigned int DIM>
    void CBFQPGenerator<T, DIM>::addRightBorderConstraintWithSlackVariables(const Vector &state,
                                                                           const Vector &target_state,
                                                                           std::size_t neighbor_idx) {
        LinearConstraint linear_constraint = cbf_operations_ptr_->rightBorderConstraint(state, target_state);
        Row slack_coefficients = Row::Zero(slack_variables_.size());
        slack_coefficients(neighbor_idx) = -1;
        addLinearConstraintForControlInputWithSlackVariables(linear_constraint, slack_coefficients);
    }

    template <typename T, unsigned int DIM>
    void CBFQPGenerator<T, DIM>::addRangeConstraintWithSlackVariables(const Vector &state,
                                                                      const Vector &target_state,
                                                                      std::size_t neighbor_idx) {
        LinearConstraint linear_constraint = cbf_operations_ptr_->rangeConstraint(state, target_state);
        Row slack_coefficients = Row::Zero(slack_variables_.size());
        slack_coefficients(neighbor_idx) = -1;
        addLinearConstraintForControlInputWithSlackVariables(linear_constraint, slack_coefficients);
    }

    template <typename T, unsigned int DIM>
    void CBFQPGenerator<T, DIM>::addControlBoundConstraint(const VectorDIM &u_min, const VectorDIM &u_max) {
        DecisionVariableBounds decision_variable_bounds = cbf_operations_ptr_->controlBoundConstraint(u_min, u_max);
        addDecisionVariableBoundsForControlInput(decision_variable_bounds);
    }

    template <typename T, unsigned int DIM>
    void
    CBFQPGenerator<T, DIM>::addCostAdditionForControlInput(const CostAddition &cost_addition) {
        constexpr T epsilon = std::numeric_limits<T>::epsilon() * T(100.0);
        size_t num_decision_variables = variables_.size();
        // number decision variable number error
        if (num_decision_variables != cost_addition.linear_term().rows() ||
            num_decision_variables != cost_addition.quadratic_term().rows() ||
            num_decision_variables != cost_addition.quadratic_term().cols()) {
            throw std::runtime_error("CBFQPGenerator::addCostAdditionForControlInput:"
                                     " number of decision variables of the controlInput does not match the "
                                     "CostAddition structure");
        }

        // add the constant term
        if (cost_addition.constant() != 0) {
            problem_.cost_function()->add_constant(cost_addition.constant());
        }

        for (size_t i = 0; i < num_decision_variables; ++i) {
            const qpcpp::Variable<T>* var1_ptr = variables_.at(i);
            // add the linear term
            if (!math::isApproximatelyEqual<T>(cost_addition.linear_term()(i), T(0.0), epsilon)) {
                problem_.cost_function()->addLinearTerm(var1_ptr, cost_addition.linear_term()(i));
            }
            // add quadratic term
            for (std::size_t j = 0; j < num_decision_variables; ++j) {
                const qpcpp::Variable<T>* var2_ptr = variables_.at(j);
                // add quadratic term
                if (!math::isApproximatelyEqual<T>(cost_addition.quadratic_term()(i,j), T(0.0), epsilon)) {
                    problem_.cost_function()->addQuadraticTerm(var1_ptr, var2_ptr, cost_addition.quadratic_term()(i,j));
                }
            }
        }
    }

    template <typename T, unsigned int DIM>
    void CBFQPGenerator<T, DIM>::addCostAdditionForSlackVariables(const CostAddition &cost_addition) {
        constexpr T epsilon = std::numeric_limits<T>::epsilon() * T(100.0);
        size_t num_slack_variables = slack_variables_.size();

        // number decision variable number error
        if (num_slack_variables != cost_addition.linear_term().rows() ||
            num_slack_variables != cost_addition.quadratic_term().rows() ||
            num_slack_variables != cost_addition.quadratic_term().cols()) {
            throw std::runtime_error("CBFQPGenerator::addCostAdditionForSlackVariables:"
                                     " number of slack variables does not match the "
                                     "CostAddition structure");
        }

        // add the constant term
        if (cost_addition.constant() != 0) {
            problem().cost_function()->add_constant(cost_addition.constant());
        }

        for (size_t i = 0; i < num_slack_variables; ++i) {
            const qpcpp::Variable<T>* var1_ptr = slack_variables_.at(i);
            // add the linear term
            if (!math::isApproximatelyEqual<T>(cost_addition.linear_term()(i), T(0.0), epsilon)) {
                problem().cost_function()->addLinearTerm(var1_ptr, cost_addition.linear_term()(i));
            }
            // add quadratic term
            for (std::size_t j = 0; j < num_slack_variables; ++j) {
                const qpcpp::Variable<T>* var2_ptr = slack_variables_.at(j);
                // add quadratic term
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
        // linear constraint structure error.
        if (num_decision_variables != linear_constraint.coefficients().cols()) {
            throw std::runtime_error("CBFQPGenerator::"
                                     "addLinearConstraintForControlInput:"
                                     " number of decision variables of the controlInput does not match the "
                                     "LinearConstraint structure");
        }

        // initialize linear constraint with lower bound and upper bound
        qpcpp::LinearConstraint<T>* qpcpp_linear_constraint = problem_.addLinearConstraint(
                linear_constraint.lower_bound(),
                linear_constraint.upper_bound());
        // set coefficients for this linear constraint
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

        // linear constraint structure error.
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

        // initialize linear constraint with lower bound and upper bound
        qpcpp::LinearConstraint<T>* qpcpp_linear_constraint = problem_.addLinearConstraint(
                linear_constraint.lower_bound(),
                linear_constraint.upper_bound());
        // set coefficients for this linear constraint
        for (std::size_t decision_variable_idx = 0;
             decision_variable_idx < num_decision_variables;
             ++decision_variable_idx) {
            const qpcpp::Variable<T>* var_ptr = variables_.at(decision_variable_idx);
            qpcpp_linear_constraint->setCoefficient(var_ptr, linear_constraint.coefficients()(decision_variable_idx));
        }
        // set coefficients for this linear constraint (slack variables part)
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
        if (num_decision_variables != decision_variable_bounds.lower_bounds().rows() ||
            num_decision_variables != decision_variable_bounds.upper_bounds().rows()) {
            throw std::invalid_argument("CBFQPGenerator::"
                                        "addDecisionVariableBoundsForControlInput:"
                                        " number of decision variables of the controlInput does not match the "
                                        "DecisionVariablesBounds structure");
        }
        for (int decision_variable_idx = 0;
             decision_variable_idx < DIM;
             ++decision_variable_idx) {
            qpcpp::Variable<T>* var_ptr =
                    variables_.at(decision_variable_idx);

            var_ptr->set_min(std::max(
                    var_ptr->min(),
                    decision_variable_bounds.lower_bounds()(decision_variable_idx)));

            var_ptr->set_max(std::min(
                    var_ptr->max(),
                    decision_variable_bounds.upper_bounds()(decision_variable_idx)));
        }
    }

    template <typename T, unsigned int DIM>
    typename CBFQPGenerator<T, DIM>::VectorDIM
    CBFQPGenerator<T, DIM>::generatorCBFControlInput() const {
        assert(DIM == cbf_operations_ptr_->numDecisionVariables());
        VectorDIM solution;
        for (size_t d = 0; d < DIM; ++d) {
            solution(d) = variables_.at(d)->solution_value();
        }
        return solution;
    }

    template class CBFQPGenerator<double, 3U>;
//    template class CBFQPGenerator<float, 3U>;

} // cbf