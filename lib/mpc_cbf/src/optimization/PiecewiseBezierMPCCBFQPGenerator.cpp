//
// Created by lishuo on 9/21/24.
//

#include <mpc_cbf/optimization/PiecewiseBezierMPCCBFQPGenerator.h>

namespace mpc_cbf {
    template <typename T, unsigned int DIM>
    void PiecewiseBezierMPCCBFQPGenerator<T, DIM>::addPiecewise(
            std::unique_ptr<PiecewiseBezierMPCCBFQPOperations> &&piecewise_mpc_cbf_operations_ptr, int num_neighbors, bool slack_mode) {
        // init the PiecewiseBezierMPCQPGenerator API
        std::unique_ptr<PiecewiseBezierMPCQPOperations> piecewise_mpc_operations_ptr = piecewise_mpc_cbf_operations_ptr->piecewise_mpc_operations_ptr();
        piecewise_mpc_qp_generator_ptr_->addPiecewise(std::move(piecewise_mpc_operations_ptr));
        piecewise_mpc_cbf_operations_ptr_ = std::move(piecewise_mpc_cbf_operations_ptr);
        // add slack variable to the problem
        if (slack_mode) {
            for (size_t i = 0; i < num_neighbors; ++i) {
                qpcpp::Variable<T>* variable_ptr = problem().addVariable(0, std::numeric_limits<T>::max()); // slack variables are non-negative
                slack_variables_.push_back(variable_ptr);
            }
        }
    }

    template <typename T, unsigned int DIM>
    qpcpp::Problem<T> &PiecewiseBezierMPCCBFQPGenerator<T, DIM>::problem() {return piecewise_mpc_qp_generator_ptr_->problem();}

    template <typename T, unsigned int DIM>
    std::shared_ptr<typename PiecewiseBezierMPCCBFQPGenerator<T, DIM>::PiecewiseBezierMPCQPGenerator>
    PiecewiseBezierMPCCBFQPGenerator<T, DIM>::piecewise_mpc_qp_generator_ptr() {return piecewise_mpc_qp_generator_ptr_;}

    template <typename T, unsigned int DIM>
    void PiecewiseBezierMPCCBFQPGenerator<T, DIM>::addSlackCost(const std::vector<double> &slack_weights) {
        CostAddition cost_addition = piecewise_mpc_cbf_operations_ptr_->slackCost(slack_weights);
        addCostAdditionForSlackVariables(cost_addition);
    }

    template <typename T, unsigned int DIM>
    void PiecewiseBezierMPCCBFQPGenerator<T, DIM>::addSafetyCBFConstraint(const State &current_state,
                                                                          const Vector &other_pos,
                                                                          T slack_value) {
        LinearConstraint linear_constraint = piecewise_mpc_cbf_operations_ptr_->safetyCBFConstraint(current_state, other_pos, slack_value);
        piecewise_mpc_qp_generator_ptr_->addLinearConstraintForPiecewise(linear_constraint);
    }

    template <typename T, unsigned int DIM>
    void PiecewiseBezierMPCCBFQPGenerator<T, DIM>::addFovLBConstraint(const State &current_state,
                                                                      const Vector &other_pos,
                                                                      T slack_value) {
        std::vector<LinearConstraint> linear_constraints = piecewise_mpc_cbf_operations_ptr_->fovLBConstraint(current_state, other_pos, slack_value);
        for (size_t i = 0; i < linear_constraints.size(); ++i) {
            piecewise_mpc_qp_generator_ptr_->addLinearConstraintForPiecewise(linear_constraints.at(i));
        }
    }

    template <typename T, unsigned int DIM>
    void PiecewiseBezierMPCCBFQPGenerator<T, DIM>::addFovRBConstraint(const State &current_state,
                                                                      const Vector &other_pos,
                                                                      T slack_value) {
        std::vector<LinearConstraint> linear_constraints = piecewise_mpc_cbf_operations_ptr_->fovRBConstraint(current_state, other_pos, slack_value);
        for (size_t i = 0; i < linear_constraints.size(); ++i) {
            piecewise_mpc_qp_generator_ptr_->addLinearConstraintForPiecewise(linear_constraints.at(i));
        }
    }

    template <typename T, unsigned int DIM>
    void PiecewiseBezierMPCCBFQPGenerator<T, DIM>::addPredSafetyCBFConstraints(const std::vector<State> &pred_states,
                                                                               const Vector &other_pos,
                                                                               const std::vector<T>& slack_values) {
        std::vector<LinearConstraint> linear_constraints = piecewise_mpc_cbf_operations_ptr_->predSafetyCBFConstraints(pred_states, other_pos);
        for (size_t i = 0; i < linear_constraints.size(); ++i) {
            piecewise_mpc_qp_generator_ptr_->addLinearConstraintForPiecewise(linear_constraints.at(i));
        }
    }

    template <typename T, unsigned int DIM>
    void PiecewiseBezierMPCCBFQPGenerator<T, DIM>::addPredFovLBConstraints(const std::vector<State> &pred_states,
                                                                           const Vector &other_pos,
                                                                           const std::vector<T>& slack_values) {
        std::vector<LinearConstraint> linear_constraints = piecewise_mpc_cbf_operations_ptr_->predFovLBConstraints(pred_states, other_pos);
        for (size_t i = 0; i < linear_constraints.size(); ++i) {
            piecewise_mpc_qp_generator_ptr_->addLinearConstraintForPiecewise(linear_constraints.at(i));
        }
    }

    template <typename T, unsigned int DIM>
    void PiecewiseBezierMPCCBFQPGenerator<T, DIM>::addPredFovRBConstraints(const std::vector<State> &pred_states,
                                                                           const Vector &other_pos,
                                                                           const std::vector<T>& slack_values) {
        std::vector<LinearConstraint> linear_constraints = piecewise_mpc_cbf_operations_ptr_->predFovRBConstraints(pred_states, other_pos);
        for (size_t i = 0; i < linear_constraints.size(); ++i) {
            piecewise_mpc_qp_generator_ptr_->addLinearConstraintForPiecewise(linear_constraints.at(i));
        }
    }

    template <typename T, unsigned int DIM>
    void PiecewiseBezierMPCCBFQPGenerator<T, DIM>::addSafetyCBFConstraintWithSlackVariables(
            const State &current_state,
            const Vector &other_pos) {
        LinearConstraint linear_constraint = piecewise_mpc_cbf_operations_ptr_->safetyCBFConstraint(current_state, other_pos);
        Row slack_coefficients = -Row::Ones(slack_variables_.size());
        addLinearConstraintForPiecewiseWithSlackVariables(linear_constraint, slack_coefficients);
    }

    template <typename T, unsigned int DIM>
    void PiecewiseBezierMPCCBFQPGenerator<T, DIM>::addFovLBConstraintWithSlackVariables(
            const State &current_state,
            const Vector &other_pos) {
        std::vector<LinearConstraint> linear_constraints = piecewise_mpc_cbf_operations_ptr_->fovLBConstraint(current_state, other_pos);
        Row slack_coefficients = -Row::Ones(slack_variables_.size());
        for (size_t i = 0; i < linear_constraints.size(); ++i) {
            addLinearConstraintForPiecewiseWithSlackVariables(linear_constraints.at(i), slack_coefficients);
        }
    }

    template <typename T, unsigned int DIM>
    void PiecewiseBezierMPCCBFQPGenerator<T, DIM>::addFovRBConstraintWithSlackVariables(
            const State &current_state,
            const Vector &other_pos) {
        std::vector<LinearConstraint> linear_constraints = piecewise_mpc_cbf_operations_ptr_->fovRBConstraint(current_state, other_pos);
        Row slack_coefficients = -Row::Ones(slack_variables_.size());
        for (size_t i = 0; i < linear_constraints.size(); ++i) {
            addLinearConstraintForPiecewiseWithSlackVariables(linear_constraints.at(i), slack_coefficients);
        }
    }

    template <typename T, unsigned int DIM>
    void PiecewiseBezierMPCCBFQPGenerator<T, DIM>::addPredFovLBConstraintsWithSlackVariables(
            const std::vector<State> &pred_states, const Vector &other_pos) {
        std::vector<LinearConstraint> linear_constraints = piecewise_mpc_cbf_operations_ptr_->predFovLBConstraints(pred_states, other_pos);
        Row slack_coefficients = -Row::Ones(slack_variables_.size());
        for (size_t i = 0; i < linear_constraints.size(); ++i) {
            addLinearConstraintForPiecewiseWithSlackVariables(linear_constraints.at(i), slack_coefficients);
        }
    }

    template <typename T, unsigned int DIM>
    void PiecewiseBezierMPCCBFQPGenerator<T, DIM>::addPredFovRBConstraintsWithSlackVariables(
            const std::vector<State> &pred_states, const Vector &other_pos) {
        std::vector<LinearConstraint> linear_constraints = piecewise_mpc_cbf_operations_ptr_->predFovRBConstraints(pred_states, other_pos);
        Row slack_coefficients = -Row::Ones(slack_variables_.size());
        for (size_t i = 0; i < linear_constraints.size(); ++i) {
            addLinearConstraintForPiecewiseWithSlackVariables(linear_constraints.at(i), slack_coefficients);
        }
    }

    template <typename T, unsigned int DIM>
    void PiecewiseBezierMPCCBFQPGenerator<T, DIM>::addCostAdditionForSlackVariables(
            const CostAddition &cost_addition) {
        constexpr T epsilon = std::numeric_limits<T>::epsilon() * T(100.0);
        size_t num_slack_variables = slack_variables_.size();

        // number decision variable number error
        if (num_slack_variables != cost_addition.linear_term().rows() ||
            num_slack_variables != cost_addition.quadratic_term().rows() ||
            num_slack_variables != cost_addition.quadratic_term().cols()) {
            throw std::runtime_error("PiecewiseBezierMPCCBFQPGenerator::addCostAdditionForSlackVariables:"
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
    void PiecewiseBezierMPCCBFQPGenerator<T, DIM>::addLinearConstraintForPiecewiseWithSlackVariables(
            const LinearConstraint &linear_constraint,
            const Row &slack_coefficients) {

        size_t num_decision_variables = piecewise_mpc_qp_generator_ptr_->variables().size();
        size_t num_slack_variables = slack_variables_.size();

        // linear constraint structure error.
        if (num_decision_variables != linear_constraint.coefficients().cols()) {
            throw std::runtime_error("PiecewiseBezierMPCCBFQPGenerator::"
                                     "addLinearConstraintForPiecewiseWithSlackVariables:"
                                     " number of decision variables of the piecewise does not match the "
                                     "LinearConstraint structure");
        }
        if (num_slack_variables != slack_coefficients.size()) {
            throw std::runtime_error("PiecewiseBezierMPCCBFQPGenerator::"
                                     "addLinearConstraintForPiecewiseWithSlackVariables:"
                                     " number of slack variables does not match the "
                                     "slack_coefficients structure");
        }

        // initialize linear constraint with lower bound and upper bound
        qpcpp::LinearConstraint<T>* qpcpp_linear_constraint = problem().addLinearConstraint(
                linear_constraint.lower_bound(),
                linear_constraint.upper_bound());
        // set coefficients for this linear constraint (curve variables part)
        for (std::size_t decision_variable_idx = 0;
             decision_variable_idx < num_decision_variables;
             ++decision_variable_idx) {
            const qpcpp::Variable<T>* var_ptr = piecewise_mpc_qp_generator_ptr_->variables().at(decision_variable_idx);
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

    template class PiecewiseBezierMPCCBFQPGenerator<double, 3U>;

} // mpc_cbf