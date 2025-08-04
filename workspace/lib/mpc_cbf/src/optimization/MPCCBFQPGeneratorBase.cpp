//
// Created on 8/4/25.
//

#include <mpc_cbf/optimization/MPCCBFQPGeneratorBase.h>

#include <limits>

namespace mpc_cbf {
    template <typename T, unsigned int DIM>
    MPCCBFQPGeneratorBase<T, DIM>::MPCCBFQPGeneratorBase() 
        : piecewise_mpc_qp_generator_ptr_(std::make_shared<PiecewiseBezierMPCQPGenerator>()),
          slack_mode_(false),
          num_neighbors_(0) {
    }

    template <typename T, unsigned int DIM>
    typename MPCCBFQPGeneratorBase<T, DIM>::Problem& 
    MPCCBFQPGeneratorBase<T, DIM>::problem() {
        return piecewise_mpc_qp_generator_ptr_->problem();
    }

    template <typename T, unsigned int DIM>
    std::shared_ptr<typename MPCCBFQPGeneratorBase<T, DIM>::PiecewiseBezierMPCQPGenerator>
    MPCCBFQPGeneratorBase<T, DIM>::piecewise_mpc_qp_generator_ptr() {
        return piecewise_mpc_qp_generator_ptr_;
    }

    template <typename T, unsigned int DIM>
    void MPCCBFQPGeneratorBase<T, DIM>::addSlackVariables(int num_neighbors) {
        num_neighbors_ = num_neighbors;
        slack_mode_ = true;
        
        // Add slack variables to the problem
        for (size_t i = 0; i < num_neighbors; ++i) {
            qpcpp::Variable<T>* variable_ptr = problem().addVariable(0, std::numeric_limits<T>::max()); // slack variables are non-negative
            slack_variables_.push_back(variable_ptr);
        }
    }

    template <typename T, unsigned int DIM>
    void MPCCBFQPGeneratorBase<T, DIM>::addCostAdditionForSlackVariables(const CostAddition& cost_addition) {
        constexpr T epsilon = std::numeric_limits<T>::epsilon() * T(100.0);
        size_t num_slack_variables = slack_variables_.size();

        // number decision variable number error
        if (num_slack_variables != cost_addition.linear_term().rows() ||
            num_slack_variables != cost_addition.quadratic_term().rows() ||
            num_slack_variables != cost_addition.quadratic_term().cols()) {
            throw std::runtime_error("MPCCBFQPGeneratorBase::addCostAdditionForSlackVariables:"
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
    void MPCCBFQPGeneratorBase<T, DIM>::addLinearConstraintForPiecewiseWithSlackVariables(
            const LinearConstraint& linear_constraint,
            const Row& slack_coefficients) {

        size_t num_decision_variables = piecewise_mpc_qp_generator_ptr_->variables().size();
        size_t num_slack_variables = slack_variables_.size();

        // linear constraint structure error.
        if (num_decision_variables != linear_constraint.coefficients().cols()) {
            throw std::runtime_error("MPCCBFQPGeneratorBase::"
                                     "addLinearConstraintForPiecewiseWithSlackVariables:"
                                     " number of decision variables of the piecewise does not match the "
                                     "LinearConstraint structure");
        }
        if (num_slack_variables != slack_coefficients.size()) {
            throw std::runtime_error("MPCCBFQPGeneratorBase::"
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

    // Explicit template instantiation
    template class MPCCBFQPGeneratorBase<double, 3U>;

} // mpc_cbf