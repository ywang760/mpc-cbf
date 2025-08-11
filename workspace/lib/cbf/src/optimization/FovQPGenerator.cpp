//
// Created on 4/10/25.
//

#include <cbf/optimization/FovQPGenerator.h>

namespace cbf {
template <typename T, unsigned int DIM>
FovQPGenerator<T, DIM>::FovQPGenerator(std::shared_ptr<FovCBF> cbf, int num_robots, bool slack_mode)
    : CBFQPGeneratorBase<T, DIM>(), cbf_(cbf) {

    // Add slack variables to the problem if in slack mode
    // Slack variables allow constraints to be violated by paying a penalty in the cost function
    slack_mode_ = slack_mode;
    if (slack_mode) {
        for (size_t i = 0; i < num_robots; ++i) {
            qpcpp::Variable<T>* variable_ptr = this->problem().addVariable(
                0, std::numeric_limits<T>::max()); // slack variables are non-negative
            this->slack_variables_.push_back(variable_ptr);
        }
    }
}

template <typename T, unsigned int DIM>
void FovQPGenerator<T, DIM>::addSafetyConstraint(const Vector& state, const Vector& target_state,
                                                 bool use_slack, std::size_t slack_idx) {
    // Get safety constraint coefficients from the CBF implementation
    Vector coefficients = -1.0 * cbf_->getSafetyConstraints(state, target_state);

    // Get the safety constraint bound from the CBF implementation
    T bound = cbf_->getSafetyBound(state, target_state);

    // Create a linear constraint with lower bound negative infinity
    LinearConstraint linear_constraint(coefficients, std::numeric_limits<T>::lowest(), bound);

    // Handle constraint with or without slack variable
    if (use_slack && !this->slack_variables_.empty()) {
        // Create a row vector for slack variable coefficients (all zeros except at slack_idx)
        Row slack_coefficients = Row::Zero(this->slack_variables_.size());
        slack_coefficients(slack_idx) = -1; // Negative coefficient allows constraint relaxation

        // Add the constraint with slack variable to the QP problem
        this->addLinearConstraintForControlInputWithSlackVariables(
            linear_constraint, slack_coefficients, this->slack_variables_);
    } else {
        // Add the constraint without slack variable to the QP problem
        this->addLinearConstraintForControlInput(linear_constraint);
    }
}

template <typename T, unsigned int DIM>
void FovQPGenerator<T, DIM>::addLeftBorderConstraint(const Vector& state,
                                                     const Vector& target_state, bool use_slack,
                                                     std::size_t slack_idx) {
    // Get left border constraint coefficients from the CBF implementation
    Vector coefficients = -1.0 * cbf_->getLBConstraints(state, target_state);

    // Get the left border constraint bound from the CBF implementation
    T bound = cbf_->getLBBound(state, target_state);

    // Create a linear constraint with lower bound negative infinity
    LinearConstraint linear_constraint(coefficients, std::numeric_limits<T>::lowest(), bound);

    // Handle constraint with or without slack variable
    if (use_slack && !this->slack_variables_.empty()) {
        // Create a row vector for slack variable coefficients (all zeros except at slack_idx)
        Row slack_coefficients = Row::Zero(this->slack_variables_.size());
        slack_coefficients(slack_idx) = -1; // Negative coefficient allows constraint relaxation

        // Add the constraint with slack variable to the QP problem
        this->addLinearConstraintForControlInputWithSlackVariables(
            linear_constraint, slack_coefficients, this->slack_variables_);
    } else {
        // Add the constraint without slack variable to the QP problem
        this->addLinearConstraintForControlInput(linear_constraint);
    }
}

template <typename T, unsigned int DIM>
void FovQPGenerator<T, DIM>::addRightBorderConstraint(const Vector& state,
                                                      const Vector& target_state, bool use_slack,
                                                      std::size_t slack_idx) {
    // Get right border constraint coefficients from the CBF implementation
    Vector coefficients = -1.0 * cbf_->getRBConstraints(state, target_state);

    // Get the right border constraint bound from the CBF implementation
    T bound = cbf_->getRBBound(state, target_state);

    // Create a linear constraint with lower bound negative infinity
    LinearConstraint linear_constraint(coefficients, std::numeric_limits<T>::lowest(), bound);

    // Handle constraint with or without slack variable
    if (use_slack && !this->slack_variables_.empty()) {
        // Create a row vector for slack variable coefficients (all zeros except at slack_idx)
        Row slack_coefficients = Row::Zero(this->slack_variables_.size());
        slack_coefficients(slack_idx) = -1; // Negative coefficient allows constraint relaxation

        // Add the constraint with slack variable to the QP problem
        this->addLinearConstraintForControlInputWithSlackVariables(
            linear_constraint, slack_coefficients, this->slack_variables_);
    } else {
        // Add the constraint without slack variable to the QP problem
        this->addLinearConstraintForControlInput(linear_constraint);
    }
}

template <typename T, unsigned int DIM>
void FovQPGenerator<T, DIM>::addRangeConstraint(const Vector& state, const Vector& target_state,
                                                bool use_slack, std::size_t slack_idx) {
    // Get range constraint coefficients from the CBF implementation
    Vector coefficients = -1.0 * cbf_->getRangeConstraints(state, target_state);

    // Get the range constraint bound from the CBF implementation
    T bound = cbf_->getRangeBound(state, target_state);

    // Create a linear constraint with lower bound negative infinity
    LinearConstraint linear_constraint(coefficients, std::numeric_limits<T>::lowest(), bound);

    // Handle constraint with or without slack variable
    if (use_slack && !this->slack_variables_.empty()) {
        // Create a row vector for slack variable coefficients (all zeros except at slack_idx)
        Row slack_coefficients = Row::Zero(this->slack_variables_.size());
        slack_coefficients(slack_idx) = -1; // Negative coefficient allows constraint relaxation

        // Add the constraint with slack variable to the QP problem
        this->addLinearConstraintForControlInputWithSlackVariables(
            linear_constraint, slack_coefficients, this->slack_variables_);
    } else {
        // Add the constraint without slack variable to the QP problem
        this->addLinearConstraintForControlInput(linear_constraint);
    }
}

template <typename T, unsigned int DIM>
void FovQPGenerator<T, DIM>::addMinVelConstraints(const Vector& state) {
    // Get minimum velocity constraint coefficients for all dimensions
    typename CBFQPGeneratorBase<T, DIM>::Matrix coefficient_matrix =
        cbf_->getMinVelContraints(state);

    // Get minimum velocity constraint bounds for all dimensions
    Vector bounds = cbf_->getMinVelBounds(state);

    // Ensure the coefficient matrix and bounds vector dimensions match
    assert(coefficient_matrix.rows() == bounds.size());

    // Create and add a linear constraint for each minimum velocity constraint
    for (size_t i = 0; i < coefficient_matrix.rows(); ++i) {
        // Extract the coefficients for this constraint
        Vector coefficients = -1.0 * coefficient_matrix.row(i);

        // Extract the bound for this constraint
        T bound = bounds(i);

        // Create a linear constraint with lower bound negative infinity
        LinearConstraint linear_constraint(coefficients, std::numeric_limits<T>::lowest(), bound);
        this->addLinearConstraintForControlInput(linear_constraint);
    }
}

template <typename T, unsigned int DIM>
void FovQPGenerator<T, DIM>::addMaxVelConstraints(const Vector& state) {
    // Get maximum velocity constraint coefficients for all dimensions
    typename CBFQPGeneratorBase<T, DIM>::Matrix coefficient_matrix =
        cbf_->getMaxVelContraints(state);

    // Get maximum velocity constraint bounds for all dimensions
    Vector bounds = cbf_->getMaxVelBounds(state);

    // Ensure the coefficient matrix and bounds vector dimensions match
    assert(coefficient_matrix.rows() == bounds.size());

    // Create and add a linear constraint for each maximum velocity constraint
    for (size_t i = 0; i < coefficient_matrix.rows(); ++i) {
        // Extract the coefficients for this constraint
        Vector coefficients = -1.0 * coefficient_matrix.row(i);

        // Extract the bound for this constraint
        T bound = bounds(i);

        // Create a linear constraint with lower bound negative infinity
        LinearConstraint linear_constraint(coefficients, std::numeric_limits<T>::lowest(), bound);
        this->addLinearConstraintForControlInput(linear_constraint);
    }
}

// Explicit template instantiation
template class FovQPGenerator<double, 3U>;

} // namespace cbf