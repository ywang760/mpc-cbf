//
// Created by lishuo on 8/31/24.
//

#include <cbf/optimization/CBFQPOperations.h>

namespace cbf {
    template <typename T, unsigned int DIM>
    CBFQPOperations<T, DIM>::CBFQPOperations(std::shared_ptr<FovCBF> cbf) : cbf_(cbf) {}

    template <typename T, unsigned int DIM>
    typename CBFQPOperations<T, DIM>::CostAddition
    CBFQPOperations<T, DIM>::desiredControlCost(const VectorDIM &desired_u) {
        Matrix quadratic_term(DIM, DIM);
        quadratic_term.setZero();
        Vector linear_term(DIM);
        linear_term.setZero();

        // Set the quadratic term to identity matrix for a standard quadratic norm
        quadratic_term.setIdentity();

        // Set the linear term as -2 * desired_u to create the expanded form of ||u - desired_u||^2
        // When expanded, ||u - desired_u||^2 = u^T * u - 2 * u^T * desired_u + desired_u^T * desired_u
        linear_term = -2.0 * desired_u;

        // Set the constant term as the squared norm of desired_u
        // This completes the expansion of ||u - desired_u||^2
        T constant_term = desired_u.transpose() * desired_u;

        return CostAddition(quadratic_term, linear_term, constant_term);
    }

    template <typename T, unsigned int DIM>
    typename CBFQPOperations<T, DIM>::CostAddition
    CBFQPOperations<T, DIM>::slackCost(const std::vector<T> &slack_weights) {
        // Initialize cost terms for slack variables
        Matrix quadratic_term(slack_weights.size(), slack_weights.size());
        quadratic_term.setZero();
        Vector linear_term(slack_weights.size());
        linear_term.setZero();

        // For slack variables, we only use linear costs to penalize their usage
        // No quadratic terms are used, making this a linear penalty
        for (std::size_t i = 0; i < slack_weights.size(); ++i) {
            linear_term(i) = slack_weights.at(i);
        }

        // Return the cost addition with zero constant term
        return CostAddition(quadratic_term, linear_term, 0);
    }

    template <typename T, unsigned int DIM>
    typename CBFQPOperations<T, DIM>::LinearConstraint
    CBFQPOperations<T, DIM>::safetyConstraint(const Vector &state,
                                              const Vector &target_state) {
        // Get safety constraint coefficients from the CBF implementation
        // The coefficients represent the gradient of the barrier function
        Vector coefficients = -1.0 * cbf_->getSafetyConstraints(state, target_state);

        // Get the safety constraint bound from the CBF implementation
        // This bound ensures the barrier function remains non-negative
        T bound = cbf_->getSafetyBound(state, target_state);

        // Create a linear constraint with lower bound negative infinity
        // This represents: coefficients^T * u <= bound
        return LinearConstraint(coefficients, std::numeric_limits<T>::lowest(), bound);
    }

    template <typename T, unsigned int DIM>
    typename CBFQPOperations<T, DIM>::LinearConstraint
    CBFQPOperations<T, DIM>::leftBorderConstraint(const Vector &state,
                                                  const Vector &target_state) {
        // Get left border constraint coefficients from the CBF implementation
        // These coefficients ensure the target stays within the left FOV boundary
        Vector coefficients = -1.0 * cbf_->getLBConstraints(state, target_state);

        // Get the left border constraint bound from the CBF implementation
        T bound = cbf_->getLBBound(state, target_state);

        // Create a linear constraint with lower bound negative infinity
        // This represents: coefficients^T * u <= bound
        return LinearConstraint(coefficients, std::numeric_limits<T>::lowest(), bound);
    }

    template <typename T, unsigned int DIM>
    typename CBFQPOperations<T, DIM>::LinearConstraint
    CBFQPOperations<T, DIM>::rightBorderConstraint(const Vector &state,
                                                   const Vector &target_state) {
        // Get right border constraint coefficients from the CBF implementation
        // These coefficients ensure the target stays within the right FOV boundary
        Vector coefficients = -1.0 * cbf_->getRBConstraints(state, target_state);

        // Get the right border constraint bound from the CBF implementation
        T bound = cbf_->getRBBound(state, target_state);

        // Create a linear constraint with lower bound negative infinity
        // This represents: coefficients^T * u <= bound
        return LinearConstraint(coefficients, std::numeric_limits<T>::lowest(), bound);
    }

    template <typename T, unsigned int DIM>
    typename CBFQPOperations<T, DIM>::LinearConstraint
    CBFQPOperations<T, DIM>::rangeConstraint(const Vector &state,
                                             const Vector &target_state) {
        // Get range constraint coefficients from the CBF implementation
        // These coefficients ensure the target stays within sensing range
        Vector coefficients = -1.0 * cbf_->getRangeConstraints(state, target_state);

        // Get the range constraint bound from the CBF implementation
        T bound = cbf_->getRangeBound(state, target_state);

        // Create a linear constraint with lower bound negative infinity
        // This represents: coefficients^T * u <= bound
        return LinearConstraint(coefficients, std::numeric_limits<T>::lowest(), bound);
    }

    template <typename T, unsigned int DIM>
    std::vector<typename CBFQPOperations<T, DIM>::LinearConstraint>
    CBFQPOperations<T, DIM>::minVelConstraints(const Vector &state) {
        // Get minimum velocity constraint coefficients for all dimensions
        Matrix coefficient_matrix = cbf_->getMinVelContraints(state);

        // Get minimum velocity constraint bounds for all dimensions
        Vector bounds = cbf_->getMinVelBounds(state);

        // Ensure the coefficient matrix and bounds vector dimensions match
        assert(coefficient_matrix.rows() == bounds.size());

        // Create a vector of linear constraints, one for each minimum velocity constraint
        std::vector<LinearConstraint> linear_constraints;
        for (size_t i = 0; i < coefficient_matrix.rows(); ++i) {
            // Extract the coefficients for this constraint
            Vector coefficients = -1.0 * coefficient_matrix.row(i);

            // Extract the bound for this constraint
            T bound = bounds(i);

            // Create a linear constraint with lower bound negative infinity
            // This represents: coefficients^T * u <= bound
            linear_constraints.push_back(LinearConstraint(coefficients, std::numeric_limits<T>::lowest(), bound));
        }
        return linear_constraints;
    }

    template <typename T, unsigned int DIM>
    std::vector<typename CBFQPOperations<T, DIM>::LinearConstraint>
    CBFQPOperations<T, DIM>::maxVelConstraints(const Vector &state) {
        // Get maximum velocity constraint coefficients for all dimensions
        Matrix coefficient_matrix = cbf_->getMaxVelContraints(state);

        // Get maximum velocity constraint bounds for all dimensions
        Vector bounds = cbf_->getMaxVelBounds(state);

        // Ensure the coefficient matrix and bounds vector dimensions match
        assert(coefficient_matrix.rows() == bounds.size());

        // Create a vector of linear constraints, one for each maximum velocity constraint
        std::vector<LinearConstraint> linear_constraints;
        for (size_t i = 0; i < coefficient_matrix.rows(); ++i) {
            // Extract the coefficients for this constraint
            Vector coefficients = -1.0 * coefficient_matrix.row(i);

            // Extract the bound for this constraint
            T bound = bounds(i);

            // Create a linear constraint with lower bound negative infinity
            // This represents: coefficients^T * u <= bound
            linear_constraints.push_back(LinearConstraint(coefficients, std::numeric_limits<T>::lowest(), bound));
        }
        return linear_constraints;
    }

    template <typename T, unsigned int DIM>
    typename CBFQPOperations<T, DIM>::LinearConstraint
    CBFQPOperations<T, DIM>::leftBorderConstraintWithSlackVar(const Vector &state,
                                                              const Vector &target_state,
                                                              const T &slack) {
        // Get left border constraint coefficients from the CBF implementation
        Vector coefficients = -1.0 * cbf_->getLBConstraints(state, target_state);

        // Get the left border constraint bound from the CBF implementation
        T bound = cbf_->getLBBound(state, target_state);

        // Create a linear constraint with lower bound negative infinity and add slack to the bound
        // This allows the constraint to be violated by up to 'slack' amount
        return LinearConstraint(coefficients, std::numeric_limits<T>::lowest(), bound+slack);
    }

    template <typename T, unsigned int DIM>
    typename CBFQPOperations<T, DIM>::LinearConstraint
    CBFQPOperations<T, DIM>::rightBorderConstraintWithSlackVar(const Vector &state,
                                                               const Vector &target_state,
                                                               const T &slack) {
        // Get right border constraint coefficients from the CBF implementation
        Vector coefficients = -1.0 * cbf_->getRBConstraints(state, target_state);

        // Get the right border constraint bound from the CBF implementation
        T bound = cbf_->getRBBound(state, target_state);

        // Create a linear constraint with lower bound negative infinity and add slack to the bound
        // This allows the constraint to be violated by up to 'slack' amount
        return LinearConstraint(coefficients, std::numeric_limits<T>::lowest(), bound+slack);
    }

    template <typename T, unsigned int DIM>
    typename CBFQPOperations<T, DIM>::LinearConstraint
    CBFQPOperations<T, DIM>::rangeConstraintWithSlackVar(const Vector &state,
                                                         const Vector &target_state,
                                                         const T &slack) {
        // Get range constraint coefficients from the CBF implementation
        Vector coefficients = -1.0 * cbf_->getRangeConstraints(state, target_state);

        // Get the range constraint bound from the CBF implementation
        T bound = cbf_->getRangeBound(state, target_state);

        // Create a linear constraint with lower bound negative infinity
        // Note: For some reason slack is not added to the bound in this case, unlike the border constraints
        return LinearConstraint(coefficients, std::numeric_limits<T>::lowest(), bound);
    }

    template <typename T, unsigned int DIM>
    typename CBFQPOperations<T, DIM>::DecisionVariableBounds
    CBFQPOperations<T, DIM>::controlBoundConstraint(const VectorDIM &u_min,
                                                    const VectorDIM &u_max) {
        // Initialize vectors for the lower and upper bounds on control inputs
        VectorDIM lower_bounds;
        VectorDIM upper_bounds;

        // Set the individual bounds for each control dimension
        for (size_t d = 0; d < DIM; ++d) {
            lower_bounds(d) = u_min(d);
            upper_bounds(d) = u_max(d);
        }

        // Return a structure containing the bounds
        return DecisionVariableBounds(lower_bounds, upper_bounds);
    }

    // Explicitly instantiate the template class for double and 3D
    template class CBFQPOperations<double, 3U>;

} // cbf