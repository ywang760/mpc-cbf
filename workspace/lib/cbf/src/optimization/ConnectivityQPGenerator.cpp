//
// Created on 4/10/25.
//

#include <cbf/optimization/ConnectivityQPGenerator.h>

namespace cbf {
    template <typename T, unsigned int DIM>
    ConnectivityQPGenerator<T, DIM>::ConnectivityQPGenerator(std::shared_ptr<ConnectivityCBF> cbf, int num_neighbors, bool slack_mode) 
        : CBFQPGeneratorBase<T, DIM>(num_neighbors, slack_mode), cbf_(cbf) {
    }

    // TODO: Implement the addConnectivityConstraint method
    // template <typename T, unsigned int DIM>
    // void ConnectivityQPGenerator<T, DIM>::addConnectivityConstraint(const Vector &state,
    //                                                 const Vector &target_state,
    //                                                 bool use_slack,
    //                                                 std::size_t slack_idx) {
    //     // Get connectivity constraint coefficients from the CBF implementation
    //     Vector coefficients = -1.0 * cbf_->getConnectivityConstraints(state, target_state);

    //     // Get the connectivity constraint bound from the CBF implementation
    //     T bound = cbf_->getConnectivityBound(state, target_state);

    //     // Create a linear constraint with lower bound negative infinity
    //     LinearConstraint linear_constraint(coefficients, std::numeric_limits<T>::lowest(), bound);
        
    //     // Handle constraint with or without slack variable
    //     if (use_slack && !this->slack_variables_.empty()) {
    //         // Create a row vector for slack variable coefficients (all zeros except at slack_idx)
    //         Row slack_coefficients = Row::Zero(this->slack_variables_.size());
    //         slack_coefficients(slack_idx) = -1; // Negative coefficient allows constraint relaxation
            
    //         // Add the constraint with slack variable to the QP problem
    //         this->addLinearConstraintForControlInputWithSlackVariables(linear_constraint, slack_coefficients);
    //     } else {
    //         // Add the constraint without slack variable to the QP problem
    //         this->addLinearConstraintForControlInput(linear_constraint);
    //     }
    // }

    template <typename T, unsigned int DIM>
    void ConnectivityQPGenerator<T, DIM>::addSafetyConstraint(const Vector &state,
                                                   const Vector &target_state,
                                                   bool use_slack,
                                                   std::size_t slack_idx) {
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
            this->addLinearConstraintForControlInputWithSlackVariables(linear_constraint, slack_coefficients);
        } else {
            // Add the constraint without slack variable to the QP problem
            this->addLinearConstraintForControlInput(linear_constraint);
        }
    }

    template <typename T, unsigned int DIM>
    void ConnectivityQPGenerator<T, DIM>::addMinVelConstraints(const Vector &state) {
        // Get minimum velocity constraint coefficients for all dimensions
        typename CBFQPGeneratorBase<T, DIM>::Matrix coefficient_matrix = cbf_->getMinVelContraints(state);

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
    void ConnectivityQPGenerator<T, DIM>::addMaxVelConstraints(const Vector &state) {
        // Get maximum velocity constraint coefficients for all dimensions
        typename CBFQPGeneratorBase<T, DIM>::Matrix coefficient_matrix = cbf_->getMaxVelContraints(state);

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
    template class ConnectivityQPGenerator<double, 3U>;

} // cbf