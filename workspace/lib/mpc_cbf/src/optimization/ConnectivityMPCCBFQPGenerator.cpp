//
// Created by yutong on 8/4/25.
//

#include <mpc_cbf/optimization/ConnectivityMPCCBFQPGenerator.h>

namespace mpc_cbf {
    template <typename T, unsigned int DIM>
    void ConnectivityMPCCBFQPGenerator<T, DIM>::addPiecewise(
            std::unique_ptr<ConnectivityMPCCBFQPOperations> &&piecewise_mpc_cbf_operations_ptr, int num_neighbors, bool slack_mode) {
        // init the PiecewiseBezierMPCQPGenerator API
        std::unique_ptr<PiecewiseBezierMPCQPOperations> piecewise_mpc_operations_ptr = piecewise_mpc_cbf_operations_ptr->piecewise_mpc_operations_ptr();
        this->piecewise_mpc_qp_generator_ptr_->addPiecewise(std::move(piecewise_mpc_operations_ptr));
        piecewise_mpc_cbf_operations_ptr_ = std::move(piecewise_mpc_cbf_operations_ptr);
        // add slack variable to the problem
        if (slack_mode) {
            this->addSlackVariables(num_neighbors);
        }
    }

    template <typename T, unsigned int DIM>
    void ConnectivityMPCCBFQPGenerator<T, DIM>::addSafetyCBFConstraint(const State &current_state,
                                                                      const Vector &other_pos,
                                                                      T slack_value) {
        LinearConstraint linear_constraint = piecewise_mpc_cbf_operations_ptr_->safetyCBFConstraint(current_state, other_pos, slack_value);
        this->piecewise_mpc_qp_generator_ptr_->addLinearConstraintForPiecewise(linear_constraint);
    }

    template <typename T, unsigned int DIM>
    void ConnectivityMPCCBFQPGenerator<T, DIM>::addConnectivityConstraint(const State &current_state,
                                                                          const Vector &other_pos,
                                                                          T slack_value) {
        LinearConstraint linear_constraint = piecewise_mpc_cbf_operations_ptr_->connectivityConstraint(current_state, other_pos, slack_value);
        this->piecewise_mpc_qp_generator_ptr_->addLinearConstraintForPiecewise(linear_constraint);
    }

    template <typename T, unsigned int DIM>
    void ConnectivityMPCCBFQPGenerator<T, DIM>::addPredSafetyCBFConstraints(const std::vector<State> &pred_states,
                                                                           const Vector &other_pos,
                                                                           const std::vector<T>& slack_values) {
        std::vector<LinearConstraint> linear_constraints = piecewise_mpc_cbf_operations_ptr_->predSafetyCBFConstraints(pred_states, other_pos);
        for (size_t i = 0; i < linear_constraints.size(); ++i) {
            this->piecewise_mpc_qp_generator_ptr_->addLinearConstraintForPiecewise(linear_constraints.at(i));
        }
    }

    template <typename T, unsigned int DIM>
    void ConnectivityMPCCBFQPGenerator<T, DIM>::addPredConnectivityConstraints(const std::vector<State> &pred_states,
                                                                              const Vector &other_pos,
                                                                              const std::vector<T>& slack_values) {
        std::vector<LinearConstraint> linear_constraints = piecewise_mpc_cbf_operations_ptr_->predConnectivityConstraints(pred_states, other_pos);
        for (size_t i = 0; i < linear_constraints.size(); ++i) {
            this->piecewise_mpc_qp_generator_ptr_->addLinearConstraintForPiecewise(linear_constraints.at(i));
        }
    }

    template <typename T, unsigned int DIM>
    void ConnectivityMPCCBFQPGenerator<T, DIM>::addSafetyCBFConstraintWithSlackVariables(
            const State &current_state,
            const Vector &other_pos,
            std::size_t neighbor_idx) {
        LinearConstraint linear_constraint = piecewise_mpc_cbf_operations_ptr_->safetyCBFConstraint(current_state, other_pos);
        
        // Create slack coefficient vector (only the neighbor_idx-th slack variable has coefficient -1)
        Row slack_coefficients = Row::Zero(this->slack_variables_.size());
        slack_coefficients(neighbor_idx) = -1.0;
        
        this->addLinearConstraintForPiecewiseWithSlackVariables(linear_constraint, slack_coefficients);
    }

    template <typename T, unsigned int DIM>
    void ConnectivityMPCCBFQPGenerator<T, DIM>::addConnectivityConstraintWithSlackVariables(
            const State &current_state,
            const Vector &other_pos,
            std::size_t neighbor_idx) {
        LinearConstraint linear_constraint = piecewise_mpc_cbf_operations_ptr_->connectivityConstraint(current_state, other_pos);
        
        // Create slack coefficient vector (only the neighbor_idx-th slack variable has coefficient -1)
        Row slack_coefficients = Row::Zero(this->slack_variables_.size());
        slack_coefficients(neighbor_idx) = -1.0;
        
        this->addLinearConstraintForPiecewiseWithSlackVariables(linear_constraint, slack_coefficients);
    }

    template <typename T, unsigned int DIM>
    void ConnectivityMPCCBFQPGenerator<T, DIM>::addPredSafetyCBFConstraintsWithSlackVariables(
            const std::vector<State> &pred_states,
            const Vector &other_pos,
            std::size_t neighbor_idx) {
        std::vector<LinearConstraint> linear_constraints = piecewise_mpc_cbf_operations_ptr_->predSafetyCBFConstraints(pred_states, other_pos);
        
        // Create slack coefficient vector (only the neighbor_idx-th slack variable has coefficient -1)
        Row slack_coefficients = Row::Zero(this->slack_variables_.size());
        slack_coefficients(neighbor_idx) = -1.0;
        
        for (const auto& constraint : linear_constraints) {
            this->addLinearConstraintForPiecewiseWithSlackVariables(constraint, slack_coefficients);
        }
    }

    template <typename T, unsigned int DIM>
    void ConnectivityMPCCBFQPGenerator<T, DIM>::addPredConnectivityConstraintsWithSlackVariables(
            const std::vector<State> &pred_states,
            const Vector &other_pos,
            std::size_t neighbor_idx) {
        std::vector<LinearConstraint> linear_constraints = piecewise_mpc_cbf_operations_ptr_->predConnectivityConstraints(pred_states, other_pos);
        
        // Create slack coefficient vector (only the neighbor_idx-th slack variable has coefficient -1)
        Row slack_coefficients = Row::Zero(this->slack_variables_.size());
        slack_coefficients(neighbor_idx) = -1.0;
        
        for (const auto& constraint : linear_constraints) {
            this->addLinearConstraintForPiecewiseWithSlackVariables(constraint, slack_coefficients);
        }
    }

    // Explicit template instantiation
    template class ConnectivityMPCCBFQPGenerator<double, 3U>;

} // mpc_cbf
