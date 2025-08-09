//
// Created by yutong on 8/4/25.
//

#include <mpc_cbf/optimization/ConnectivityMPCCBFQPGenerator.h>

namespace mpc_cbf {
template <typename T, unsigned int DIM>
ConnectivityMPCCBFQPGenerator<T, DIM>::ConnectivityMPCCBFQPGenerator(
    std::unique_ptr<ConnectivityMPCCBFQPOperations>&& piecewise_mpc_cbf_operations_ptr,
    int num_neighbors, bool slack_mode) {
    // init the PiecewiseBezierMPCQPGenerator API
    std::unique_ptr<PiecewiseBezierMPCQPOperations> piecewise_mpc_operations_ptr =
        piecewise_mpc_cbf_operations_ptr->piecewise_mpc_operations_ptr();
    this->piecewise_mpc_qp_generator_ptr_->addPiecewise(std::move(piecewise_mpc_operations_ptr));
    // setup the fields for ConnectivityMPCCBFQPGenerator
    piecewise_mpc_cbf_operations_ptr_ = std::move(piecewise_mpc_cbf_operations_ptr);
    slack_mode_ = slack_mode;
    // add slack variable to the problem
    if (slack_mode_) {
        this->addSlackVariables(num_neighbors);
    }
}

template <typename T, unsigned int DIM>
void ConnectivityMPCCBFQPGenerator<T, DIM>::addSafetyCBFConstraint(const Vector& current_state,
                                                                   const Vector& neighbor_state,
                                                                   std::size_t neighbor_idx,
                                                                   T slack_value) {
    LinearConstraint linear_constraint = piecewise_mpc_cbf_operations_ptr_->safetyCBFConstraint(
        current_state, neighbor_state, slack_value);

    if (this->slack_mode_) {
        // Create slack coefficient vector (only the neighbor_idx-th slack
        // variable has coefficient -1)
        Row slack_coefficients = Row::Zero(this->slack_variables_.size());
        slack_coefficients(neighbor_idx) = -1.0;

        this->addLinearConstraintForPiecewiseWithSlackVariables(linear_constraint,
                                                                slack_coefficients);
    } else {
        // If not in slack mode, just add the linear constraint directly
        this->piecewise_mpc_qp_generator_ptr_->addLinearConstraintForPiecewise(linear_constraint);
    }
}

template <typename T, unsigned int DIM>
void ConnectivityMPCCBFQPGenerator<T, DIM>::addConnectivityConstraint(
    const Eigen::MatrixXd& robot_states, size_t self_idx, T slack_value) {
    LinearConstraint linear_constraint = piecewise_mpc_cbf_operations_ptr_->connectivityConstraint(
        robot_states, self_idx, slack_value);

    if (this->slack_mode_) {
        // Create slack coefficient vector (only the self_idx-th slack
        // variable has coefficient -1)
        Row slack_coefficients = Row::Zero(this->slack_variables_.size());
        slack_coefficients(self_idx) = -1.0;

        this->addLinearConstraintForPiecewiseWithSlackVariables(linear_constraint,
                                                                slack_coefficients);
    } else {
        this->piecewise_mpc_qp_generator_ptr_->addLinearConstraintForPiecewise(linear_constraint);
    }
}

template <typename T, unsigned int DIM>
void ConnectivityMPCCBFQPGenerator<T, DIM>::addPredSafetyCBFConstraints(
    const std::vector<State>& pred_states, const Vector& neighbor_state, std::size_t neighbor_idx) {
    std::vector<LinearConstraint> linear_constraints =
        piecewise_mpc_cbf_operations_ptr_->predSafetyCBFConstraints(pred_states, neighbor_state);

    if (this->slack_mode_) {
        // Create slack coefficient vector (only the neighbor_idx-th slack
        // variable has coefficient -1)
        Row slack_coefficients = Row::Zero(this->slack_variables_.size());
        slack_coefficients(neighbor_idx) = -1.0;

        for (const auto& constraint : linear_constraints) {
            this->addLinearConstraintForPiecewiseWithSlackVariables(constraint, slack_coefficients);
        }
    } else {
        for (size_t i = 0; i < linear_constraints.size(); ++i) {
            this->piecewise_mpc_qp_generator_ptr_->addLinearConstraintForPiecewise(
                linear_constraints.at(i));
        }
    }
}

template <typename T, unsigned int DIM>
void ConnectivityMPCCBFQPGenerator<T, DIM>::addPredConnectivityConstraints(
    const std::vector<State>& pred_states, const Eigen::MatrixXd& robot_states, size_t self_idx,
    const std::vector<T>& slack_values) {
    std::vector<LinearConstraint> linear_constraints =
        piecewise_mpc_cbf_operations_ptr_->predConnectivityConstraints(pred_states, robot_states,
                                                                       self_idx);
    if (this->slack_mode_) {
        // Create slack coefficient vector (only the self_idx-th slack
        // variable has coefficient -1)
        Row slack_coefficients = Row::Zero(this->slack_variables_.size());
        slack_coefficients(self_idx) = -1.0;
        for (size_t i = 0; i < linear_constraints.size(); ++i) {
            this->addLinearConstraintForPiecewiseWithSlackVariables(linear_constraints.at(i),
                                                                    slack_coefficients);
        }
    } else {
        for (size_t i = 0; i < linear_constraints.size(); ++i) {
            this->piecewise_mpc_qp_generator_ptr_->addLinearConstraintForPiecewise(
                linear_constraints.at(i));
        }
    }
}

// Explicit template instantiation
template class ConnectivityMPCCBFQPGenerator<double, 3U>;

} // namespace mpc_cbf
