//
// Created by yutong on 8/4/25.
//
#include <spdlog/spdlog.h>
#include <mpc_cbf/optimization/ConnectivityMPCCBFQPGenerator.h>
#include <limits>

namespace mpc_cbf {
template <typename T, unsigned int DIM>
ConnectivityMPCCBFQPGenerator<T, DIM>::ConnectivityMPCCBFQPGenerator(
    std::unique_ptr<ConnectivityMPCCBFQPOperations>&& piecewise_mpc_cbf_operations_ptr,
    int num_neighbors, const cbf::SlackConfig& slack_config) {
    // init the PiecewiseBezierMPCQPGenerator API
    std::unique_ptr<PiecewiseBezierMPCQPOperations> piecewise_mpc_operations_ptr =
        piecewise_mpc_cbf_operations_ptr->piecewise_mpc_operations_ptr();
    this->piecewise_mpc_qp_generator_ptr_->addPiecewise(std::move(piecewise_mpc_operations_ptr));
    // setup the fields for ConnectivityMPCCBFQPGenerator
    piecewise_mpc_cbf_operations_ptr_ = std::move(piecewise_mpc_cbf_operations_ptr);
    slack_config_ = slack_config;
    
    // Initialize separate slack variable pools based on configuration
    if (slack_config_.safety_slack) {
        for (int i = 0; i < num_neighbors; ++i) {
            auto* var = this->problem().addVariable(0, std::numeric_limits<T>::max());
            safety_slack_variables_.push_back(var);
        }
    }
    
    if (slack_config_.clf_slack) {
        for (int i = 0; i < num_neighbors; ++i) {
            auto* var = this->problem().addVariable(0, std::numeric_limits<T>::max());
            clf_slack_variables_.push_back(var);
        }
    }
    
    if (slack_config_.connectivity_slack) {
        auto* var = this->problem().addVariable(0, std::numeric_limits<T>::max());
        connectivity_slack_variables_.push_back(var);
    }
    
    // Legacy slack variables no longer needed - using separate pools
}

template <typename T, unsigned int DIM>
void ConnectivityMPCCBFQPGenerator<T, DIM>::addSafetyCBFConstraint(const Vector& current_state,
                                                                   const Vector& neighbor_state,
                                                                   std::size_t neighbor_idx,
                                                                   T slack_value) {
    LinearConstraint linear_constraint = piecewise_mpc_cbf_operations_ptr_->safetyCBFConstraint(
        current_state, neighbor_state, slack_value);

    if (slack_config_.safety_slack && !safety_slack_variables_.empty() &&
        neighbor_idx < safety_slack_variables_.size()) {
        // Use safety-specific slack variables
        Row slack_coefficients = Row::Zero(safety_slack_variables_.size());
        slack_coefficients(neighbor_idx) = -1.0;
        this->addLinearConstraintForPiecewiseWithSlackVariables(
            linear_constraint, slack_coefficients, safety_slack_variables_);
    } else {
        // No slack variables
        this->piecewise_mpc_qp_generator_ptr_->addLinearConstraintForPiecewise(linear_constraint);
    }
}

template <typename T, unsigned int DIM>
void ConnectivityMPCCBFQPGenerator<T, DIM>::addConnectivityConstraint(
    const Vector& x_self,
    const std::vector<VectorDIM>& other_positions,
    T slack_value) {
    LinearConstraint linear_constraint = piecewise_mpc_cbf_operations_ptr_->connectivityConstraint(
        x_self, other_positions, slack_value);
        if (slack_config_.connectivity_slack && !connectivity_slack_variables_.empty()) {
            // Use connectivity-specific slack variables
            Row slack_coefficients = Row::Zero(connectivity_slack_variables_.size());
            slack_coefficients(0) = -1.0;  // Only one connectivity slack variable
            this->addLinearConstraintForPiecewiseWithSlackVariables(
                linear_constraint, slack_coefficients, connectivity_slack_variables_);
        } else {
            // No slack variables
            this->piecewise_mpc_qp_generator_ptr_->addLinearConstraintForPiecewise(linear_constraint);
        }
    }

template <typename T, unsigned int DIM>
void ConnectivityMPCCBFQPGenerator<T, DIM>::addCLFConstraint(const Vector& current_state,
                                                             const Vector& neighbor_state,
                                                             std::size_t neighbor_idx,
                                                             T slack_value) {
    LinearConstraint linear_constraint = piecewise_mpc_cbf_operations_ptr_->clfConstraint(
        current_state, neighbor_state, slack_value);

    if (slack_config_.clf_slack && !clf_slack_variables_.empty() &&
        neighbor_idx < clf_slack_variables_.size()) {
        // Use CLF-specific slack variables
        Row slack_coefficients = Row::Zero(clf_slack_variables_.size());
        slack_coefficients(neighbor_idx) = -1.0;
        this->addLinearConstraintForPiecewiseWithSlackVariables(
            linear_constraint, slack_coefficients, clf_slack_variables_);
    } else {
        // No slack variables
        this->piecewise_mpc_qp_generator_ptr_->addLinearConstraintForPiecewise(linear_constraint);
    }
}

template <typename T, unsigned int DIM>
void ConnectivityMPCCBFQPGenerator<T, DIM>::addPredSafetyCBFConstraints(
    const std::vector<State>& pred_states, const Vector& neighbor_state, std::size_t neighbor_idx) {
    std::vector<LinearConstraint> linear_constraints =
        piecewise_mpc_cbf_operations_ptr_->predSafetyCBFConstraints(pred_states, neighbor_state);

    if (slack_config_.safety_slack && !safety_slack_variables_.empty() &&
        neighbor_idx < safety_slack_variables_.size()) {
        // Use safety-specific slack variables
        Row slack_coefficients = Row::Zero(safety_slack_variables_.size());
        slack_coefficients(neighbor_idx) = -1.0;

        for (const auto& constraint : linear_constraints) {
            this->addLinearConstraintForPiecewiseWithSlackVariables(
                constraint, slack_coefficients, safety_slack_variables_);
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
    const std::vector<State>& pred_states, const std::vector<VectorDIM>& other_positions,
    T slack_value) {
    std::vector<LinearConstraint> linear_constraints =
        piecewise_mpc_cbf_operations_ptr_->predConnectivityConstraints(pred_states, other_positions);
        if (slack_config_.connectivity_slack && !connectivity_slack_variables_.empty()) {
            // Use connectivity-specific slack variables
            Row slack_coefficients = Row::Zero(connectivity_slack_variables_.size());
            slack_coefficients(0) = -1.0;  // Only one connectivity slack variable
            for (size_t i = 0; i < linear_constraints.size(); ++i) {
                this->addLinearConstraintForPiecewiseWithSlackVariables(
                    linear_constraints.at(i), slack_coefficients, connectivity_slack_variables_);
            }
        } else {
            for (size_t i = 0; i < linear_constraints.size(); ++i) {
                this->piecewise_mpc_qp_generator_ptr_->addLinearConstraintForPiecewise(
                    linear_constraints.at(i));
            }
        }
}

template <typename T, unsigned int DIM>
void ConnectivityMPCCBFQPGenerator<T, DIM>::addPredCLFConstraints(
    const std::vector<State>& pred_states, const Vector& neighbor_state, std::size_t neighbor_idx) {
    std::vector<LinearConstraint> linear_constraints =
        piecewise_mpc_cbf_operations_ptr_->predCLFConstraints(pred_states, neighbor_state);

    if (slack_config_.clf_slack && !clf_slack_variables_.empty() &&
        neighbor_idx < clf_slack_variables_.size()) {
        // Use CLF-specific slack variables
        Row slack_coefficients = Row::Zero(clf_slack_variables_.size());
        slack_coefficients(neighbor_idx) = -1.0;

        for (const auto& constraint : linear_constraints) {
            this->addLinearConstraintForPiecewiseWithSlackVariables(
                constraint, slack_coefficients, clf_slack_variables_);
        }
    } else {
        for (size_t i = 0; i < linear_constraints.size(); ++i) {
            this->piecewise_mpc_qp_generator_ptr_->addLinearConstraintForPiecewise(
                linear_constraints.at(i));
        }
    }
}


template <typename T, unsigned int DIM>
std::shared_ptr<
    typename ConnectivityMPCCBFQPGenerator<T, DIM>::ConnectivityMPCCBFQPOperations::ConnectivityCBF>
ConnectivityMPCCBFQPGenerator<T, DIM>::connectivityCBF() const {
    return piecewise_mpc_cbf_operations_ptr_->connectivityCBF();
}

// Explicit template instantiation
template class ConnectivityMPCCBFQPGenerator<double, 3U>;

} // namespace mpc_cbf
