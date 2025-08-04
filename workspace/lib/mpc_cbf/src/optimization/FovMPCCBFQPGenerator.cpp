//
// Created by lishuo on 9/21/24.
//

#include <mpc_cbf/optimization/FovMPCCBFQPGenerator.h>

namespace mpc_cbf {
template <typename T, unsigned int DIM>
void FovMPCCBFQPGenerator<T, DIM>::addPiecewise(
    std::unique_ptr<FovMPCCBFQPOperations>&& piecewise_mpc_cbf_operations_ptr, int num_neighbors,
    bool slack_mode) {
    // init the PiecewiseBezierMPCQPGenerator API
    std::unique_ptr<PiecewiseBezierMPCQPOperations> piecewise_mpc_operations_ptr =
        piecewise_mpc_cbf_operations_ptr->piecewise_mpc_operations_ptr();
    this->piecewise_mpc_qp_generator_ptr_->addPiecewise(std::move(piecewise_mpc_operations_ptr));
    piecewise_mpc_cbf_operations_ptr_ = std::move(piecewise_mpc_cbf_operations_ptr);
    // add slack variable to the problem
    if (slack_mode) {
        this->addSlackVariables(num_neighbors);
    }
}

template <typename T, unsigned int DIM>
void FovMPCCBFQPGenerator<T, DIM>::addSafetyCBFConstraint(const State& current_state,
                                                          const Vector& other_pos, T slack_value) {
    LinearConstraint linear_constraint = piecewise_mpc_cbf_operations_ptr_->safetyCBFConstraint(
        current_state, other_pos, slack_value);
    this->piecewise_mpc_qp_generator_ptr_->addLinearConstraintForPiecewise(linear_constraint);
}

template <typename T, unsigned int DIM>
void FovMPCCBFQPGenerator<T, DIM>::addFovLBConstraint(const State& current_state,
                                                      const Vector& other_pos, T slack_value) {
    std::vector<LinearConstraint> linear_constraints =
        piecewise_mpc_cbf_operations_ptr_->fovLBConstraint(current_state, other_pos, slack_value);
    for (size_t i = 0; i < linear_constraints.size(); ++i) {
        this->piecewise_mpc_qp_generator_ptr_->addLinearConstraintForPiecewise(
            linear_constraints.at(i));
    }
}

template <typename T, unsigned int DIM>
void FovMPCCBFQPGenerator<T, DIM>::addFovRBConstraint(const State& current_state,
                                                      const Vector& other_pos, T slack_value) {
    std::vector<LinearConstraint> linear_constraints =
        piecewise_mpc_cbf_operations_ptr_->fovRBConstraint(current_state, other_pos, slack_value);
    for (size_t i = 0; i < linear_constraints.size(); ++i) {
        this->piecewise_mpc_qp_generator_ptr_->addLinearConstraintForPiecewise(
            linear_constraints.at(i));
    }
}

template <typename T, unsigned int DIM>
void FovMPCCBFQPGenerator<T, DIM>::addRangeCBFConstraint(const State& current_state,
                                                         const Vector& other_pos, T slack_value) {
    LinearConstraint linear_constraint = piecewise_mpc_cbf_operations_ptr_->rangeCBFConstraint(
        current_state, other_pos, slack_value);
    this->piecewise_mpc_qp_generator_ptr_->addLinearConstraintForPiecewise(linear_constraint);
}

template <typename T, unsigned int DIM>
void FovMPCCBFQPGenerator<T, DIM>::addPredSafetyCBFConstraints(
    const std::vector<State>& pred_states, const Vector& other_pos,
    const std::vector<T>& slack_values) {
    std::vector<LinearConstraint> linear_constraints =
        piecewise_mpc_cbf_operations_ptr_->predSafetyCBFConstraints(pred_states, other_pos);
    for (size_t i = 0; i < linear_constraints.size(); ++i) {
        this->piecewise_mpc_qp_generator_ptr_->addLinearConstraintForPiecewise(
            linear_constraints.at(i));
    }
}

template <typename T, unsigned int DIM>
void FovMPCCBFQPGenerator<T, DIM>::addPredFovLBConstraints(const std::vector<State>& pred_states,
                                                           const Vector& other_pos,
                                                           const std::vector<T>& slack_values) {
    std::vector<LinearConstraint> linear_constraints =
        piecewise_mpc_cbf_operations_ptr_->predFovLBConstraints(pred_states, other_pos);
    for (size_t i = 0; i < linear_constraints.size(); ++i) {
        this->piecewise_mpc_qp_generator_ptr_->addLinearConstraintForPiecewise(
            linear_constraints.at(i));
    }
}

template <typename T, unsigned int DIM>
void FovMPCCBFQPGenerator<T, DIM>::addPredFovRBConstraints(const std::vector<State>& pred_states,
                                                           const Vector& other_pos,
                                                           const std::vector<T>& slack_values) {
    std::vector<LinearConstraint> linear_constraints =
        piecewise_mpc_cbf_operations_ptr_->predFovRBConstraints(pred_states, other_pos);
    for (size_t i = 0; i < linear_constraints.size(); ++i) {
        this->piecewise_mpc_qp_generator_ptr_->addLinearConstraintForPiecewise(
            linear_constraints.at(i));
    }
}

template <typename T, unsigned int DIM>
void FovMPCCBFQPGenerator<T, DIM>::addPredRangeCBFConstraints(const std::vector<State>& pred_states,
                                                              const Vector& other_pos,
                                                              const std::vector<T>& slack_values) {
    std::vector<LinearConstraint> linear_constraints =
        piecewise_mpc_cbf_operations_ptr_->predRangeCBFConstraints(pred_states, other_pos);
    for (size_t i = 0; i < linear_constraints.size(); ++i) {
        this->piecewise_mpc_qp_generator_ptr_->addLinearConstraintForPiecewise(
            linear_constraints.at(i));
    }
}

template <typename T, unsigned int DIM>
void FovMPCCBFQPGenerator<T, DIM>::addSafetyCBFConstraintWithSlackVariables(
    const State& current_state, const Vector& other_pos, std::size_t neighbor_idx) {
    LinearConstraint linear_constraint =
        piecewise_mpc_cbf_operations_ptr_->safetyCBFConstraint(current_state, other_pos);
    Row slack_coefficients = Row::Zero(this->slack_variables_.size());
    slack_coefficients(neighbor_idx) = -1;
    this->addLinearConstraintForPiecewiseWithSlackVariables(linear_constraint, slack_coefficients);
}

template <typename T, unsigned int DIM>
void FovMPCCBFQPGenerator<T, DIM>::addFovLBConstraintWithSlackVariables(const State& current_state,
                                                                        const Vector& other_pos,
                                                                        std::size_t neighbor_idx) {
    std::vector<LinearConstraint> linear_constraints =
        piecewise_mpc_cbf_operations_ptr_->fovLBConstraint(current_state, other_pos);
    Row slack_coefficients = Row::Zero(this->slack_variables_.size());
    slack_coefficients(neighbor_idx) = -1;
    for (size_t i = 0; i < linear_constraints.size(); ++i) {
        this->addLinearConstraintForPiecewiseWithSlackVariables(linear_constraints.at(i),
                                                                slack_coefficients);
    }
}

template <typename T, unsigned int DIM>
void FovMPCCBFQPGenerator<T, DIM>::addFovRBConstraintWithSlackVariables(const State& current_state,
                                                                        const Vector& other_pos,
                                                                        std::size_t neighbor_idx) {
    std::vector<LinearConstraint> linear_constraints =
        piecewise_mpc_cbf_operations_ptr_->fovRBConstraint(current_state, other_pos);
    Row slack_coefficients = Row::Zero(this->slack_variables_.size());
    slack_coefficients(neighbor_idx) = -1;
    for (size_t i = 0; i < linear_constraints.size(); ++i) {
        this->addLinearConstraintForPiecewiseWithSlackVariables(linear_constraints.at(i),
                                                                slack_coefficients);
    }
}

template <typename T, unsigned int DIM>
void FovMPCCBFQPGenerator<T, DIM>::addRangeCBFConstraintWithSlackVariables(
    const State& current_state, const Vector& other_pos, std::size_t neighbor_idx) {
    LinearConstraint linear_constraint =
        piecewise_mpc_cbf_operations_ptr_->rangeCBFConstraint(current_state, other_pos);
    Row slack_coefficients = Row::Zero(this->slack_variables_.size());
    slack_coefficients(neighbor_idx) = -1;
    this->addLinearConstraintForPiecewiseWithSlackVariables(linear_constraint, slack_coefficients);
}

template <typename T, unsigned int DIM>
void FovMPCCBFQPGenerator<T, DIM>::addPredSafetyCBFConstraintsWithSlackVariables(
    const std::vector<State>& pred_states, const Vector& other_pos, std::size_t neighbor_idx) {
    std::vector<LinearConstraint> linear_constraints =
        piecewise_mpc_cbf_operations_ptr_->predSafetyCBFConstraints(pred_states, other_pos);
    Row slack_coefficients = Row::Zero(this->slack_variables_.size());
    slack_coefficients(neighbor_idx) = -1;
    for (size_t i = 0; i < linear_constraints.size(); ++i) {
        this->addLinearConstraintForPiecewiseWithSlackVariables(linear_constraints.at(i),
                                                                slack_coefficients);
    }
}

template <typename T, unsigned int DIM>
void FovMPCCBFQPGenerator<T, DIM>::addPredFovLBConstraintsWithSlackVariables(
    const std::vector<State>& pred_states, const Vector& other_pos, std::size_t neighbor_idx) {
    std::vector<LinearConstraint> linear_constraints =
        piecewise_mpc_cbf_operations_ptr_->predFovLBConstraints(pred_states, other_pos);
    Row slack_coefficients = Row::Zero(this->slack_variables_.size());
    slack_coefficients(neighbor_idx) = -1;
    for (size_t i = 0; i < linear_constraints.size(); ++i) {
        this->addLinearConstraintForPiecewiseWithSlackVariables(linear_constraints.at(i),
                                                                slack_coefficients);
    }
}

template <typename T, unsigned int DIM>
void FovMPCCBFQPGenerator<T, DIM>::addPredFovRBConstraintsWithSlackVariables(
    const std::vector<State>& pred_states, const Vector& other_pos, std::size_t neighbor_idx) {
    std::vector<LinearConstraint> linear_constraints =
        piecewise_mpc_cbf_operations_ptr_->predFovRBConstraints(pred_states, other_pos);
    Row slack_coefficients = Row::Zero(this->slack_variables_.size());
    slack_coefficients(neighbor_idx) = -1;
    for (size_t i = 0; i < linear_constraints.size(); ++i) {
        this->addLinearConstraintForPiecewiseWithSlackVariables(linear_constraints.at(i),
                                                                slack_coefficients);
    }
}

template <typename T, unsigned int DIM>
void FovMPCCBFQPGenerator<T, DIM>::addPredRangeCBFConstraintsWithSlackVariables(
    const std::vector<State>& pred_states, const Vector& other_pos, std::size_t neighbor_idx) {
    std::vector<LinearConstraint> linear_constraints =
        piecewise_mpc_cbf_operations_ptr_->predRangeCBFConstraints(pred_states, other_pos);
    Row slack_coefficients = Row::Zero(this->slack_variables_.size());
    slack_coefficients(neighbor_idx) = -1;
    for (size_t i = 0; i < linear_constraints.size(); ++i) {
        this->addLinearConstraintForPiecewiseWithSlackVariables(linear_constraints.at(i),
                                                                slack_coefficients);
    }
}

template class FovMPCCBFQPGenerator<double, 3U>;

} // namespace mpc_cbf