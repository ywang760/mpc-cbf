//
// Created by yutong on 8/4/25.
//

#include <mpc_cbf/optimization/ConnectivityMPCCBFQPOperations.h>

namespace mpc_cbf {
template <typename T, unsigned int DIM>
ConnectivityMPCCBFQPOperations<T, DIM>::ConnectivityMPCCBFQPOperations(
    Params& p, std::shared_ptr<DoubleIntegrator> model_ptr,
    std::shared_ptr<ConnectivityCBF> connectivity_cbf_ptr)
    : MPCCBFQPOperationsBase<T, DIM>(model_ptr), connectivity_cbf_ptr_(connectivity_cbf_ptr) {
    // mpc operations
    typename PiecewiseBezierMPCQPOperations::Params bezier_mpc_p = {p.piecewise_bezier_params,
                                                                    p.mpc_params};
    this->piecewise_mpc_operations_ptr_ =
        std::make_unique<PiecewiseBezierMPCQPOperations>(bezier_mpc_p, model_ptr);

    // mpc params
    this->h_ = p.mpc_params.h_;
    this->k_hor_ = p.mpc_params.k_hor_;
    this->mpc_tuning_ = p.mpc_params.tuning_;
    // control input predict
    this->U_basis_ =
        this->piecewise_mpc_operations_ptr_->U_basis(); // [3K, num_piece*dim*num_control_pts]
}

template <typename T, unsigned int DIM>
typename ConnectivityMPCCBFQPOperations<T, DIM>::LinearConstraint
ConnectivityMPCCBFQPOperations<T, DIM>::safetyCBFConstraint(const Vector& current_state,
                                                            const Vector& neighbor_state,
                                                            T slack_value) {
    // Use connectivity_cbf_ptr_ for proper constraint generation
    Vector a = connectivity_cbf_ptr_->getSafetyConstraints(current_state, neighbor_state);
    T b = connectivity_cbf_ptr_->getSafetyBound(current_state, neighbor_state);

    Vector A0 = Vector::Zero(DIM * this->k_hor_); // [3K, 1]
    A0.segment(0, DIM) = a;
    Row A_control_pts = -1.0 * A0.transpose() * this->U_basis_;

    return LinearConstraint(A_control_pts, std::numeric_limits<T>::lowest(), b + slack_value);
}

template <typename T, unsigned int DIM>
typename ConnectivityMPCCBFQPOperations<T, DIM>::LinearConstraint
ConnectivityMPCCBFQPOperations<T, DIM>::connectivityConstraint(const Vector& current_state,
                                                               const Vector& neighbor_state,
                                                               T slack_value) {
    // FIXME: this is not working
    // For connectivity constraint, implement a simple distance-based upper
    // bound Since getConnConstraints requires robot_states and eigenvec not
    // available here, we use a simplified approach based on maximum
    // distance constraint
    Vector a = Vector::Zero(DIM);
    VectorDIM current_pos = current_state.segment(0, DIM);
    VectorDIM neighbor_pos = neighbor_state.segment(0, DIM);
    VectorDIM relative_pos = current_pos - neighbor_pos;
    T distance = relative_pos.norm();
    if (distance > 0.01) {                            // Avoid division by zero
        a.segment(0, DIM) = -relative_pos / distance; // Negative for upper bound constraint
    }

    T b = 10.0; // Connectivity bound (distance <= 10.0) - should match
                // dmax parameter
    Vector A0 = Vector::Zero(DIM * this->k_hor_);
    A0.segment(0, DIM) = a;
    Row A_control_pts = -1.0 * A0.transpose() * this->U_basis_;

    return LinearConstraint(A_control_pts, std::numeric_limits<T>::lowest(), b + slack_value);
}

template <typename T, unsigned int DIM>
std::vector<typename ConnectivityMPCCBFQPOperations<T, DIM>::LinearConstraint>
ConnectivityMPCCBFQPOperations<T, DIM>::predSafetyCBFConstraints(
    const std::vector<State>& pred_states, const Vector& neighbor_state) {
    std::vector<LinearConstraint> linear_constraints;
    for (size_t k = 0; k < pred_states.size(); ++k) {
        const State& pred_state = pred_states.at(k);
        Vector state(2 * DIM);
        state << pred_state.pos_, pred_state.vel_;

        // Use connectivity_cbf_ptr_ for proper constraint generation
        Vector ak = connectivity_cbf_ptr_->getSafetyConstraints(state, neighbor_state);
        T bk = connectivity_cbf_ptr_->getSafetyBound(state, neighbor_state);

        Vector Ak = Vector::Zero(DIM * this->k_hor_);
        Ak.segment(k * DIM, DIM) = ak;
        Row Ak_control_pts = -1.0 * Ak.transpose() * this->U_basis_;
        linear_constraints.push_back(
            LinearConstraint(Ak_control_pts, std::numeric_limits<T>::lowest(), bk));
    }
    return linear_constraints;
}

template <typename T, unsigned int DIM>
std::vector<typename ConnectivityMPCCBFQPOperations<T, DIM>::LinearConstraint>
ConnectivityMPCCBFQPOperations<T, DIM>::predConnectivityConstraints(
    const std::vector<State>& pred_states, const Vector& neighbor_state) {
    std::vector<LinearConstraint> linear_constraints;
    for (size_t k = 0; k < pred_states.size(); ++k) {
        const State& pred_state = pred_states.at(k);

        // FIXME: this is not working
        // For connectivity constraint, implement a simple distance-based
        // upper bound Since getConnConstraints requires robot_states and
        // eigenvec not available here, we use a simplified approach based
        // on maximum distance constraint
        Vector ak = Vector::Zero(DIM);
        VectorDIM neighbor_pos = neighbor_state.segment(0, DIM);
        VectorDIM relative_pos = pred_state.pos_ - neighbor_pos;
        T distance = relative_pos.norm();
        if (distance > 0.01) {
            ak.segment(0, DIM) = -relative_pos / distance; // Negative for upper bound
        }
        T bk = 10.0; // Connectivity bound (distance <= 10.0) - should
                     // match dmax parameter

        Vector Ak = Vector::Zero(DIM * this->k_hor_);
        Ak.segment(k * DIM, DIM) = ak;
        Row Ak_control_pts = -1.0 * Ak.transpose() * this->U_basis_;
        linear_constraints.push_back(
            LinearConstraint(Ak_control_pts, std::numeric_limits<T>::lowest(), bk));
    }
    return linear_constraints;
}

template <typename T, unsigned int DIM>
std::unique_ptr<typename ConnectivityMPCCBFQPOperations<T, DIM>::PiecewiseBezierMPCQPOperations>
ConnectivityMPCCBFQPOperations<T, DIM>::piecewise_mpc_operations_ptr() {
    return std::move(this->piecewise_mpc_operations_ptr_);
}

// Explicit template instantiation
template class ConnectivityMPCCBFQPOperations<double, 3U>;

} // namespace mpc_cbf
