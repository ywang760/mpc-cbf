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
ConnectivityMPCCBFQPOperations<T, DIM>::connectivityConstraint(const Eigen::MatrixXd& robot_states,
                                                               size_t self_idx, T slack_value) {
    // Extract current robot state
    Vector current_state = robot_states.row(self_idx);

    // Extract position information for lambda2 calculation
    const auto robot_positions = robot_states.leftCols(2); // Extract only position columns (x, y)

    // Calculate lambda2 and eigenvector for connectivity
    auto [lambda2_val, eigenvec] = connectivity_cbf_ptr_->getLambda2(robot_positions);
    double epsilon = 0.1;             // lambda2_min for connectivity CBF
    double h = lambda2_val - epsilon; // barrier function: h = λ₂ - λ₂_min

    // Initialize connectivity CBF with actual number of robots and self index
    const int num_robots = robot_states.rows();
    connectivity_cbf_ptr_->initConnCBF(num_robots, self_idx);

    // Get actual connectivity constraints from CBF
    Vector a = connectivity_cbf_ptr_->getConnConstraints(current_state, robot_states, eigenvec);
    T b = connectivity_cbf_ptr_->getConnBound(current_state, robot_states, eigenvec, h);

    Vector A0 = Vector::Zero(DIM * this->k_hor_);
    A0.segment(0, DIM) = a;
    Row A_control_pts = -1.0 * A0.transpose() * this->U_basis_;

    return LinearConstraint(A_control_pts, std::numeric_limits<T>::lowest(), b + slack_value);
}

template <typename T, unsigned int DIM>
typename ConnectivityMPCCBFQPOperations<T, DIM>::LinearConstraint
ConnectivityMPCCBFQPOperations<T, DIM>::clfConstraint(const Vector& current_state,
                                                      const Vector& neighbor_state, T slack_value) {
    Vector a = connectivity_cbf_ptr_->getCLFConstraints(current_state, neighbor_state);
    T b = -1.0 * connectivity_cbf_ptr_->getCLFBound(current_state, neighbor_state);

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
    const std::vector<State>& pred_states, const Eigen::MatrixXd& robot_states, size_t self_idx) {
    std::vector<LinearConstraint> linear_constraints;

    // Extract position information for lambda2 calculation (using current robot states)
    const auto robot_positions = robot_states.leftCols(2); // Extract only position columns (x, y)

    // Calculate lambda2 and eigenvector for connectivity
    auto [lambda2_val, eigenvec] = connectivity_cbf_ptr_->getLambda2(robot_positions);
    double epsilon = 0.1;             // lambda2_min for connectivity CBF
    double h = lambda2_val - epsilon; // barrier function: h = λ₂ - λ₂_min

    // Initialize connectivity CBF with actual number of robots and self index
    const int num_robots = robot_states.rows();
    connectivity_cbf_ptr_->initConnCBF(num_robots, self_idx);

    for (size_t k = 0; k < pred_states.size(); ++k) {
        const State& pred_state = pred_states.at(k);
        Vector current_pred_state(2 * DIM);
        current_pred_state << pred_state.pos_, pred_state.vel_;

        // Create updated robot_states matrix with predicted state for current robot
        Eigen::MatrixXd updated_robot_states = robot_states;
        updated_robot_states.row(self_idx) = current_pred_state.transpose();

        // Get actual connectivity constraints from CBF
        Vector ak = connectivity_cbf_ptr_->getConnConstraints(current_pred_state,
                                                              updated_robot_states, eigenvec);
        T bk = connectivity_cbf_ptr_->getConnBound(current_pred_state, updated_robot_states,
                                                   eigenvec, h);

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
ConnectivityMPCCBFQPOperations<T, DIM>::predCLFConstraints(const std::vector<State>& pred_states,
                                                           const Vector& neighbor_state) {
    std::vector<LinearConstraint> linear_constraints;
    for (size_t k = 0; k < pred_states.size(); ++k) {
        const State& pred_state = pred_states.at(k);
        Vector state(2 * DIM);
        state << pred_state.pos_, pred_state.vel_;

        Vector ak = connectivity_cbf_ptr_->getCLFConstraints(state, neighbor_state);
        T bk = -1.0 * connectivity_cbf_ptr_->getCLFBound(state, neighbor_state);

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

template <typename T, unsigned int DIM>
std::shared_ptr<typename ConnectivityMPCCBFQPOperations<T, DIM>::ConnectivityCBF>
ConnectivityMPCCBFQPOperations<T, DIM>::connectivityCBF() const {
    return connectivity_cbf_ptr_;
}

// Explicit template instantiation
template class ConnectivityMPCCBFQPOperations<double, 3U>;

} // namespace mpc_cbf
