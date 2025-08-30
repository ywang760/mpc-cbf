//
// Created by yutong on 8/4/25.
//

#include <mpc_cbf/optimization/ConnectivityMPCCBFQPOperations.h>
#include <spdlog/spdlog.h>

namespace mpc_cbf {
auto logger = spdlog::default_logger();
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
    // === 调试缓存 ===
    last_safety_a_ = a;
    last_safety_b_ = b;

    // === 可选：调试输出（建议只在 debug 级别打印） ===
    {
        const Eigen::IOFormat oneLine(Eigen::StreamPrecision, Eigen::DontAlignCols,
                                    " ", " ", "", "", "", "");
        std::ostringstream oss;
        oss << a.transpose().format(oneLine);
        logger->debug("[safetyCBFConstraint] a = [{}], b = {}", oss.str(), b);
    }                                                            
    Vector A0 = Vector::Zero(DIM * this->k_hor_); // [3K, 1]
    A0.segment(0, DIM) = a;
    Row A_control_pts = -1.0 * A0.transpose() * this->U_basis_;
    // logger->info("this->U_basis_: {}x{}", this->U_basis_.rows(), this->U_basis_.cols());
    // for (int r = 0; r < std::min(6, (int)this->U_basis_.rows()); ++r) {
    //     std::ostringstream oss;
    //     for (int c = 0; c < std::min(10, (int)this->U_basis_.cols()); ++c) oss << this->U_basis_(r, c) << " ";
    //     logger->info("[U] row{}: {}", r, oss.str());
    // }

    return LinearConstraint(A_control_pts, std::numeric_limits<T>::lowest(), b + slack_value);
}

template <typename T, unsigned int DIM>
typename ConnectivityMPCCBFQPOperations<T, DIM>::LinearConstraint
ConnectivityMPCCBFQPOperations<T, DIM>::connectivityConstraint(const Vector& x_self, const std::vector<VectorDIM>& other_positions, T slack_value)
{
    // === 1) 组装 robot_states: [self; others] ===
    const int N = 1 + static_cast<int>(other_positions.size());
    Eigen::MatrixXd robot_states(N, 6);
    robot_states.setZero();
    // self: [pos(DIM), vel(DIM)]
    for (unsigned j = 0; j < DIM; ++j) {
        robot_states(0, j)         = static_cast<double>(x_self(j));           // pos j
        robot_states(0, DIM + j)   = static_cast<double>(x_self(DIM + j));     // vel j
    }
    // 其余机器人的位置（只需 px, py）
    for (int i = 0; i < static_cast<int>(other_positions.size()); ++i) {
        robot_states.block(i + 1, 0, 1, DIM) =
            other_positions[i].head(DIM).template cast<double>().transpose();
    }    

    // === 2) 一次性从 CBF 拿回 Ac, Bc（内部已依据 robot_states/self_idx 计算 λ2、eigenvec、h 等）===
    auto [Ac, Bc] = connectivity_cbf_ptr_->initConnCBF(robot_states, x_self, 0);
    last_conn_ac_ = Ac;
    last_conn_bc_ = Bc;

    Vector A0 = Vector::Zero(DIM * this->k_hor_);
    A0.segment(0, DIM) = Ac;                          // 注意 CONTROL_VARS == DIM 时直接放
    Row A_control_pts = -1.0 * A0.transpose() * this->U_basis_;

    return LinearConstraint(A_control_pts, std::numeric_limits<T>::lowest(), Bc + slack_value);
    }

template <typename T, unsigned int DIM>
typename ConnectivityMPCCBFQPOperations<T, DIM>::LinearConstraint
ConnectivityMPCCBFQPOperations<T, DIM>::clfConstraint(const Vector& current_state,
                                                      const Vector& neighbor_state, T slack_value) {
    Vector a = connectivity_cbf_ptr_->getCLFConstraints(current_state, neighbor_state);
    T b = -1.0 * connectivity_cbf_ptr_->getCLFBound(current_state, neighbor_state);
    last_clf_a_ = a;
    last_clf_b_ = b;

    Vector A0 = Vector::Zero(DIM * this->k_hor_);
    A0.segment(0, DIM) = a;
    Row A_control_pts = 1.0 * A0.transpose() * this->U_basis_;
    // ===== 调试输出 =====
    {
        const Eigen::IOFormat oneLine(Eigen::StreamPrecision, Eigen::DontAlignCols,
                                      " ", " ", "", "", "", "");
        std::ostringstream oss;
        oss << A0.transpose().format(oneLine);  // 打印成行向量更清晰
        logger->debug("[clfConstraint] A0 \n{}", oss.str());
    }
    logger->debug("[clfConstraint] b={}", b);
    // ===================
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
    last_pred_safety_Ak_.clear();
    last_pred_safety_bk_.clear();
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
         // === 调试缓存 ===
        last_pred_safety_Ak_.push_back(Ak);
        last_pred_safety_bk_.push_back(bk);

        // === 可选：调试输出 ===
        {
            const Eigen::IOFormat oneLine(Eigen::StreamPrecision, Eigen::DontAlignCols,
                                        " ", " ", "", "", "", "");
            std::ostringstream oss;
            oss << ak.transpose().format(oneLine);
            logger->debug("[predSafetyCBFConstraints] k={} ak=[{}], bk={}", k, oss.str(), bk);
        }
        linear_constraints.push_back(
            LinearConstraint(Ak_control_pts, std::numeric_limits<T>::lowest(), bk));
    }
    return linear_constraints;
}

template <typename T, unsigned int DIM>
std::vector<typename ConnectivityMPCCBFQPOperations<T, DIM>::LinearConstraint>
ConnectivityMPCCBFQPOperations<T, DIM>::predConnectivityConstraints(
    const std::vector<State>& pred_states,
    const std::vector<VectorDIM>& other_positions)
{
    std::vector<LinearConstraint> linear_constraints;
    for (size_t k = 0; k < pred_states.size(); ++k) {
        const State& s = pred_states[k];

        // 1) 组装 robot_states: [ self(pred_k); others ]
        const int N = 1 + static_cast<int>(other_positions.size());
        Eigen::MatrixXd robot_states(N, 6);
        robot_states.setZero();

        // self row: [pos(DIM), vel(DIM)]
        for (unsigned j = 0; j < DIM; ++j) {
            robot_states(0, j)         = static_cast<double>(s.pos_(j));
            robot_states(0, DIM + j)   = static_cast<double>(s.vel_(j));
        }

        // neighbors: positions in DIM dims
        for (int i = 0; i < static_cast<int>(other_positions.size()); ++i) {
            robot_states.block(i + 1, 0, 1, DIM) =
                other_positions[i].template head<DIM>().template cast<double>().transpose();
        }

        // x_self_k for initConnCBF
        Vector x_self_k(2 * DIM);
        x_self_k << s.pos_, s.vel_;

        // 2) 一次性拿 (Ac, Bc)
        auto [Ac, Bc] = connectivity_cbf_ptr_->initConnCBF(robot_states, x_self_k, 0 /* self_idx */);
        last_conn_ac_ = Ac;
        last_conn_bc_ = Bc;

        Vector Ak = Vector::Zero(DIM * this->k_hor_);
        Ak.segment(k * DIM, DIM) = Ac;
        Row Ak_control_pts = -1.0 * Ak.transpose() * this->U_basis_;
        linear_constraints.push_back(
            LinearConstraint(Ak_control_pts, std::numeric_limits<T>::lowest(), Bc));
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
        last_clf_a_ = ak;
        last_clf_b_ = bk;

        Vector Ak = Vector::Zero(DIM * this->k_hor_);
        Ak.segment(k * DIM, DIM) = ak;
        Row Ak_control_pts = 1.0 * Ak.transpose() * this->U_basis_;
        // ===== 调试输出=====  
        {
            const Eigen::IOFormat oneLine(Eigen::StreamPrecision, Eigen::DontAlignCols,
                                          " ", " ", "", "", "", "");
            std::ostringstream oss;
            oss << Ak.transpose().format(oneLine);
            logger->debug("[predCLFConstraints] k={} Ak (horizon)\n{}", k, oss.str());
        }    
        logger->debug("[predCLFConstraints] k={} bk={}", k, bk);
        
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
