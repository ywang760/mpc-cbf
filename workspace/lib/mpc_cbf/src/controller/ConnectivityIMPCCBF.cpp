//
// Created by yutong on 8/4/25.
//

#include <mpc_cbf/controller/ConnectivityIMPCCBF.h>
#include <spdlog/spdlog.h>
namespace mpc_cbf {
template <typename T, unsigned int DIM>
ConnectivityIMPCCBF<T, DIM>::ConnectivityIMPCCBF(
    Params& p, std::shared_ptr<DoubleIntegrator> model_ptr,
    std::shared_ptr<ConnectivityCBF> connectivity_cbf_ptr,
    std::shared_ptr<const CollisionShape> collision_shape_ptr, int num_neighbors)
    : bezier_continuity_upto_degree_(
          p.mpc_cbf_params.piecewise_bezier_params.bezier_continuity_upto_degree_),
      collision_shape_ptr_(collision_shape_ptr),
      qp_generator_(std::make_unique<ConnectivityMPCCBFQPOperations>(p.mpc_cbf_params, model_ptr,
                                                                     connectivity_cbf_ptr),
                    num_neighbors, p.impc_params.slack_config_) {
    // impc params
    impc_iter_ = p.impc_params.impc_iter_;
    cbf_horizon_ = p.impc_params.cbf_horizon_;
    slack_config_ = p.impc_params.slack_config_;

    // load mpc tuning params
    mpc_tuning_ = p.mpc_cbf_params.mpc_params.tuning_;
    v_min_ = p.mpc_cbf_params.mpc_params.limits_.v_min_;
    v_max_ = p.mpc_cbf_params.mpc_params.limits_.v_max_;
    a_min_ = p.mpc_cbf_params.mpc_params.limits_.a_min_;
    a_max_ = p.mpc_cbf_params.mpc_params.limits_.a_max_;

    const T h = p.mpc_cbf_params.mpc_params.h_;
    const T Ts = p.mpc_cbf_params.mpc_params.Ts_;

    // control sequence U control point coefficient
    const int u_interp = static_cast<int>(h / Ts);
    ts_samples_ = Vector::LinSpaced(u_interp, 0, h - Ts);

    // for position prediction
    h_ = h;
    k_hor_ = p.mpc_cbf_params.mpc_params.k_hor_;
    h_samples_ = Vector::LinSpaced(k_hor_, 0, static_cast<T>(k_hor_ - 1) * h_);
}

template <typename T, unsigned int DIM>
bool ConnectivityIMPCCBF<T, DIM>::optimize(
    std::vector<SingleParameterPiecewiseCurve>& result_curves,
    const std::vector<State>& current_states, size_t self_idx, const Vector& ref_positions) {
    bool success = true;
    assert(impc_iter_ >= 1);

    const State& current_state = current_states.at(self_idx);
    const VectorDIM& current_pos = current_state.pos_;
    const VectorDIM& current_vel = current_state.vel_;

    // Create robot_states matrix with all robot states (including self)
    Eigen::MatrixXd robot_states(current_states.size(), 2 * DIM);
    // A) 全状态（给 Safety/CLF 用，排除自己）
    std::vector<Vector> other_robot_states;
    other_robot_states.reserve(current_states.size() - 1);
    // B) 位置（给 Connectivity 用，排除自己）
    std::vector<VectorDIM> other_positions;
    other_positions.reserve(current_states.size() - 1);
    for (size_t i = 0; i < current_states.size(); ++i) {
        Vector robot_state(2 * DIM);
        robot_state << current_states[i].pos_, current_states[i].vel_;
        robot_states.row(i) = robot_state.transpose();

        if (i != self_idx) {
            other_robot_states.push_back(robot_state);                    // 全状态
            other_positions.push_back(robot_state.template head<DIM>());  // 仅位置（定长）
        }
    }
    Vector state(2 * DIM);
    state << current_pos, current_vel;

    // Compute slack weights for different constraint types based on distance
    const size_t num_neighbors = other_robot_states.size();
    std::vector<size_t> idx(num_neighbors);
    std::iota(idx.begin(), idx.end(), 0);

        // Pre-compute distances to avoid repeated calculations
        std::vector<T> distances(num_neighbors);
        for (size_t i = 0; i < num_neighbors; ++i) {
            // Use spatial dimensions only (ignore orientation) from other_robot_states
            constexpr size_t spatial_dims = 3; // FIXME: adapt to 3D
            VectorDIM other_robot_pos = other_robot_states[i].head(DIM);
            distances[i] =
                (other_robot_pos.head(spatial_dims) - current_pos.head(spatial_dims)).norm();
        }

    // Sort by pre-computed distances (closer robots get priority)
    std::sort(idx.begin(), idx.end(),
              [&distances](size_t a, size_t b) { return distances[a] < distances[b]; });

    // Prepare separate slack weights for each constraint type
    std::vector<double> safety_weights(num_neighbors);
    std::vector<double> clf_weights(num_neighbors);
    std::vector<double> connectivity_weights = {slack_config_.connectivity_slack_cost};

    // Apply distance-based weights with decay for safety and CLF constraints
    for (size_t i = 0; i < num_neighbors; ++i) {
        safety_weights[idx[i]] = slack_config_.safety_slack_cost *
                                 std::pow(slack_config_.slack_decay_rate, static_cast<T>(i));
        clf_weights[idx[i]] = slack_config_.clf_slack_cost *
                              std::pow(slack_config_.slack_decay_rate, static_cast<T>(i));
    }

    for (size_t iter = 0; iter < impc_iter_; ++iter) {
        // reset the problem
        resetProblem();

        size_t num_pieces = qp_generator_.piecewise_mpc_qp_generator_ptr()->numPieces();
        // add the position error cost
        qp_generator_.piecewise_mpc_qp_generator_ptr()->addPositionErrorPenaltyCost(current_state,
                                                                                    ref_positions);

        // minimize the control effort for the curve up to specified degree
        for (size_t d = 1; d <= bezier_continuity_upto_degree_; ++d) {
            qp_generator_.piecewise_mpc_qp_generator_ptr()->addIntegratedSquaredDerivativeCost(
                d, mpc_tuning_.w_u_eff_);
        }

        // Add separate slack costs based on configuration
        if (slack_config_.safety_slack && !qp_generator_.safety_slack_variables_.empty()) {
            qp_generator_.addSlackCost(safety_weights, qp_generator_.safety_slack_variables_);
        }

        if (slack_config_.clf_slack && !qp_generator_.clf_slack_variables_.empty()) {
            qp_generator_.addSlackCost(clf_weights, qp_generator_.clf_slack_variables_);
        }

        if (slack_config_.connectivity_slack &&
            !qp_generator_.connectivity_slack_variables_.empty()) {
            qp_generator_.addSlackCost(connectivity_weights,
                                       qp_generator_.connectivity_slack_variables_);
        }

        // add the current state constraints
        qp_generator_.piecewise_mpc_qp_generator_ptr()->addEvalConstraint(0, 0, current_pos);
        qp_generator_.piecewise_mpc_qp_generator_ptr()->addEvalConstraint(0, 1, current_vel);

        // add the continuity between pieces
        for (size_t piece_index = 0; piece_index < num_pieces - 1; ++piece_index) {
            for (size_t d = 0; d <= bezier_continuity_upto_degree_; ++d) {
                qp_generator_.piecewise_mpc_qp_generator_ptr()->addContinuityConstraint(piece_index,
                                                                                        d);
            }
        }

        // connectivity cbf constraints
        if (iter == 0) {
            // TODO: currently all slack_value are set as 0a
            // Set to different values for priority scheduling
            T slack_value = 0;
            for (size_t i = 0; i < num_neighbors; ++i) {
                qp_generator_.addSafetyCBFConstraint(state, other_robot_states[i], i, slack_value);
            }

            //Evaluate lambda2 and conditionally add connectivity or CLF constraints
            const auto robot_positions = robot_states.leftCols(DIM); // Extract only position columns (x, y)
            const double Rs = qp_generator_.connectivityCBF()->getDmax();            
            const double sigma = std::pow(Rs, 4) / std::log(2.0);                    
            auto [lambda2, eigenvec] = mpc_cbf::getLambda2FromL(robot_positions, Rs, sigma);
            last_lambda2_ = static_cast<T>(lambda2); // 记录最近一次 λ2
            // TODO: this 0.1 threshold is arbitrary - should be configurable
            if (lambda2 > 0.1) {
            //     // Use single connectivity constraint when graph is well-connected
                qp_generator_.addConnectivityConstraint(state, other_positions, 0);
                current_mode_ = ConstraintMode::Connectivity;
            } else {
                spdlog::info("enter CLF");
                current_mode_ = ConstraintMode::CLF;
                // Use pairwise CLF constraints when graph connectivity is poor
                for (size_t i = 0; i < num_neighbors; ++i) {
                    qp_generator_.addCLFConstraint(state, other_robot_states[i], i, slack_value);
                }
            }
        } else if (iter > 0 && success) {
            // pred the robot's position in the horizon, use for the CBF
            // constraints
            std::vector<State> pred_states;
            pred_states.reserve(cbf_horizon_);
            for (size_t k = 0; k < cbf_horizon_; ++k) {
                State pred_state;
                pred_state.pos_ = result_curves.back().eval(h_samples_(k), 0);
                pred_state.vel_ = result_curves.back().eval(h_samples_(k), 1);
                pred_states.push_back(pred_state);
            }
            // Create slack values once per neighbor (not growing with horizon)
            // TODO: currently all slack_values are set as 0
            // Set to different values for priority scheduling
            const std::vector<T> slack_values(cbf_horizon_, 0);
            for (size_t i = 0; i < num_neighbors; ++i) {
                qp_generator_.addPredSafetyCBFConstraints(pred_states, other_robot_states[i], i);
            }

            // Evaluate lambda2 and conditionally add predicted connectivity or CLF constraints
            const auto robot_positions = robot_states.leftCols(3); // Extract only position columns (x, y)
            const double Rs = qp_generator_.connectivityCBF()->getDmax();            
            const double sigma = std::pow(Rs, 4) / std::log(2.0);                    
            auto [lambda2, eigenvec] = mpc_cbf::getLambda2FromL(robot_positions, Rs, sigma);
            last_lambda2_ = static_cast<T>(lambda2); // 记录最近一次 λ2
            // TODO: this 0.1 threshold is arbitrary - should be configurable
            if (lambda2 > 0.1) {
                current_mode_ = ConstraintMode::Connectivity;
                // Use single connectivity constraint when graph is well-connected
                qp_generator_.addPredConnectivityConstraints(pred_states, other_positions, 0);
            } else {
                spdlog::info("enter CLF");
                current_mode_ = ConstraintMode::CLF;
                // Use pairwise CLF constraints when graph connectivity is poor
                for (size_t i = 0; i < num_neighbors; ++i) {
                    qp_generator_.addPredCLFConstraints(pred_states, other_robot_states[i], i);
                }
            }
        }

        // dynamics constraints
        qp_generator_.piecewise_mpc_qp_generator_ptr()->addEvalBoundConstraints(2, a_min_, a_max_);
        qp_generator_.piecewise_mpc_qp_generator_ptr()->addEvalBoundConstraints(1, v_min_, v_max_);

        // solve QP
        Problem& problem = qp_generator_.problem();
        CPLEXSolver cplex_solver;
        SolveStatus solve_status = cplex_solver.solve(problem);
        if (solve_status == SolveStatus::OPTIMAL) {
            success = true;
            // generate bezier
            result_curves.push_back(
                qp_generator_.piecewise_mpc_qp_generator_ptr()->generateCurveFromSolution());
        } else {
            success = false;
            // === 统计失败归因 ===
            ++fail_cnt_total_;
            switch (current_mode_) {
                case ConstraintMode::Connectivity: {
                    ++fail_cnt_connectivity_;
                    spdlog::warn("[QP FAIL] mode=Connectivity, lambda2={:.6f}",
                        static_cast<double>(last_lambda2_));
           // 仅收集，不打印
                    if (samples_conn_.size() < kMaxSamples) {
                        samples_conn_.emplace_back(
                            qp_generator_.debugLastConn_ac(),
                            qp_generator_.debugLastConn_bc()
                        );
                    }
                        break;
                }
                case ConstraintMode::CLF: {
                    ++fail_cnt_clf_;
                    spdlog::warn("[QP FAIL] mode=CLF, lambda2={:.6f}",
                        static_cast<double>(last_lambda2_));
                    // 仅收集，不打印
                    if (samples_clf_.size() < kMaxSamples) {
                        samples_clf_.emplace_back(
                            qp_generator_.debugLastCLF_a(),
                            qp_generator_.debugLastCLF_b()
                        );
                    }
                    break;
                }
                case ConstraintMode::Safety: {
                    ++fail_cnt_safety_;
                    spdlog::warn("[QP FAIL] mode=Safety, lambda2={:.6f}", static_cast<double>(last_lambda2_));

                    if (samples_safety_.size() < kMaxSamples) {
                        samples_safety_.emplace_back(
                            qp_generator_.debugLastSafety_a(),   // 就是我们刚加的接口
                            qp_generator_.debugLastSafety_b()
                        );
                    }
                    break;
                }
                default:
                    spdlog::warn("[QP FAIL] mode=None (not set), lambda2={:.6f}",
                                 static_cast<double>(last_lambda2_));
                    break;
            }
            
            break;
        }
    }

    return success;
}

template <typename T, unsigned int DIM>
void ConnectivityIMPCCBF<T, DIM>::resetProblem() {
    qp_generator_.problem().resetProblem();
    current_mode_ = ConstraintMode::Safety; 
}

template class ConnectivityIMPCCBF<double, 3U>;
} // namespace mpc_cbf