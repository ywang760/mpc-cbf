//
// Created by yutong on 8/4/25.
//

#include <mpc_cbf/controller/ConnectivityIMPCCBF.h>

namespace mpc_cbf {
template <typename T, unsigned int DIM>
ConnectivityIMPCCBF<T, DIM>::ConnectivityIMPCCBF(
    Params& p, std::shared_ptr<DoubleIntegrator> model_ptr,
    std::shared_ptr<ConnectivityCBF> connectivity_cbf_ptr, uint64_t bezier_continuity_upto_degree,
    std::shared_ptr<const CollisionShape> collision_shape_ptr, int num_neighbors)
    : bezier_continuity_upto_degree_(bezier_continuity_upto_degree),
      collision_shape_ptr_(collision_shape_ptr) {
    // impc params
    impc_iter_ = p.impc_params.impc_iter_;
    cbf_horizon_ = p.impc_params.cbf_horizon_;
    slack_cost_ = p.impc_params.slack_cost_;
    slack_decay_rate_ = p.impc_params.slack_decay_rate_;
    slack_mode_ = p.impc_params.slack_mode_;
    // Initiate the qp_generator
    std::unique_ptr<ConnectivityMPCCBFQPOperations> connectivity_mpc_cbf_operations_ptr =
        std::make_unique<ConnectivityMPCCBFQPOperations>(p.mpc_cbf_params, model_ptr,
                                                         connectivity_cbf_ptr);
    qp_generator_.addPiecewise(std::move(connectivity_mpc_cbf_operations_ptr), num_neighbors,
                               slack_mode_);

    // load mpc tuning params
    mpc_tuning_ = p.mpc_cbf_params.mpc_params.tuning_;
    v_min_ = p.mpc_cbf_params.mpc_params.limits_.v_min_;
    v_max_ = p.mpc_cbf_params.mpc_params.limits_.v_max_;
    a_min_ = p.mpc_cbf_params.mpc_params.limits_.a_min_;
    a_max_ = p.mpc_cbf_params.mpc_params.limits_.a_max_;

    const T h = p.mpc_cbf_params.mpc_params.h_;
    const T Ts = p.mpc_cbf_params.mpc_params.Ts_;
    assert(Ts <= h);
    assert(std::abs(h - static_cast<int>(h / Ts) * Ts) <
           1e-10); // Use epsilon for floating point comparison

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

    // Extract other robot states (both position and velocity) from
    // current_states (excluding self)
    std::vector<VectorDIM> other_robot_positions;
    std::vector<Vector> other_robot_states; // Full state vectors for CBF
    for (size_t i = 0; i < current_states.size(); ++i) {
        if (i != self_idx) {
            other_robot_positions.push_back(current_states[i].pos_);
            // Create full state vector (position + velocity)
            Vector neighbor_state_vec(2 * DIM);
            neighbor_state_vec << current_states[i].pos_, current_states[i].vel_;
            other_robot_states.push_back(neighbor_state_vec);
        }
    }

    Vector state(2 * DIM);
    state << current_pos, current_vel;

    // if slack_mode, compute the slack weights
    const size_t num_neighbors = other_robot_states.size();
    std::vector<double> slack_weights(num_neighbors);
    if (slack_mode_) {
        std::vector<size_t> idx(num_neighbors);
        std::iota(idx.begin(), idx.end(), 0);

        // Pre-compute distances to avoid repeated calculations
        std::vector<T> distances(num_neighbors);
        for (size_t i = 0; i < num_neighbors; ++i) {
            // Use spatial dimensions only (ignore orientation)
            constexpr size_t spatial_dims = 2; // FIXME: adapt to 3D
            distances[i] =
                (other_robot_positions[i].head(spatial_dims) - current_pos.head(spatial_dims))
                    .norm();
        }

        // Sort by pre-computed distances (closer robots get priority)
        std::sort(idx.begin(), idx.end(),
                  [&distances](size_t a, size_t b) { return distances[a] < distances[b]; });

        // Define slack weights with decay (closer robots get higher slack cost)
        const T w_init = slack_cost_;
        const T decay_factor = slack_decay_rate_;
        for (size_t i = 0; i < num_neighbors; ++i) {
            slack_weights[idx[i]] = w_init * std::pow(decay_factor, static_cast<T>(i));
        }
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

        // add the slack cost, to give priority to neighbors
        if (slack_mode_) {
            qp_generator_.addSlackCost(slack_weights);
        }

        // add the current state constraints
        qp_generator_.piecewise_mpc_qp_generator_ptr()->addEvalConstraint(0, 0, current_pos);
        qp_generator_.piecewise_mpc_qp_generator_ptr()->addEvalConstraint(0, 1, current_vel);

        // add the continuity between pieces
        for (size_t piece_index = 0; piece_index < num_pieces - 1; ++piece_index) {
            for (size_t d = 0; d < bezier_continuity_upto_degree_; ++d) {
                qp_generator_.piecewise_mpc_qp_generator_ptr()->addContinuityConstraint(piece_index,
                                                                                        d);
            }
        }

        // add the collision avoidance constraints
        AlignedBox robot_bbox_at_zero = collision_shape_ptr_->boundingBox(VectorDIM::Zero());
        for (size_t i = 0; i < num_neighbors; ++i) {
            const VectorDIM& other_robot_pos = other_robot_positions[i];

            // TODO: this is model specific, here the last dimension is
            // orientation data
            VectorDIM current_xy = current_pos;
            current_xy(DIM - 1) = 0;
            VectorDIM other_robot_xy = other_robot_pos;
            other_robot_xy(DIM - 1) = 0;

            Hyperplane hyperplane =
                separating_hyperplanes::voronoi<T, DIM>(current_xy, other_robot_xy);
            // shifted hyperplane
            Hyperplane shifted_hyperplane =
                math::shiftHyperplane<T, DIM>(hyperplane, robot_bbox_at_zero);
            qp_generator_.piecewise_mpc_qp_generator_ptr()->addHyperplaneConstraintForPiece(
                0,
                shifted_hyperplane); // add the collision avoidance
                                     // constraint on the first segment
        }

        // connectivity cbf constraints
        if (iter == 0) {
            for (size_t i = 0; i < num_neighbors; ++i) {
                T slack_value = 0; // TODO: pass in param
                if (!slack_mode_) {
                    qp_generator_.addSafetyCBFConstraint(state, other_robot_states[i], slack_value);
                    // FIXME: Implement connectivity-specific constraint methods
                    // qp_generator_.addConnectivityConstraint(state,
                    // other_robot_states[i], slack_value);
                } else {
                    qp_generator_.addSafetyCBFConstraintWithSlackVariables(
                        state, other_robot_states[i], i);
                    // FIXME: Implement connectivity-specific constraint methods
                    // with slack
                    // qp_generator_.addConnectivityConstraintWithSlackVariables(state,
                    // other_robot_states[i], i);
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
            const std::vector<T> slack_values(cbf_horizon_,
                                              0); // TODO: pass in param

            for (size_t i = 0; i < num_neighbors; ++i) {
                if (!slack_mode_) {
                    qp_generator_.addPredSafetyCBFConstraints(pred_states, other_robot_states[i],
                                                              slack_values);
                    // FIXME: Implement predicted connectivity constraint methods
                    // qp_generator_.addPredConnectivityConstraints(pred_states,
                    // other_robot_states[i], slack_values);
                } else {
                    qp_generator_.addPredSafetyCBFConstraintsWithSlackVariables(
                        pred_states, other_robot_states[i], i);
                    // FIXME: Implement predicted connectivity constraint methods
                    // with slack
                    // qp_generator_.addPredConnectivityConstraintsWithSlackVariables(pred_states,
                    // other_robot_states[i], i);
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
            break;
        }
    }

    return success;
}

template <typename T, unsigned int DIM>
void ConnectivityIMPCCBF<T, DIM>::resetProblem() {
    qp_generator_.problem().resetProblem();
}

template class ConnectivityIMPCCBF<double, 3U>;
} // namespace mpc_cbf