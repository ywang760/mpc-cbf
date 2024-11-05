//
// Created by lishuo on 9/22/24.
//

#include <mpc_cbf/controller/BezierIMPCCBF.h>

namespace mpc_cbf {
    template <typename T, unsigned int DIM>
    BezierIMPCCBF<T, DIM>::BezierIMPCCBF(Params &p, std::shared_ptr<DoubleIntegrator> model_ptr,
                                         std::shared_ptr<FovCBF> fov_cbf_ptr, uint64_t bezier_continuity_upto_degree,
                                         std::shared_ptr<const CollisionShape> collision_shape_ptr,
                                         int impc_iter)
    : bezier_continuity_upto_degree_(bezier_continuity_upto_degree),
    collision_shape_ptr_(collision_shape_ptr),
    impc_iter_(impc_iter) {
        // Initiate the qp_generator
        std::unique_ptr<PiecewiseBezierMPCCBFQPOperations> piecewise_mpc_cbf_operations_ptr =
                std::make_unique<PiecewiseBezierMPCCBFQPOperations>(p, model_ptr, fov_cbf_ptr);
        qp_generator_.addPiecewise(std::move(piecewise_mpc_cbf_operations_ptr));

        // load mpc tuning params
        mpc_tuning_ = p.mpc_params.tuning_;

        T h = p.mpc_params.h_;
        T Ts = p.mpc_params.Ts_;
        assert(Ts <= h);
        assert((h - (int)(h/Ts) * Ts) == 0);
        // control sequence U control point coefficient
        int u_interp = int(h / Ts);
        ts_samples_ = Vector::LinSpaced(u_interp, 0, h - Ts); // [I,1]

        // for position pred
        h_ = p.mpc_params.h_;
        k_hor_ = p.mpc_params.k_hor_;
        h_samples_ = Vector::LinSpaced(k_hor_, 0, (k_hor_-1)*h_);
    }

    template <typename T, unsigned int DIM>
    bool BezierIMPCCBF<T, DIM>::optimize(std::vector<SingleParameterPiecewiseCurve> &result_curves,
                                         const State &current_state,
                                         const std::vector<VectorDIM> &other_robot_positions,
                                         const Vector &ref_positions) {
        bool success = true;
        assert(impc_iter_ >= 1);

        for (size_t iter = 0; iter < impc_iter_; ++iter) {
            // reset the problem
            qp_generator_.problem().resetProblem();

            size_t num_pieces = qp_generator_.piecewise_mpc_qp_generator_ptr()->numPieces();
            // add the position error cost
            qp_generator_.piecewise_mpc_qp_generator_ptr()->addPositionErrorPenaltyCost(current_state, ref_positions);
//            qp_generator_.piecewise_mpc_qp_generator_ptr()->addEvalPositionErrorPenaltyCost(ref_positions);
            // minimize the control effort for the curve up to specified degree
            for (size_t d = 1; d <= bezier_continuity_upto_degree_; ++d) {
                qp_generator_.piecewise_mpc_qp_generator_ptr()->addIntegratedSquaredDerivativeCost(d,
                                                                                                   mpc_tuning_.w_u_eff_);
            }

            // add the current state constraints
            VectorDIM current_pos = current_state.pos_;
            VectorDIM current_vel = current_state.vel_;
            qp_generator_.piecewise_mpc_qp_generator_ptr()->addEvalConstraint(0, 0, current_pos);
            qp_generator_.piecewise_mpc_qp_generator_ptr()->addEvalConstraint(0, 1, current_vel);

            // add the continuity between pieces
            for (size_t piece_index = 0; piece_index < num_pieces - 1; ++piece_index) {
                for (size_t d = 0; d < bezier_continuity_upto_degree_; ++d) {
                    qp_generator_.piecewise_mpc_qp_generator_ptr()->addContinuityConstraint(piece_index, d);
                }
            }

            // cbf constraints
            if (iter == 0) {
                for (size_t i = 0; i < other_robot_positions.size(); ++i) {
                    Vector other_xy(2);
                    other_xy << other_robot_positions.at(i)(0), other_robot_positions.at(i)(1);
                    qp_generator_.addSafetyCBFConstraint(current_state, other_xy);
                    qp_generator_.addFovLBConstraint(current_state, other_xy);
                    qp_generator_.addFovRBConstraint(current_state, other_xy);
                }
            } else if (iter > 0 && success) {
                // pred the robot's position in the horizon, use for the CBF constraints
                std::vector<State> pred_states;
                for (size_t k = 0; k < 4; ++k) {
                    State pred_state;
                    pred_state.pos_ = result_curves.back().eval(h_samples_(k), 0);
                    pred_state.vel_ = result_curves.back().eval(h_samples_(k), 1);
                    pred_states.push_back(pred_state);
                }
                for (size_t i = 0; i < other_robot_positions.size(); ++i) {
                    Vector other_xy(2);
                    other_xy << other_robot_positions.at(i)(0), other_robot_positions.at(i)(1);
//                    qp_generator_.addPredSafetyCBFConstraints(pred_states, other_xy);
                    qp_generator_.addSafetyCBFConstraint(current_state, other_xy);
                    qp_generator_.addPredFovLBConstraints(pred_states, other_xy);
                    qp_generator_.addPredFovRBConstraints(pred_states, other_xy);
                }
            }

            // dynamics constraints
            T a_max = 1;
            VectorDIM a_min_vec = {-a_max, -a_max, -0.5*a_max};
            VectorDIM a_max_vec = {a_max, a_max, 0.5*a_max};
            AlignedBox derivative_bbox(a_min_vec, a_max_vec);
            qp_generator_.piecewise_mpc_qp_generator_ptr()->addBoundingBoxConstraintAll(derivative_bbox, 1);
            qp_generator_.piecewise_mpc_qp_generator_ptr()->addBoundingBoxConstraintAll(derivative_bbox, 2);

            // solve QP
            Problem &problem = qp_generator_.problem();
            CPLEXSolver cplex_solver;
            SolveStatus solve_status = cplex_solver.solve(problem);
            if (solve_status == SolveStatus::OPTIMAL) {
                success = true;
            } else {
                success = false;
            }

            // generate bezier
            result_curves.push_back(qp_generator_.piecewise_mpc_qp_generator_ptr()->generateCurveFromSolution());
        }

        return success;
    }

    template <typename T, unsigned int DIM>
    void BezierIMPCCBF<T, DIM>::resetProblem() {
        qp_generator_.problem().resetProblem();
    }

    template <typename T, unsigned int DIM>
    typename BezierIMPCCBF<T, DIM>::Vector
    BezierIMPCCBF<T, DIM>::generatorDerivativeControlInputs(uint64_t derivative_degree) {
        Matrix U_basis = qp_generator_.piecewise_mpc_qp_generator_ptr()->piecewise_operations_ptr()->evalSamplingBasisMatrix(ts_samples_, derivative_degree); // [3I, num_piece*dim*num_control_pts]
        Vector control_pts_variables_value = qp_generator_.piecewise_mpc_qp_generator_ptr()->getVariablesValue(); // [num_piece*dim*num_control_pts, 1]
        return U_basis * control_pts_variables_value;
    }

    template class BezierIMPCCBF<double, 3U>;
} // mpc_cbf