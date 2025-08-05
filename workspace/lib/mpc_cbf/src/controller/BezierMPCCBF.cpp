//
// Created by lishuo on 9/21/24.
//

#include <mpc_cbf/controller/BezierMPCCBF.h>

namespace mpc_cbf {
template <typename T, unsigned int DIM>
BezierMPCCBF<T, DIM>::BezierMPCCBF(Params& p, std::shared_ptr<DoubleIntegrator> model_ptr,
                                   std::shared_ptr<FovCBF> fov_cbf_ptr,
                                   uint64_t bezier_continuity_upto_degree,
                                   std::shared_ptr<const CollisionShape> collision_shape_ptr)
    : bezier_continuity_upto_degree_(bezier_continuity_upto_degree),
      collision_shape_ptr_(collision_shape_ptr) {
    // Initiate the qp_generator
    std::unique_ptr<PiecewiseBezierMPCCBFQPOperations> piecewise_mpc_cbf_operations_ptr =
        std::make_unique<PiecewiseBezierMPCCBFQPOperations>(p, model_ptr, fov_cbf_ptr);
    qp_generator_.addPiecewise(std::move(piecewise_mpc_cbf_operations_ptr));

    // load mpc tuning params
    mpc_tuning_ = p.mpc_params.tuning_;

    T h = p.mpc_params.h_;
    T Ts = p.mpc_params.Ts_;
    assert(Ts <= h);
    assert((h - (int) (h / Ts) * Ts) == 0);
    // control sequence U control point coefficient
    int u_interp = int(h / Ts);
    ts_samples_ = Vector::LinSpaced(u_interp, 0, h - Ts); // [I,1]
}

template <typename T, unsigned int DIM>
bool BezierMPCCBF<T, DIM>::optimize(SingleParameterPiecewiseCurve& result_curve,
                                    const State& current_state,
                                    const std::vector<VectorDIM>& other_robot_positions,
                                    const Vector& ref_positions) {
    size_t num_pieces = qp_generator_.piecewise_mpc_qp_generator_ptr()->numPieces();
    // add the position error cost
    qp_generator_.piecewise_mpc_qp_generator_ptr()->addPositionErrorPenaltyCost(current_state,
                                                                                ref_positions);
    // minimize the control effort for the curve up to specified degree
    for (size_t d = 1; d <= bezier_continuity_upto_degree_; ++d) {
        qp_generator_.piecewise_mpc_qp_generator_ptr()->addIntegratedSquaredDerivativeCost(
            d, mpc_tuning_.w_u_eff_);
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
    for (size_t i = 0; i < other_robot_positions.size(); ++i) {
        Vector other_xy(2);
        other_xy << other_robot_positions.at(i)(0), other_robot_positions.at(i)(1);
        qp_generator_.addSafetyCBFConstraint(current_state, other_xy);
        qp_generator_.addFovLBConstraint(current_state, other_xy);
        qp_generator_.addFovRBConstraint(current_state, other_xy);
    }

    // solve QP
    Problem& problem = qp_generator_.problem();
    CPLEXSolver cplex_solver;
    SolveStatus solve_status = cplex_solver.solve(problem);
    bool success;
    if (solve_status == SolveStatus::OPTIMAL) {
        success = true;
    } else {
        success = false;
    }

    // generate bezier
    result_curve = qp_generator_.piecewise_mpc_qp_generator_ptr()->generateCurveFromSolution();

    return success;
}

template <typename T, unsigned int DIM>
typename BezierMPCCBF<T, DIM>::Vector
BezierMPCCBF<T, DIM>::generatorDerivativeControlInputs(uint64_t derivative_degree) {
    Matrix U_basis = qp_generator_.piecewise_mpc_qp_generator_ptr()
                         ->piecewise_operations_ptr()
                         ->evalSamplingBasisMatrix(
                             ts_samples_, derivative_degree); // [3I, num_piece*dim*num_control_pts]
    Vector control_pts_variables_value =
        qp_generator_.piecewise_mpc_qp_generator_ptr()
            ->getVariablesValue(); // [num_piece*dim*num_control_pts, 1]
    return U_basis * control_pts_variables_value;
}

template class BezierMPCCBF<double, 3U>;
} // namespace mpc_cbf