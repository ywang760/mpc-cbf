//
// Created by lishuo on 8/22/24.
//

#include <mpc/controller/BezierMPC.h>

namespace mpc {
    template <typename T, unsigned int DIM>
    BezierMPC<T, DIM>::BezierMPC(Params &p, std::shared_ptr<DoubleIntegrator> model_ptr, uint64_t bezier_continuity_upto_degree)
    : model_ptr_(model_ptr),
    bezier_continuity_upto_degree_(bezier_continuity_upto_degree) {
        // Initiate the qp_generator
        std::unique_ptr<PiecewiseBezierMPCQPOperation> piecewise_bezier_qp_operation =
                std::make_unique<PiecewiseBezierMPCQPOperation>(p, model_ptr);
        qp_generator_.addPiecewise(std::move(piecewise_bezier_qp_operation));

        // load mpc tuning params
        mpc_tuning_ = p.mpc_params.tuning_;

        // mpc execution
        T h = p.mpc_params.h_;
        T Ts = p.mpc_params.Ts_;
        assert(Ts < h);
        assert((h - (int)(h/Ts) * Ts) == 0);

        // model predict
        int u_interp = int(h / Ts);
        A0_ = model_ptr_->get_A0(u_interp); // A0.pos_: [3I, 6], A0.vel_: [3I, 6],
        Lambda_ = model_ptr_->get_lambda(u_interp); // Lambda_.pos_: [3I, 3I], Lambda_.vel_: [3I, 3I]
        // control sequence U control point coefficient
        ts_samples_ = Vector::LinSpaced(u_interp, 0, h - Ts); // [I,1]
    }

    template <typename T, unsigned int DIM>
    bool BezierMPC<T, DIM>::optimize(SingleParameterPiecewiseCurve &result_curve,
                                     const State &current_state, const std::vector<VectorDIM>& other_robot_positions,
                                     const VectorDIM &target) {

        size_t num_pieces = qp_generator_.numPieces();

        // add the position error cost
        qp_generator_.addPositionErrorPenaltyCost(current_state, target);
        // minimize the control effort for the curve up to specified degree
        for (size_t d = 1; d <= bezier_continuity_upto_degree_; ++d) {
            qp_generator_.addIntegratedSquaredDerivativeCost(d, mpc_tuning_.w_u_eff_); // TODO pass in this weight from params
        }
        std::cout << "error after cost\n";

        // add the current state constraints
        VectorDIM current_pos = current_state.pos_;
        VectorDIM current_vel = current_state.vel_;
        qp_generator_.addEvalConstraint(0, 0, current_pos);
        qp_generator_.addEvalConstraint(0, 1, current_vel);

        // add the continuity between pieces
        for (size_t piece_index = 0; piece_index < num_pieces-1; ++piece_index) {
            for (size_t d = 0; d < bezier_continuity_upto_degree_; ++d) {
                qp_generator_.addContinuityConstraint(piece_index, d);
            }
        }

        // add the collision avoidance constraints
        for (size_t i = 0; i < other_robot_positions.size(); ++i) {
            const VectorDIM& other_robot_pos = other_robot_positions.at(i);
            Hyperplane hyperplane = separating_hyperplanes::voronoi<T, DIM>(current_pos, other_robot_pos);
            qp_generator_.addHyperplaneConstraintForPiece(0, hyperplane); // add the collision avoidance constraint on the first segment
        }

        std::cout << "error after constraints\n";

        // solve QP
        Problem &problem = qp_generator_.problem();
        CPLEXSolver cplex_solver;
        SolveStatus solve_status = cplex_solver.solve(problem);
        bool success;
        if (solve_status == SolveStatus::OPTIMAL) {
            success = true;
        } else {
            success = false;
        }

        // generate bezier
        result_curve = qp_generator_.generateCurveFromSolution();

        return success;
    }

    template <typename T, unsigned int DIM>
    typename BezierMPC<T, DIM>::Vector
    BezierMPC<T, DIM>::generatorDerivativeControlInputs(uint64_t derivative_degree) {
        Matrix U_basis = qp_generator_.piecewise_operations_ptr()->evalSamplingBasisMatrix(ts_samples_, derivative_degree); // [3I, num_piece*dim*num_control_pts]
        Vector control_pts_variables_value = qp_generator_.getVariablesValue(); // [num_piece*dim*num_control_pts, 1]
        return U_basis * control_pts_variables_value;
    }

    template class BezierMPC<double, 3U>;
    template class BezierMPC<float, 3U>;
    template class BezierMPC<double, 2U>;
    template class BezierMPC<float, 2U>;
} // mpc