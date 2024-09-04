//
// Created by lishuo on 9/1/24.
//

#include <cbf/controller/CBFControl.h>

namespace cbf {
    template <typename T, unsigned int DIM>
    CBFControl<T, DIM>::CBFControl(std::shared_ptr<FovCBF> cbf) {
        std::unique_ptr<CBFQPOperations> cbf_operations = std::make_unique<CBFQPOperations>(cbf);
        qp_generator_.addCBFOperations(std::move(cbf_operations));
    }

    template <typename T, unsigned int DIM>
    bool CBFControl<T, DIM>::optimize(VectorDIM& cbf_u,
                                      const VectorDIM &desired_u,
                                      const Vector &state,
                                      const Vector &target_state,
                                      const VectorDIM& u_min,
                                      const VectorDIM& u_max) {
        // add cost
        qp_generator_.addDesiredControlCost(desired_u);
        // add constraints
        qp_generator_.addSafetyConstraint(state, target_state);
        qp_generator_.addLeftBorderConstraint(state, target_state);
        qp_generator_.addRightBorderConstraint(state, target_state);
//        qp_generator_.addRangeConstraint(state, target_state);
        qp_generator_.addControlBoundConstraint(u_min, u_max);

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

        cbf_u = qp_generator_.generatorCBFControlInput();
        return success;
    }

    template class CBFControl<double, 3U>;
//    template class CBFControl<float, 3U>;
} // cbf