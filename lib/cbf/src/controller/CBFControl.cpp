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
                                      const State &current_state,
                                      const std::vector<VectorDIM> &other_robot_positions,
                                      const VectorDIM& u_min,
                                      const VectorDIM& u_max) {
        // add cost
        qp_generator_.addDesiredControlCost(desired_u);
        // add constraints
        Vector state(2*DIM);
        state << current_state.pos_, current_state.vel_;
        for (size_t i = 0; i < other_robot_positions.size(); ++i) {
            Vector other_xy(2);
            other_xy << other_robot_positions.at(i)(0), other_robot_positions.at(i)(1);
            qp_generator_.addSafetyConstraint(state, other_xy);
            qp_generator_.addLeftBorderConstraint(state, other_xy);
            qp_generator_.addRightBorderConstraint(state, other_xy);
            qp_generator_.addRangeConstraint(state, other_xy);
        }
        qp_generator_.addMinVelConstraints(state);
        qp_generator_.addMaxVelConstraints(state);
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

    template <typename T, unsigned int DIM>
    bool CBFControl<T, DIM>::optimizeWithSlackVariables(VectorDIM& cbf_u,
                                                        const VectorDIM &desired_u,
                                                        const Vector &state,
                                                        const std::vector<VectorDIM> &other_robots_states,
                                                        const std::vector<T> &slacks,
                                                        const VectorDIM& u_min,
                                                        const VectorDIM& u_max) {

        // add cost
        qp_generator_.addDesiredControlCost(desired_u);
        // add constraints
        for (int i = 0; i < other_robots_states.size(); ++i)
        {
            qp_generator_.addSafetyConstraint(state, other_robots_states.at(i));
            qp_generator_.addLeftBorderConstraintWithSlackVar(state, other_robots_states.at(i), slacks.at(i));
            qp_generator_.addRightBorderConstraintWithSlackVar(state, other_robots_states.at(i), slacks.at(i));
            // qp_generator_.addRangeConstraintWithSlackVar(state, other_robots_states.at(i), slacks.at(i));
        }
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