//
// Created on 4/10/25.
//

#include <cbf/controller/ConnectivityControl.h>

namespace cbf {
    template <typename T, unsigned int DIM>
    ConnectivityControl<T, DIM>::ConnectivityControl(std::shared_ptr<ConnectivityCBF> cbf, int num_robots, bool slack_mode, T slack_cost, T slack_decay_rate)
        : qp_generator_(cbf, num_robots, slack_mode), slack_mode_(slack_mode), slack_cost_(slack_cost), slack_decay_rate_(slack_decay_rate), num_robots_(num_robots)
    {
    }

    template <typename T, unsigned int DIM>
    bool ConnectivityControl<T, DIM>::optimize(VectorDIM &cbf_u,
                                               const VectorDIM &desired_u,
                                               std::vector<State> current_states,
                                               size_t self_idx,
                                               const VectorDIM &u_min,
                                               const VectorDIM &u_max)
    {

        State current_state = current_states.at(self_idx);
        Vector state(2 * DIM);
        state << current_state.pos_, current_state.vel_;

        // if slack_mode, compute the slack weights
        std::vector<double> slack_weights(num_robots_);
        if (slack_mode_) {
            // Define slack weights with decay
            T w_init = slack_cost_;
            T decay_factor = slack_decay_rate_;
            for (size_t i = 0; i < num_robots_; ++i) {
                slack_weights.at(i) = w_init * pow(decay_factor, i);
            }
        }

        // add cost
        qp_generator_.addDesiredControlCost(desired_u);
        if (slack_mode_) {
            qp_generator_.addSlackCost(slack_weights);
        }
        
        // add constraints

        // Add safety constraints
        for (size_t i = 0; i < num_robots_ - 1; ++i)
        {
            Vector neighbor_state(2 * DIM);
            neighbor_state << current_states.at(i + (i >= self_idx ? 1 : 0)).pos_,
                              current_states.at(i + (i >= self_idx ? 1 : 0)).vel_;

            if (!slack_mode_) {
                qp_generator_.addSafetyConstraint(state, neighbor_state);
            } else {
                qp_generator_.addSafetyConstraint(state, neighbor_state, true, i);
            }
        }

        // Add velocity and acceleration constraints
        qp_generator_.addMinVelConstraints(state);
        qp_generator_.addMaxVelConstraints(state);
        // qp_generator_.addControlBoundConstraint(u_min, u_max);

        // Add connectivity constraint
        Eigen::MatrixXd robot_states(current_states.size(), 2 * DIM);
        for (size_t i = 0; i < current_states.size(); ++i) {
            robot_states.row(i).head(DIM) = current_states[i].pos_.transpose();
            robot_states.row(i).tail(DIM) = current_states[i].vel_.transpose();
        }
        if (slack_mode_) {
            // Use the last slack variable for connectivity constraint
            qp_generator_.addConnConstraint(state, robot_states, self_idx, true);
        } else {
            qp_generator_.addConnConstraint(state, robot_states, self_idx);
        }

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

    // Explicit template instantiation
    template class ConnectivityControl<double, 3U>;

} // cbf