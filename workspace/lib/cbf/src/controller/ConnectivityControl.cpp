//
// Created on 4/10/25.
//

#include <cbf/controller/ConnectivityControl.h>

namespace cbf {
    auto logger = spdlog::default_logger();
    template <typename T, unsigned int DIM>
    ConnectivityControl<T, DIM>::ConnectivityControl(std::shared_ptr<ConnectivityCBF> cbf, int number_neighbors, bool slack_mode, T slack_cost, T slack_decay_rate)
        : cbf_(cbf), qp_generator_(cbf, number_neighbors, slack_mode), slack_mode_(slack_mode), slack_cost_(slack_cost), slack_decay_rate_(slack_decay_rate)
    {
    }

    template <typename T, unsigned int DIM>
    bool ConnectivityControl<T, DIM>::optimize(VectorDIM &cbf_u,
                                               const VectorDIM &desired_u,
                                               std::vector<State> current_states,
                                               size_t ego_robot_idx,
                                               const VectorDIM &u_min,
                                               const VectorDIM &u_max)
    {

        State current_state = current_states.at(ego_robot_idx);

        // For backward compatibility
        std::vector<VectorDIM> other_robot_positions;
        for (size_t i = 0; i < current_states.size(); ++i)
        {
            if (i != ego_robot_idx)
            {

                VectorDIM &pos = current_states.at(i).pos_;
                pos(2) = 0; // Set z-coordinate to zero for 2D control
                other_robot_positions.push_back(pos);
            }
        }

        // if slack_mode, compute the slack weights
        int num_neighbors = current_states.size() - 1; // Exclude the ego robot itself
        std::vector<double> slack_weights(num_neighbors);
        if (slack_mode_) {
            // Define slack weights with decay
            T w_init = slack_cost_;
            T decay_factor = slack_decay_rate_;
            for (size_t i = 0; i < num_neighbors; ++i) {
                slack_weights.at(i) = w_init * pow(decay_factor, i);
            }
        }

        // add cost
        qp_generator_.addDesiredControlCost(desired_u);
        if (slack_mode_) {
            qp_generator_.addSlackCost(slack_weights);
        }
        
        // add constraints
        Vector state(2*DIM);
        state << current_state.pos_, current_state.vel_;

        for (size_t i = 0; i < num_neighbors; ++i)
        {
            Vector neighbor_state(6);
            neighbor_state << current_states.at(i + (i >= ego_robot_idx ? 1 : 0)).pos_, 
                              current_states.at(i + (i >= ego_robot_idx ? 1 : 0)).vel_;

            if (!slack_mode_) {
                qp_generator_.addSafetyConstraint(state, neighbor_state);
            } else {
                qp_generator_.addSafetyConstraint(state, neighbor_state, true, i);
            }
        }

        // Add velocity constraints
        qp_generator_.addMinVelConstraints(state);
        qp_generator_.addMaxVelConstraints(state);
        qp_generator_.addControlBoundConstraint(u_min, u_max);

        Eigen::MatrixXd positions(current_states.size(), 3);
        for (size_t i = 0; i < current_states.size(); ++i) {
            positions.row(i) = current_states.at(i).pos_;
        }
        double Rs_value = cbf_->getDmax();
        double sigma_val = std::pow(Rs_value, 4) / std::log(2.0);
        auto [lambda2, lambda2_vec] = getLambda2FromL(positions, Rs_value, sigma_val);
        logger->info("Current λ₂ value: {}", lambda2);
        // Add connectivity constraint
        if (lambda2 > 0.1) {
            qp_generator_.addConnConstraint(state, other_robot_positions, false, 0);
        } else {
            int local_slack_idx = 0;
            for (size_t i = 0; i < current_states.size(); ++i)
            {
                if (i == ego_robot_idx) continue;
                Vector neighbor_state(6);
                neighbor_state.segment(0, 3) = current_states.at(i).pos_;
                neighbor_state.segment(3, 3) = current_states.at(i).vel_;
                if (!slack_mode_) {
                    qp_generator_.addCLFConstraint(state, neighbor_state, false, 0);
                } else {
                    qp_generator_.addCLFConstraint(state, neighbor_state, true, local_slack_idx);
                    ++local_slack_idx;
                }
            }
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