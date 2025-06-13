//
// Created on 4/10/25.
//

#include <cbf/controller/ConnectivityControl.h>

namespace cbf {
    template <typename T, unsigned int DIM>
    ConnectivityControl<T, DIM>::ConnectivityControl(std::shared_ptr<ConnectivityCBF> cbf, int number_neighbors, bool slack_mode, T slack_cost, T slack_decay_rate)
        : qp_generator_(cbf, number_neighbors, slack_mode), slack_mode_(slack_mode), slack_cost_(slack_cost), slack_decay_rate_(slack_decay_rate)
    {
    }

    template <typename T, unsigned int DIM>
    bool ConnectivityControl<T, DIM>::optimize(VectorDIM& cbf_u,
                                      const VectorDIM &desired_u,
                                      const State &current_state,
                                      const std::vector<VectorDIM> &other_robot_positions,
                                      const std::vector<Matrix> &other_robot_covs,
                                      const VectorDIM& u_min,
                                      const VectorDIM& u_max) {

        VectorDIM current_pos = current_state.pos_;
        VectorDIM current_vel = current_state.vel_;
        
        // if slack_mode, compute the slack weights
        int num_neighbors = other_robot_positions.size();
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
        
        for (size_t i = 0; i < other_robot_positions.size(); ++i) {
            Vector other_xy(2);
            other_xy << other_robot_positions.at(i)(0), other_robot_positions.at(i)(1);
            
            if (!slack_mode_) {
                qp_generator_.addSafetyConstraint(state, other_xy);
            } else {
                qp_generator_.addSafetyConstraint(state, other_xy, true, i);
            }
        }
        
        // Add velocity constraints
        qp_generator_.addMinVelConstraints(state);
        qp_generator_.addMaxVelConstraints(state);
        qp_generator_.addControlBoundConstraint(u_min, u_max);

        // Convert from VectorDIM (Vector3d) to Vector (VectorXd) // TODO: is this necessary?
        std::vector<Vector> other_positions_dyn;
        for (const auto& p : other_robot_positions) {
            Vector v(p.size());
            v = p;
            other_positions_dyn.push_back(v);
        }

        // Add connectivity constraints
        if (slack_mode_) {
            qp_generator_.addConnectivityConstraint(state, other_positions_dyn, true, num_neighbors);
        } else {
            qp_generator_.addConnectivityConstraint(state, other_positions_dyn);
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