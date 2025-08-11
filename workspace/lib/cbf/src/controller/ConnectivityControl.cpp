//
// Created on 4/10/25.
//

#include <cbf/controller/ConnectivityControl.h>

namespace cbf {
auto logger = common::initializeLogging();

template <typename T, unsigned int DIM>
ConnectivityControl<T, DIM>::ConnectivityControl(std::shared_ptr<ConnectivityCBF> cbf,
                                                 const SlackConfig& slack_config,
                                                 size_t num_neighbors)
    : qp_generator_(cbf, slack_config, num_neighbors),
      slack_config_(slack_config),
      num_neighbors_(num_neighbors),
      cbf_(cbf) {}

template <typename T, unsigned int DIM>
bool ConnectivityControl<T, DIM>::optimize(VectorDIM& cbf_u, const VectorDIM& desired_u,
                                           std::vector<State> current_states, size_t self_idx,
                                           const VectorDIM& u_min, const VectorDIM& u_max) {

    State current_state = current_states.at(self_idx);
    Vector state(2 * DIM);
    state << current_state.pos_, current_state.vel_;

    // Add cost terms
    qp_generator_.addDesiredControlCost(desired_u);
    
    // Add separate slack costs based on configuration
    if (slack_config_.safety_slack) {
        std::vector<T> safety_costs(num_neighbors_, slack_config_.safety_slack_cost);
        // Apply decay over horizon if needed
        for (size_t i = 0; i < safety_costs.size(); ++i) {
            safety_costs[i] *= pow(slack_config_.slack_decay_rate, i);
        }
        qp_generator_.addSlackCost(safety_costs, qp_generator_.safety_slack_variables_);
    }
    
    if (slack_config_.clf_slack) {
        std::vector<T> clf_costs(num_neighbors_, slack_config_.clf_slack_cost);
        for (size_t i = 0; i < clf_costs.size(); ++i) {
            clf_costs[i] *= pow(slack_config_.slack_decay_rate, i);
        }
        qp_generator_.addSlackCost(clf_costs, qp_generator_.clf_slack_variables_);
    }
    
    if (slack_config_.connectivity_slack) {
        std::vector<T> connectivity_costs = {slack_config_.connectivity_slack_cost};
        qp_generator_.addSlackCost(connectivity_costs, qp_generator_.connectivity_slack_variables_);
    }

    // add constraints

    // Add safety constraints
    for (size_t i = 0; i < num_neighbors_; ++i) {
        Vector neighbor_state(2 * DIM);
        neighbor_state << current_states.at(i + (i >= self_idx ? 1 : 0)).pos_,
            current_states.at(i + (i >= self_idx ? 1 : 0)).vel_;
        qp_generator_.addSafetyConstraint(state, neighbor_state, i);
    }

    // Add velocity and acceleration constraints
    qp_generator_.addMinVelConstraints(state);
    qp_generator_.addMaxVelConstraints(state);
    // qp_generator_.addControlBoundConstraint(u_min, u_max);

    Eigen::MatrixXd robot_states(current_states.size(), 2 * DIM);
    for (size_t i = 0; i < current_states.size(); ++i) {
        robot_states.row(i).head(DIM) = current_states[i].pos_.transpose();
        robot_states.row(i).tail(DIM) = current_states[i].vel_.transpose();
    }
    const auto robot_positions =
        robot_states.leftCols(2); // Extract only the position columns (x, y)
    auto [lambda2, eigenvec] = cbf_->getLambda2(robot_positions);
    // TODO: this 0.1 threshold is arbitrary
    if (lambda2 > 0.1) {
        qp_generator_.addConnConstraint(state, robot_states, self_idx);
    } else {
        // Add CLF constraints with proper indexing
        for (size_t i = 0; i < num_neighbors_; ++i) {
            Vector neighbor_state(2 * DIM);
            neighbor_state << current_states.at(i + (i >= self_idx ? 1 : 0)).pos_,
                current_states.at(i + (i >= self_idx ? 1 : 0)).vel_;
            qp_generator_.addCLFConstraint(state, neighbor_state, i);
        }
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

    cbf_u = qp_generator_.generatorCBFControlInput();
    return success;
}

// Explicit template instantiation
template class ConnectivityControl<double, 3U>;

} // namespace cbf