//
// Created by lishuo on 9/1/24.
//

#include <cbf/controller/FovControl.h>

namespace cbf {
    template <typename T, unsigned int DIM>
    FovControl<T, DIM>::FovControl(std::shared_ptr<FovCBF> cbf, int number_neighbors, bool slack_mode, T slack_cost, T slack_decay_rate)
        : qp_generator_(cbf, number_neighbors, slack_mode), slack_mode_(slack_mode), slack_cost_(slack_cost), slack_decay_rate_(slack_decay_rate)
    {
    }

    template <typename T, unsigned int DIM>
    bool FovControl<T, DIM>::optimize(VectorDIM& cbf_u,
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
            std::vector<std::pair<VectorDIM, Matrix>> other_pos_covs(num_neighbors);
            for (size_t i = 0; i < num_neighbors; ++i) {
                other_pos_covs[i] = std::make_pair(other_robot_positions.at(i), other_robot_covs.at(i));
            }
            std::vector<size_t> idx(num_neighbors);
            std::iota(idx.begin(), idx.end(), 0);
            // sort
            std::sort(idx.begin(), idx.end(), [this, current_pos, other_pos_covs](size_t a_idx, size_t b_idx) {
                return compareDist(current_pos, other_pos_covs.at(a_idx), other_pos_covs.at(b_idx));
            });
            // define slack weights
            T w_init = slack_cost_;
            T decay_factor = slack_decay_rate_;
            for (size_t i = 0; i < num_neighbors; ++i) {
                size_t sort_idx = idx.at(i);
                slack_weights.at(i) = w_init * pow(decay_factor, sort_idx);
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
                qp_generator_.addLeftBorderConstraint(state, other_xy);
                qp_generator_.addRightBorderConstraint(state, other_xy);
                qp_generator_.addRangeConstraint(state, other_xy);
            } else {
                qp_generator_.addSafetyConstraint(state, other_xy, true, i);
                qp_generator_.addLeftBorderConstraint(state, other_xy, true, i);
                qp_generator_.addRightBorderConstraint(state, other_xy, true, i);
                qp_generator_.addRangeConstraint(state, other_xy, true, i);
            }
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
    T FovControl<T, DIM>::distanceToEllipse(const VectorDIM &robot_pos,
                                            const Vector &target_mean,
                                            const Matrix &target_cov) {
        assert(DIM == 3); // DIM other than 3 is not implemented yet.

        if (!isinf(target_cov(0, 0))) {
            Eigen::EigenSolver<Matrix> es(target_cov.block(0, 0, DIM-1, DIM-1));
            Vector eigenvalues = es.eigenvalues().real();
            Matrix eigenvectors = es.eigenvectors().real();

            // s = 4.605 for 90% confidence interval
            // s = 5.991 for 95% confidence interval
            // s = 9.210 for 99% confidence interval
            T s = 4.605;
            T a = sqrt(s * eigenvalues(0)); // major axis
            T b = sqrt(s * eigenvalues(1)); // minor axis

            // a could be smaller than b, so swap them
            if (a < b)
            {
                T temp = a;
                a = b;
                b = temp;
            }

            int m = 0; // higher eigenvalue index
            int l = 1; // lower eigenvalue index
            if (eigenvalues(1) > eigenvalues(0))
            {
                m = 1;
                l = 0;
            }

            T theta = atan2(eigenvectors(1, m), eigenvectors(0, m)); // angle of the major axis wrt positive x-asis (ccw rotation)
            if (theta < 0.0) {
                theta += M_PI;
            } // angle in [0, 2pi]

            T slope = atan2(-target_mean(1) + robot_pos(1), -target_mean(0) + robot_pos(0));
            T x_n = target_mean(0) + a * cos(slope - theta) * cos(theta) - b * sin(slope - theta) * sin(theta);
            T y_n = target_mean(1) + a * cos(slope - theta) * sin(theta) + b * sin(slope - theta) * cos(theta);

            Vector p_near(2);
            p_near << x_n, y_n;

            T dist = sqrt(pow(p_near(0) - robot_pos(0), 2) + pow(p_near(1) - robot_pos(1), 2));

            if (isnan(dist)) {
                dist = 5.0;
                return dist;
            }

            // Check if robot is inside ellipse
            T d = sqrt(pow(target_mean(0) - robot_pos(0), 2) + pow(target_mean(1) - robot_pos(1), 2));
            T range = sqrt(pow(target_mean(0) - p_near(0), 2) + pow(target_mean(1) - p_near(1), 2));

            if (d < range) {
                return -dist;
            } else {
                return dist;
            }
        }

        return -5.0;
    }


    template <typename T, unsigned int DIM>
    bool FovControl<T, DIM>::compareDist(const VectorDIM& p_current,
                                         const std::pair<VectorDIM, Matrix>& a,
                                         const std::pair<VectorDIM, Matrix>& b) {
        return distanceToEllipse(p_current, a.first, a.second) < distanceToEllipse(p_current, b.first, b.second);
    }

    template class FovControl<double, 3U>;

} // cbf