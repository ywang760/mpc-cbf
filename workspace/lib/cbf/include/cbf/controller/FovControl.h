//
// Created by lishuo on 9/1/24.
//

#ifndef CBF_FOVCONTROL_H
#define CBF_FOVCONTROL_H

#include <cbf/detail/FovCBF.h>
#include <cbf/optimization/FovQPGenerator.h>
#include <model/DoubleIntegrator.h>
#include <qpcpp/solvers/CPLEX.h>
#include <numeric>

namespace cbf {
    template <typename T, unsigned int DIM>
    class FovControl
    {
    public:
        using QPGenerator = cbf::FovQPGenerator<T, DIM>;
        using Problem = qpcpp::Problem<T>;
        using CPLEXSolver = qpcpp::CPLEXSolver<T>;
        using SolveStatus = qpcpp::SolveStatus;
        using Vector = math::Vector<T>;
        using VectorDIM = math::VectorDIM<T, DIM>;
        using Matrix = math::Matrix<T>;
        using State = model::State<T, DIM>;

        FovControl(std::shared_ptr<FovCBF> cbf, int num_robots = 0, bool slack_mode = false, T slack_cost = 1000, T slack_decay_rate = 0.2);
        ~FovControl() = default;
        bool optimize(VectorDIM &cbf_u, const VectorDIM &desired_u,
                      const State &current_state,
                      const std::vector<VectorDIM> &other_robot_positions,
                      const std::vector<Matrix> &other_robot_covs,
                      const VectorDIM& u_min, const VectorDIM& u_max);

        T distanceToEllipse(const VectorDIM& robot_position, const Vector& target_mean, const Matrix& target_cov);
        bool compareDist(const VectorDIM& p_current, const std::pair<VectorDIM, Matrix>& a, const std::pair<VectorDIM, Matrix>& b);

    private:
        QPGenerator qp_generator_;
        bool slack_mode_;
        T slack_cost_;
        T slack_decay_rate_;
    };

} // cbf

#endif //CBF_FOVCONTROL_H
