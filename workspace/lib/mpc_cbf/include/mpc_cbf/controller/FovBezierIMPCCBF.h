//
// Created by lishuo on 9/22/24.
//

#ifndef MPC_CBF_FOVBEZIERIMPCCBF_H
#define MPC_CBF_FOVBEZIERIMPCCBF_H

#include <math/Helpers.h>
#include <math/collision_shapes/CollisionShape.h>
#include <mpc_cbf/optimization/FovMPCCBFQPGenerator.h>
#include <numeric>
#include <qpcpp/Problem.h>
#include <qpcpp/solvers/CPLEX.h>
#include <separating_hyperplanes/Voronoi.h>

namespace mpc_cbf {
template <typename T, unsigned int DIM>
class FovBezierIMPCCBF {
  public:
    using MPCCBFParams = typename mpc_cbf::FovMPCCBFQPOperations<T, DIM>::Params;
    using DoubleIntegrator = typename mpc_cbf::FovMPCCBFQPOperations<T, DIM>::DoubleIntegrator;
    using FovMPCCBFQPOperations = mpc_cbf::FovMPCCBFQPOperations<T, DIM>;
    using FovMPCCBFQPGenerator = mpc_cbf::FovMPCCBFQPGenerator<T, DIM>;
    using FovCBF = typename mpc_cbf::FovMPCCBFQPOperations<T, DIM>::FovCBF;
    using State = typename mpc_cbf::FovMPCCBFQPOperations<T, DIM>::State;
    using TuningParams = mpc::TuningParams<T>;

    using SingleParameterPiecewiseCurve = splines::SingleParameterPiecewiseCurve<T, DIM>;
    using VectorDIM = math::VectorDIM<T, DIM>;
    using Vector = math::Vector<T>;
    using AlignedBox = math::AlignedBox<T, DIM>;
    using Hyperplane = math::Hyperplane<T, DIM>;
    using Matrix = math::Matrix<T>;
    using CollisionShape = math::CollisionShape<T, DIM>;
    using Problem = qpcpp::Problem<T>;
    using CPLEXSolver = qpcpp::CPLEXSolver<T>;
    using SolveStatus = qpcpp::SolveStatus;

    struct IMPCParams {
        int cbf_horizon_;
        int impc_iter_;
        T slack_cost_;
        T slack_decay_rate_;
        bool slack_mode_;
    };

    struct Params {
        MPCCBFParams& mpc_cbf_params;
        IMPCParams& impc_params;
    };

    FovBezierIMPCCBF(Params& p, std::shared_ptr<DoubleIntegrator> model_ptr,
                     std::shared_ptr<FovCBF> fov_cbf_ptr, uint64_t bezier_continuity_upto_degree,
                     std::shared_ptr<const CollisionShape> collision_shape_ptr,
                     int num_neighbors = 0);
    ~FovBezierIMPCCBF() = default;

    bool optimize(std::vector<SingleParameterPiecewiseCurve>& result_curve,
                  const State& current_state, const std::vector<VectorDIM>& other_robot_positions,
                  const std::vector<Matrix>& other_robot_covs, const Vector& ref_positions);

    T distanceToEllipse(const VectorDIM& robot_position, const Vector& target_mean,
                        const Matrix& target_cov);

    bool compareDist(const VectorDIM& p_current, const std::pair<VectorDIM, Matrix>& a,
                     const std::pair<VectorDIM, Matrix>& b);

    void resetProblem();

    Vector generatorDerivativeControlInputs(uint64_t derivative_degree);

  private:
    TuningParams mpc_tuning_;
    VectorDIM v_min_;
    VectorDIM v_max_;
    VectorDIM a_min_;
    VectorDIM a_max_;

    Vector ts_samples_;
    T h_;
    int k_hor_;
    Vector h_samples_;

    FovMPCCBFQPGenerator qp_generator_;
    // the derivative degree that the resulting trajectory must be
    // continuous upto.
    uint64_t bezier_continuity_upto_degree_;
    std::shared_ptr<const CollisionShape> collision_shape_ptr_;
    int impc_iter_;
    int cbf_horizon_;
    T slack_cost_;
    T slack_decay_rate_;
    bool slack_mode_;
};

} // namespace mpc_cbf

#endif // MPC_CBF_BEZIERIMPCCBF_H
