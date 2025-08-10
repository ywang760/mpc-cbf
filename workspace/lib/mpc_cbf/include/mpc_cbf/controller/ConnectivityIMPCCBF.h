//
// Created by yutong on 8/4/25.
//

#ifndef MPC_CBF_CONNECTIVITYIMPCCBF_H
#define MPC_CBF_CONNECTIVITYIMPCCBF_H

#include <math/Helpers.h>
#include <math/collision_shapes/CollisionShape.h>
#include <mpc_cbf/optimization/ConnectivityMPCCBFQPGenerator.h>
#include <numeric>
#include <qpcpp/Problem.h>
#include <qpcpp/solvers/CPLEX.h>

namespace mpc_cbf {
template <typename T, unsigned int DIM>
class ConnectivityIMPCCBF {
  public:
    using MPCCBFParams = typename mpc_cbf::ConnectivityMPCCBFQPOperations<T, DIM>::Params;
    using DoubleIntegrator =
        typename mpc_cbf::ConnectivityMPCCBFQPOperations<T, DIM>::DoubleIntegrator;
    using ConnectivityMPCCBFQPOperations = mpc_cbf::ConnectivityMPCCBFQPOperations<T, DIM>;
    using ConnectivityMPCCBFQPGenerator = mpc_cbf::ConnectivityMPCCBFQPGenerator<T, DIM>;
    using ConnectivityCBF =
        typename mpc_cbf::ConnectivityMPCCBFQPOperations<T, DIM>::ConnectivityCBF;
    using State = typename mpc_cbf::ConnectivityMPCCBFQPOperations<T, DIM>::State;
    using TuningParams = mpc::TuningParams<T>;

    using SingleParameterPiecewiseCurve = splines::SingleParameterPiecewiseCurve<T, DIM>;
    using VectorDIM = math::VectorDIM<T, DIM>;
    using Vector = math::Vector<T>;

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

    ConnectivityIMPCCBF(Params& p, std::shared_ptr<DoubleIntegrator> model_ptr,
                        std::shared_ptr<ConnectivityCBF> connectivity_cbf_ptr,
                        std::shared_ptr<const CollisionShape> collision_shape_ptr,
                        int num_neighbors = 0);
    ~ConnectivityIMPCCBF() = default;

    bool optimize(std::vector<SingleParameterPiecewiseCurve>& result_curve,
                  const std::vector<State>& current_states, size_t self_idx,
                  const Vector& ref_positions);

    void resetProblem();

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

    ConnectivityMPCCBFQPGenerator qp_generator_;
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

#endif // MPC_CBF_CONNECTIVITYIMPCCBF_H