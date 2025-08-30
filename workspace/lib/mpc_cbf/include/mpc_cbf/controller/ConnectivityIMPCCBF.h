//
// Created by yutong on 8/4/25.
//

#ifndef MPC_CBF_CONNECTIVITYIMPCCBF_H
#define MPC_CBF_CONNECTIVITYIMPCCBF_H

#include <cbf/detail/ConnectivityCBF.h>
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
        cbf::SlackConfig slack_config_;
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
    // ===[ 统计结果 Getter ]===
    size_t failTotal() const        { return fail_cnt_total_; }
    size_t failConn()  const        { return fail_cnt_connectivity_; }
    size_t failCLF()   const        { return fail_cnt_clf_; }
    T      lastLambda2() const      { return last_lambda2_; }
    const std::vector<std::pair<Vector, T>>& samples_conn() const { return samples_conn_; }
    const std::vector<std::pair<Vector, T>>& samples_clf()  const { return samples_clf_;  }
    std::size_t failSafety() const { return fail_cnt_safety_; }

    // 和 samples_conn()/samples_clf() 相同风格，返回 (a, b)
    const std::vector<std::pair<Vector, T>>& samples_safety() const { return samples_safety_; }

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

    // ===[ A 法：失败归因统计 ]=========================================
    enum class ConstraintMode {
        None = 0,
        Safety,
        Connectivity,
        CLF
    };

    ConstraintMode current_mode_ = ConstraintMode::None;

    // 累计统计
    size_t fail_cnt_total_        = 0;
    size_t fail_cnt_connectivity_ = 0;
    size_t fail_cnt_clf_          = 0;
    std::size_t fail_cnt_safety_{0};
    std::vector<std::pair<Vector, T>> samples_safety_;
    T last_lambda2_ = T(0);
    static constexpr size_t kMaxSamples = 5;

    // 简单容器：只存 ac/bc；如需同时存 robot_idx/sim_t/loop_idx，可扩展结构体
    std::vector<std::pair<Vector, T>> samples_conn_; // 前 5 条 connectivity 失败的 (ac, bc)
    std::vector<std::pair<Vector, T>> samples_clf_;  // 前 5 条 CLF 失败的 (a,  b)
    cbf::SlackConfig slack_config_;
};

} // namespace mpc_cbf

#endif // MPC_CBF_CONNECTIVITYIMPCCBF_H