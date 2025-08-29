//
// Created by yutong on 8/4/25.
//

#ifndef MPC_CBF_CONNECTIVITYMPCQPOPERATIONS_H
#define MPC_CBF_CONNECTIVITYMPCQPOPERATIONS_H

#include <mpc_cbf/detail/ConnectivityCBF.h>
#include <mpc/optimization/PiecewiseBezierMPCQPOperations.h>
#include <mpc_cbf/optimization/MPCCBFQPOperationsBase.h>

namespace mpc_cbf {
template <typename T, unsigned int DIM>
class ConnectivityMPCCBFQPOperations : public MPCCBFQPOperationsBase<T, DIM> {
  public:
    using Base = MPCCBFQPOperationsBase<T, DIM>;
    using ConnectivityCBF = mpc_cbf::ConnectivityCBF;
    using PiecewiseBezierMPCQPOperations = mpc::PiecewiseBezierMPCQPOperations<T, DIM>;
    using DoubleIntegrator = typename PiecewiseBezierMPCQPOperations::DoubleIntegrator;
    using typename Base::CostAddition;
    using typename Base::LinearConstraint;
    using typename Base::Matrix;
    using typename Base::QPOperation;
    using typename Base::State;
    using typename Base::Vector;
    using typename Base::VectorDIM;
    using Row = math::Row<T>;
    const Vector& debugLastConn_ac() const { return last_conn_ac_; }
    T              debugLastConn_bc() const { return last_conn_bc_; }
    const Vector& debugLastCLF_a()   const { return last_clf_a_;   }
    T              debugLastCLF_b()   const { return last_clf_b_;   }


    struct Params {
        mpc::PiecewiseBezierParams<T, DIM>& piecewise_bezier_params;
        mpc::MPCParams<T>& mpc_params;
        mpc_cbf::ConnectivityCBFParams<T>& connectivity_cbf_params;
    };

    ConnectivityMPCCBFQPOperations(Params& p, std::shared_ptr<DoubleIntegrator> model_ptr,
                                   std::shared_ptr<ConnectivityCBF> connectivity_cbf_ptr);

    // Connectivity-specific constraint methods
    LinearConstraint safetyCBFConstraint(const Vector& current_state, const Vector& neighbor_state,
                                         T slack_value = 0);
    LinearConstraint connectivityConstraint(const Vector& x_self, const std::vector<VectorDIM>& other_positions,
                                          T slack_value = 0);
    LinearConstraint clfConstraint(const Vector& current_state, const Vector& neighbor_state,
                                   T slack_value = 0);

    // Predicted constraints
    std::vector<LinearConstraint> predSafetyCBFConstraints(const std::vector<State>& pred_states,
                                                           const Vector& neighbor_state);
    std::vector<LinearConstraint> predConnectivityConstraints(const std::vector<State>& pred_states,
                                                            const std::vector<VectorDIM>& other_positions);
    std::vector<LinearConstraint> predCLFConstraints(const std::vector<State>& pred_states,
                                                     const Vector& neighbor_state);

    // Accessors
    std::unique_ptr<PiecewiseBezierMPCQPOperations> piecewise_mpc_operations_ptr();
    std::shared_ptr<ConnectivityCBF> connectivityCBF() const;
    const Vector& debugLastSafety_a() const { return last_safety_a_; }
    T             debugLastSafety_b() const { return last_safety_b_; }

    // 预测版（可选）
    const std::vector<Vector>& debugLastPredSafety_Ak() const { return last_pred_safety_Ak_; }
    const std::vector<T>&      debugLastPredSafety_bk() const { return last_pred_safety_bk_; }

  private:
    std::shared_ptr<ConnectivityCBF> connectivity_cbf_ptr_;
    Vector last_conn_ac_;  // 最近一次 connectivity 的 Ac
    T      last_conn_bc_{0};
    Vector last_clf_a_;    // 最近一次 CLF 的 a
    T      last_clf_b_{0};
    Vector last_safety_a_;
    T      last_safety_b_{};

    // 预测版（逐步）可选
    std::vector<Vector> last_pred_safety_Ak_; // 每个 k 的 Ak（维度 DIM*k_hor）
    std::vector<T>      last_pred_safety_bk_; // 每个 k 的 bk
};

} // namespace mpc_cbf

#endif // MPC_CBF_CONNECTIVITYMPCQPOPERATIONS_H
