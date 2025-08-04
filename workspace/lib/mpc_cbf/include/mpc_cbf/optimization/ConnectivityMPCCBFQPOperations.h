//
// Created by yutong on 8/4/25.
//

#ifndef MPC_CBF_CONNECTIVITYMPCQPOPERATIONS_H
#define MPC_CBF_CONNECTIVITYMPCQPOPERATIONS_H

#include <cbf/detail/ConnectivityCBF.h>
#include <mpc/optimization/PiecewiseBezierMPCQPOperations.h>
#include <mpc_cbf/optimization/MPCCBFQPOperationsBase.h>

namespace mpc_cbf {
template <typename T, unsigned int DIM>
class ConnectivityMPCCBFQPOperations : public MPCCBFQPOperationsBase<T, DIM> {
  public:
    using Base = MPCCBFQPOperationsBase<T, DIM>;
    using ConnectivityCBF = cbf::ConnectivityCBF;
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

    struct Params {
        mpc::PiecewiseBezierParams<T, DIM>& piecewise_bezier_params;
        mpc::MPCParams<T>& mpc_params;
        cbf::ConnectivityCBFParams<T>& connectivity_cbf_params;
    };

    ConnectivityMPCCBFQPOperations(Params& p, std::shared_ptr<DoubleIntegrator> model_ptr,
                                   std::shared_ptr<ConnectivityCBF> connectivity_cbf_ptr);

    // Connectivity-specific constraint methods
    LinearConstraint safetyCBFConstraint(const Vector& current_state, const Vector& neighbor_state,
                                         T slack_value = 0);
    LinearConstraint connectivityConstraint(const Vector& current_state,
                                            const Vector& neighbor_state, T slack_value = 0);

    // Predicted constraints
    std::vector<LinearConstraint> predSafetyCBFConstraints(const std::vector<State>& pred_states,
                                                           const Vector& neighbor_state);
    std::vector<LinearConstraint> predConnectivityConstraints(const std::vector<State>& pred_states,
                                                              const Vector& neighbor_state);

    // Accessors
    std::unique_ptr<PiecewiseBezierMPCQPOperations> piecewise_mpc_operations_ptr();

  private:
    std::shared_ptr<ConnectivityCBF> connectivity_cbf_ptr_;
};

} // namespace mpc_cbf

#endif // MPC_CBF_CONNECTIVITYMPCQPOPERATIONS_H
