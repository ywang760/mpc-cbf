//
// Created by yutong on 8/4/25.
//

#ifndef MPC_CBF_CONNECTIVITYMPCQPGENERATOR_H
#define MPC_CBF_CONNECTIVITYMPCQPGENERATOR_H

#include <mpc/optimization/PiecewiseBezierMPCQPGenerator.h>
#include <mpc_cbf/optimization/ConnectivityMPCCBFQPOperations.h>
#include <mpc_cbf/optimization/MPCCBFQPGeneratorBase.h>

namespace mpc_cbf {
template <typename T, unsigned int DIM>
class ConnectivityMPCCBFQPGenerator : public MPCCBFQPGeneratorBase<T, DIM> {
  public:
    using Base = MPCCBFQPGeneratorBase<T, DIM>;
    using ConnectivityMPCCBFQPOperations = mpc_cbf::ConnectivityMPCCBFQPOperations<T, DIM>;
    using PiecewiseBezierMPCQPGenerator = mpc::PiecewiseBezierMPCQPGenerator<T, DIM>;
    using PiecewiseBezierMPCQPOperations = mpc::PiecewiseBezierMPCQPOperations<T, DIM>;
    using typename Base::CostAddition;
    using typename Base::LinearConstraint;
    using typename Base::Row;
    using typename Base::State;
    using typename Base::Vector;
    using typename Base::VectorDIM;

    // Constructor
    ConnectivityMPCCBFQPGenerator(
        std::unique_ptr<ConnectivityMPCCBFQPOperations>&& piecewise_mpc_cbf_operations_ptr,
        int num_neighbors, bool slack_mode);

    // Connectivity-specific constraint methods
    void addSafetyCBFConstraint(const Vector& current_state, const Vector& neighbor_state,
                                std::size_t neighbor_idx, T slack_value = 0);
    void addConnectivityConstraint(const Vector& x_self, const std::vector<VectorDIM>& other_positions,
                                  T slack_value = 0);
    void addCLFConstraint(const Vector& current_state, const Vector& neighbor_state,
                          std::size_t neighbor_idx, T slack_value = 0);

    // Connectivity-specific predicted constraint methods
    void addPredSafetyCBFConstraints(const std::vector<State>& pred_states,
                                     const Vector& neighbor_state, std::size_t neighbor_idx);
    void addPredConnectivityConstraints(const std::vector<State>& pred_states,
                                      const std::vector<VectorDIM>& other_positions,
                                      T slack_value = 0);
    void addPredCLFConstraints(const std::vector<State>& pred_states, const Vector& neighbor_state,
                               std::size_t neighbor_idx);

    // Accessor for connectivity CBF
    std::shared_ptr<typename ConnectivityMPCCBFQPOperations::ConnectivityCBF>
    connectivityCBF() const;

  private:
    std::unique_ptr<ConnectivityMPCCBFQPOperations> piecewise_mpc_cbf_operations_ptr_;
    bool slack_mode_ = false;
};

} // namespace mpc_cbf

#endif // MPC_CBF_CONNECTIVITYMPCQPGENERATOR_H
