//
// Created by yutong on 8/4/25.
//

#ifndef MPC_CBF_CONNECTIVITYMPCQPGENERATOR_H
#define MPC_CBF_CONNECTIVITYMPCQPGENERATOR_H

#include <cbf/detail/ConnectivityCBF.h>
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

    // Constructor
    ConnectivityMPCCBFQPGenerator(
        std::unique_ptr<ConnectivityMPCCBFQPOperations>&& piecewise_mpc_cbf_operations_ptr,
        int num_neighbors, const cbf::SlackConfig& slack_config);

    // Connectivity-specific constraint methods
    void addSafetyCBFConstraint(const Vector& current_state, const Vector& neighbor_state,
                                std::size_t neighbor_idx, T slack_value = 0);
    void addConnectivityConstraint(const Eigen::MatrixXd& robot_states, size_t self_idx,
                                   T slack_value = 0);
    void addCLFConstraint(const Vector& current_state, const Vector& neighbor_state,
                          std::size_t neighbor_idx, T slack_value = 0);

    // Connectivity-specific predicted constraint methods
    void addPredSafetyCBFConstraints(const std::vector<State>& pred_states,
                                     const Vector& neighbor_state, std::size_t neighbor_idx);
    void addPredConnectivityConstraints(const std::vector<State>& pred_states,
                                        const Eigen::MatrixXd& robot_states, size_t self_idx,
                                        const std::vector<T>& slack_values);
    void addPredCLFConstraints(const std::vector<State>& pred_states, const Vector& neighbor_state,
                               std::size_t neighbor_idx);

    // Accessor for connectivity CBF
    std::shared_ptr<typename ConnectivityMPCCBFQPOperations::ConnectivityCBF>
    connectivityCBF() const;

  public:
    // Public slack variable pools for derived class access
    cbf::SlackConfig slack_config_;
    std::vector<qpcpp::Variable<T>*> safety_slack_variables_;    // n-1 variables
    std::vector<qpcpp::Variable<T>*> clf_slack_variables_;       // n-1 variables  
    std::vector<qpcpp::Variable<T>*> connectivity_slack_variables_; // 1 variable

  private:
    std::unique_ptr<ConnectivityMPCCBFQPOperations> piecewise_mpc_cbf_operations_ptr_;
};

} // namespace mpc_cbf

#endif // MPC_CBF_CONNECTIVITYMPCQPGENERATOR_H
