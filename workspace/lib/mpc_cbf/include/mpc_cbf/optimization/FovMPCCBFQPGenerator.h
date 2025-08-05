//
// Created by lishuo on 9/21/24.
//

#ifndef MPC_CBF_FOVMPCCBFQPGENERATOR_H
#define MPC_CBF_FOVMPCCBFQPGENERATOR_H

#include <mpc/optimization/PiecewiseBezierMPCQPGenerator.h>
#include <mpc_cbf/optimization/FovMPCCBFQPOperations.h>
#include <mpc_cbf/optimization/MPCCBFQPGeneratorBase.h>

namespace mpc_cbf {
template <typename T, unsigned int DIM>
class FovMPCCBFQPGenerator : public MPCCBFQPGeneratorBase<T, DIM> {
  public:
    using Base = MPCCBFQPGeneratorBase<T, DIM>;
    using FovMPCCBFQPOperations = mpc_cbf::FovMPCCBFQPOperations<T, DIM>;
    using PiecewiseBezierMPCQPGenerator = mpc::PiecewiseBezierMPCQPGenerator<T, DIM>;
    using PiecewiseBezierMPCQPOperations = mpc::PiecewiseBezierMPCQPOperations<T, DIM>;
    using typename Base::CostAddition;
    using typename Base::LinearConstraint;
    using typename Base::Row;
    using typename Base::State;
    using typename Base::Vector;

    // Initialization method
    void addPiecewise(std::unique_ptr<FovMPCCBFQPOperations>&& piecewise_mpc_cbf_operations_ptr,
                      int num_neighbors, bool slack_mode);

    // FoV-specific constraint methods
    void addSafetyCBFConstraint(const State& current_state, const Vector& other_pos,
                                T slack_value = 0);
    void addFovLBConstraint(const State& current_state, const Vector& other_pos, T slack_value = 0);
    void addFovRBConstraint(const State& current_state, const Vector& other_pos, T slack_value = 0);
    void addRangeCBFConstraint(const State& current_state, const Vector& other_pos,
                               T slack_value = 0);

    // FoV-specific predicted constraint methods
    void addPredSafetyCBFConstraints(const std::vector<State>& pred_states, const Vector& other_pos,
                                     const std::vector<T>& slack_values);
    void addPredFovLBConstraints(const std::vector<State>& pred_states, const Vector& other_pos,
                                 const std::vector<T>& slack_values);
    void addPredFovRBConstraints(const std::vector<State>& pred_states, const Vector& other_pos,
                                 const std::vector<T>& slack_values);
    void addPredRangeCBFConstraints(const std::vector<State>& pred_states, const Vector& other_pos,
                                    const std::vector<T>& slack_values);

    // FoV-specific constraint methods with slack variables
    void addSafetyCBFConstraintWithSlackVariables(const State& current_state,
                                                  const Vector& other_pos,
                                                  std::size_t neighbor_idx);
    void addFovLBConstraintWithSlackVariables(const State& current_state, const Vector& other_pos,
                                              std::size_t neighbor_idx);
    void addFovRBConstraintWithSlackVariables(const State& current_state, const Vector& other_pos,
                                              std::size_t neighbor_idx);
    void addRangeCBFConstraintWithSlackVariables(const State& current_state,
                                                 const Vector& other_pos, std::size_t neighbor_idx);

    // FoV-specific predicted constraint methods with slack variables
    void addPredSafetyCBFConstraintsWithSlackVariables(const std::vector<State>& pred_states,
                                                       const Vector& other_pos,
                                                       std::size_t neighbor_idx);
    void addPredFovLBConstraintsWithSlackVariables(const std::vector<State>& pred_states,
                                                   const Vector& other_pos,
                                                   std::size_t neighbor_idx);
    void addPredFovRBConstraintsWithSlackVariables(const std::vector<State>& pred_states,
                                                   const Vector& other_pos,
                                                   std::size_t neighbor_idx);
    void addPredRangeCBFConstraintsWithSlackVariables(const std::vector<State>& pred_states,
                                                      const Vector& other_pos,
                                                      std::size_t neighbor_idx);

  private:
    std::unique_ptr<FovMPCCBFQPOperations> piecewise_mpc_cbf_operations_ptr_;
};

} // namespace mpc_cbf

#endif // MPC_CBF_FOVMPCCBFQPGENERATOR_H
