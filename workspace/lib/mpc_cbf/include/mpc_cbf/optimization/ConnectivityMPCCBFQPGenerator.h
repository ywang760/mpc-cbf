//
// Created by yutong on 8/4/25.
//

#ifndef MPC_CBF_CONNECTIVITYMPCQPGENERATOR_H
#define MPC_CBF_CONNECTIVITYMPCQPGENERATOR_H

#include <mpc_cbf/optimization/MPCCBFQPGeneratorBase.h>
#include <mpc_cbf/optimization/ConnectivityMPCCBFQPOperations.h>
#include <mpc/optimization/PiecewiseBezierMPCQPGenerator.h>

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
        using typename Base::State;
        using typename Base::Vector;
        using typename Base::Row;
        
        // Initialization method
        void addPiecewise(std::unique_ptr<ConnectivityMPCCBFQPOperations>&& piecewise_mpc_cbf_operations_ptr, int num_neighbors, bool slack_mode);
        
        // Connectivity-specific constraint methods
        void addSafetyCBFConstraint(const State& current_state, const Vector& other_pos, T slack_value=0);
        void addConnectivityConstraint(const State& current_state, const Vector& other_pos, T slack_value=0);

        // Connectivity-specific predicted constraint methods
        void addPredSafetyCBFConstraints(const std::vector<State>& pred_states, const Vector& other_pos, const std::vector<T>& slack_values);
        void addPredConnectivityConstraints(const std::vector<State>& pred_states, const Vector& other_pos, const std::vector<T>& slack_values);

        // Connectivity-specific constraint methods with slack variables
        void addSafetyCBFConstraintWithSlackVariables(const State& current_state, const Vector& other_pos, std::size_t neighbor_idx);
        void addConnectivityConstraintWithSlackVariables(const State& current_state, const Vector& other_pos, std::size_t neighbor_idx);

        // Connectivity-specific predicted constraint methods with slack variables
        void addPredSafetyCBFConstraintsWithSlackVariables(const std::vector<State>& pred_states, const Vector& other_pos, std::size_t neighbor_idx);
        void addPredConnectivityConstraintsWithSlackVariables(const std::vector<State>& pred_states, const Vector& other_pos, std::size_t neighbor_idx);

    private:
        std::unique_ptr<ConnectivityMPCCBFQPOperations> piecewise_mpc_cbf_operations_ptr_;
    };

} // mpc_cbf

#endif //MPC_CBF_CONNECTIVITYMPCQPGENERATOR_H
