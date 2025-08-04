//
// Created by lishuo on 9/20/24.
//

#ifndef MPC_CBF_FOVMPCCBFQPOPERATIONS_H
#define MPC_CBF_FOVMPCCBFQPOPERATIONS_H

#include <mpc_cbf/optimization/MPCCBFQPOperationsBase.h>
#include <mpc/optimization/PiecewiseBezierMPCQPOperations.h>
#include <cbf/detail/FovCBF.h>

namespace mpc_cbf {
    template <typename T, unsigned int DIM>
    class FovMPCCBFQPOperations : public MPCCBFQPOperationsBase<T, DIM> {
    public:
        using Base = MPCCBFQPOperationsBase<T, DIM>;
        using PiecewiseBezierMPCQPOperations = mpc::PiecewiseBezierMPCQPOperations<T, DIM>;
        using DoubleIntegrator = typename PiecewiseBezierMPCQPOperations::DoubleIntegrator;
        using FovCBF = cbf::FovCBF;
        using QPOperation = qpcpp::QPOperations<T>;
        using CostAddition = typename QPOperation::CostAddition;
        using LinearConstraint = typename QPOperation::LinearConstraint;
        using typename Base::VectorDIM;
        using typename Base::Vector;
        using typename Base::State;
        using Row = math::Row<T>;
        using Matrix = math::Matrix<T>;


        struct Params {
            mpc::PiecewiseBezierParams<T, DIM> &piecewise_bezier_params;
            mpc::MPCParams<T> &mpc_params;
            cbf::FoVCBFParams<T> &fov_cbf_params;
        };

        FovMPCCBFQPOperations(Params &p, std::shared_ptr<DoubleIntegrator> model_ptr, std::shared_ptr<FovCBF> fov_cbf_ptr);

        // FoV-specific constraint methods
        LinearConstraint safetyCBFConstraint(const State& current_state, const Vector& other_pos, T slack_value=0);
        std::vector<LinearConstraint> fovLBConstraint(const State& current_state, const Vector& other_pos, T slack_value=0);
        std::vector<LinearConstraint> fovRBConstraint(const State& current_state, const Vector& other_pos, T slack_value=0);
        LinearConstraint rangeCBFConstraint(const State& current_state, const Vector& other_pos, T slack_value=0);

        // Predicted constraints
        std::vector<LinearConstraint> predSafetyCBFConstraints(const std::vector<State>& pred_states, const Vector& other_pos);
        std::vector<LinearConstraint> predFovLBConstraints(const std::vector<State>& pred_states, const Vector& other_pos);
        std::vector<LinearConstraint> predFovRBConstraints(const std::vector<State>& pred_states, const Vector& other_pos);
        std::vector<LinearConstraint> predRangeCBFConstraints(const std::vector<State>& pred_states, const Vector& other_pos);

    private:
        std::shared_ptr<FovCBF> fov_cbf_ptr_;
    };

} // mpc_cbf

#endif //MPC_CBF_FOVMPCCBFQPOPERATIONS_H
