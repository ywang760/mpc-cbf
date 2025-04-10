//
// Created by lishuo on 9/20/24.
//

#ifndef MPC_PIECEWISEBEZIERMPCCBFQPOPERATIONS_H
#define MPC_PIECEWISEBEZIERMPCCBFQPOPERATIONS_H

#include <mpc/optimization/PiecewiseBezierMPCQPOperations.h>
#include <cbf/detail/FovCBF.h>

namespace mpc_cbf {
    template <typename T, unsigned int DIM>
    class PiecewiseBezierMPCCBFQPOperations {
    public:
        using PiecewiseBezierMPCQPOperations = mpc::PiecewiseBezierMPCQPOperations<T, DIM>;
        using DoubleIntegrator = typename PiecewiseBezierMPCQPOperations::DoubleIntegrator;
        using FovCBF = cbf::FovCBF;
        using QPOperation = qpcpp::QPOperations<T>;
        using CostAddition = typename QPOperation::CostAddition;
        using LinearConstraint = typename QPOperation::LinearConstraint;
        using VectorDIM = math::VectorDIM<T, DIM>;
        using Vector = math::Vector<T>;
        using Row = math::Row<T>;
        using Matrix = math::Matrix<T>;
        using State = model::State<T, DIM>;


        struct Params {
            mpc::PiecewiseBezierParams<T, DIM> &piecewise_bezier_params;
            mpc::MPCParams<T> &mpc_params;
            cbf::FoVCBFParams<T> &fov_cbf_params;
        };

        PiecewiseBezierMPCCBFQPOperations(Params &p, std::shared_ptr<DoubleIntegrator> model_ptr, std::shared_ptr<FovCBF> fov_cbf_ptr);
        std::unique_ptr<PiecewiseBezierMPCQPOperations> piecewise_mpc_operations_ptr();

        CostAddition slackCost(const std::vector<double>& slack_weights);
        LinearConstraint safetyCBFConstraint(const State& current_state, const Vector& other_pos, T slack_value=0);
        std::vector<LinearConstraint> fovLBConstraint(const State& current_state, const Vector& other_pos, T slack_value=0);
        std::vector<LinearConstraint> fovRBConstraint(const State& current_state, const Vector& other_pos, T slack_value=0);
        LinearConstraint rangeCBFConstraint(const State& current_state, const Vector& other_pos, T slack_value=0);

        std::vector<LinearConstraint> predSafetyCBFConstraints(const std::vector<State>& pred_states, const Vector& other_pos);
        std::vector<LinearConstraint> predFovLBConstraints(const std::vector<State>& pred_states, const Vector& other_pos);
        std::vector<LinearConstraint> predFovRBConstraints(const std::vector<State>& pred_states, const Vector& other_pos);
        std::vector<LinearConstraint> predRangeCBFConstraints(const std::vector<State>& pred_states, const Vector& other_pos);

    private:
        std::unique_ptr<PiecewiseBezierMPCQPOperations> piecewise_mpc_operations_ptr_;
        std::shared_ptr<DoubleIntegrator> model_ptr_;
        std::shared_ptr<FovCBF> fov_cbf_ptr_;

        // mpc params
        T h_;
        int k_hor_;
        mpc::TuningParams<T> mpc_tuning_;
        // model control predict
        Matrix U_basis_;
    };

} // mpc_cbf

#endif //MPC_PIECEWISEBEZIERMPCCBFQPOPERATIONS_H
