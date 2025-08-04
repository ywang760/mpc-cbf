//
// Created on 8/4/25.
//

#ifndef MPC_CBF_MPCCBFQPOPERATIONSBASE_H
#define MPC_CBF_MPCCBFQPOPERATIONSBASE_H

#include <mpc/optimization/PiecewiseBezierMPCQPOperations.h>
#include <memory>

namespace mpc_cbf {
    template <typename T, unsigned int DIM>
    class MPCCBFQPOperationsBase {
    public:
        // Common type definitions
        using PiecewiseBezierMPCQPOperations = mpc::PiecewiseBezierMPCQPOperations<T, DIM>;
        using DoubleIntegrator = typename PiecewiseBezierMPCQPOperations::DoubleIntegrator;
        using QPOperation = qpcpp::QPOperations<T>;
        using CostAddition = typename QPOperation::CostAddition;
        using LinearConstraint = typename QPOperation::LinearConstraint;
        using State = model::State<T, DIM>;
        using Vector = math::Vector<T>;
        using VectorDIM = math::VectorDIM<T, DIM>;
        using Matrix = math::Matrix<T>;
        using Row = math::Row<T>;
        
        // Virtual destructor for proper inheritance
        virtual ~MPCCBFQPOperationsBase() = default;
        
        // Common interface (concrete implementations)
        CostAddition slackCost(const std::vector<double>& slack_weights);
        std::unique_ptr<PiecewiseBezierMPCQPOperations> piecewise_mpc_operations_ptr();
        
    protected:
        // Constructor for derived classes
        MPCCBFQPOperationsBase(std::shared_ptr<DoubleIntegrator> model_ptr);
        
        // Common members accessible to derived classes
        std::unique_ptr<PiecewiseBezierMPCQPOperations> piecewise_mpc_operations_ptr_;
        std::shared_ptr<DoubleIntegrator> model_ptr_;
        T h_;
        int k_hor_;
        mpc::TuningParams<T> mpc_tuning_;
        Matrix U_basis_;
    };

} // mpc_cbf

#endif //MPC_CBF_MPCCBFQPOPERATIONSBASE_H