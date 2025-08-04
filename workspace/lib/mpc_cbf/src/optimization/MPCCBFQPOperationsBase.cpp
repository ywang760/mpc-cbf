//
// Created on 8/4/25.
//

#include <mpc_cbf/optimization/MPCCBFQPOperationsBase.h>

namespace mpc_cbf {
    template <typename T, unsigned int DIM>
    MPCCBFQPOperationsBase<T, DIM>::MPCCBFQPOperationsBase(std::shared_ptr<DoubleIntegrator> model_ptr)
        : model_ptr_(model_ptr) {
        // Base class constructor - derived classes will initialize piecewise_mpc_operations_ptr_
    }

    template <typename T, unsigned int DIM>
    std::unique_ptr<typename MPCCBFQPOperationsBase<T, DIM>::PiecewiseBezierMPCQPOperations>
    MPCCBFQPOperationsBase<T, DIM>::piecewise_mpc_operations_ptr() {
        return std::move(piecewise_mpc_operations_ptr_);
    }

    template <typename T, unsigned int DIM>
    typename MPCCBFQPOperationsBase<T, DIM>::CostAddition
    MPCCBFQPOperationsBase<T, DIM>::slackCost(const std::vector<double>& slack_weights) {
        Matrix quadratic_term(slack_weights.size(), slack_weights.size());
        quadratic_term.setZero();
        Vector linear_term(slack_weights.size());
        linear_term.setZero();

        for (std::size_t i = 0; i < slack_weights.size(); ++i) {
            linear_term(i) = slack_weights.at(i);
        }
        return CostAddition(quadratic_term, linear_term, 0);
    }

    // Explicit template instantiation
    template class MPCCBFQPOperationsBase<double, 3U>;

} // mpc_cbf