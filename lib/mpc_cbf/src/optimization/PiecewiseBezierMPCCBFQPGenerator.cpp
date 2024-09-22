//
// Created by lishuo on 9/21/24.
//

#include <mpc_cbf/optimization/PiecewiseBezierMPCCBFQPGenerator.h>

namespace mpc_cbf {
    template <typename T, unsigned int DIM>
    void PiecewiseBezierMPCCBFQPGenerator<T, DIM>::addPiecewise(
            std::unique_ptr<PiecewiseBezierMPCCBFQPOperations> &&piecewise_mpc_cbf_operations_ptr) {
        // init the PiecewiseBezierMPCQPGenerator API
        std::unique_ptr<PiecewiseBezierMPCQPOperations> piecewise_mpc_operations_ptr = piecewise_mpc_cbf_operations_ptr->piecewise_mpc_operations_ptr();
        piecewise_mpc_qp_generator_ptr_->addPiecewise(std::move(piecewise_mpc_operations_ptr));
        piecewise_mpc_cbf_operations_ptr_ = std::move(piecewise_mpc_cbf_operations_ptr);
    }

    template <typename T, unsigned int DIM>
    qpcpp::Problem<T> &PiecewiseBezierMPCCBFQPGenerator<T, DIM>::problem() {return piecewise_mpc_qp_generator_ptr_->problem();}

    template <typename T, unsigned int DIM>
    std::shared_ptr<typename PiecewiseBezierMPCCBFQPGenerator<T, DIM>::PiecewiseBezierMPCQPGenerator>
    PiecewiseBezierMPCCBFQPGenerator<T, DIM>::piecewise_mpc_qp_generator_ptr() {return piecewise_mpc_qp_generator_ptr_;}

    template <typename T, unsigned int DIM>
    void PiecewiseBezierMPCCBFQPGenerator<T, DIM>::addSafetyCBFConstraint(const State &current_state,
                                                                          const Vector &other_pos) {
        LinearConstraint linear_constraint = piecewise_mpc_cbf_operations_ptr_->safetyCBFConstraint(current_state, other_pos);
        piecewise_mpc_qp_generator_ptr_->addLinearConstraintForPiecewise(linear_constraint);
    }

    template <typename T, unsigned int DIM>
    void PiecewiseBezierMPCCBFQPGenerator<T, DIM>::addFovLBConstraint(const State &current_state,
                                                                      const Vector &other_pos) {
        LinearConstraint linear_constraint = piecewise_mpc_cbf_operations_ptr_->fovLBConstraint(current_state, other_pos);
        piecewise_mpc_qp_generator_ptr_->addLinearConstraintForPiecewise(linear_constraint);
    }

    template <typename T, unsigned int DIM>
    void PiecewiseBezierMPCCBFQPGenerator<T, DIM>::addFovRBConstraint(const State &current_state,
                                                                      const Vector &other_pos) {
        LinearConstraint linear_constraint = piecewise_mpc_cbf_operations_ptr_->fovRBConstraint(current_state, other_pos);
        piecewise_mpc_qp_generator_ptr_->addLinearConstraintForPiecewise(linear_constraint);
    }

    template class PiecewiseBezierMPCCBFQPGenerator<double, 3U>;

} // mpc_cbf