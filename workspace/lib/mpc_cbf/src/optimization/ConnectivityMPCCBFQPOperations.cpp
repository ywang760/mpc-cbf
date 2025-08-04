//
// Created by yutong on 8/4/25.
//

#include <mpc_cbf/optimization/ConnectivityMPCCBFQPOperations.h>

namespace mpc_cbf {
    template <typename T, unsigned int DIM>
    ConnectivityMPCCBFQPOperations<T, DIM>::ConnectivityMPCCBFQPOperations(Params &p, std::shared_ptr<DoubleIntegrator> model_ptr, std::shared_ptr<ConnectivityCBF> connectivity_cbf_ptr)
            : MPCCBFQPOperationsBase<T, DIM>(model_ptr), connectivity_cbf_ptr_(connectivity_cbf_ptr) {
                // mpc operations
                typename PiecewiseBezierMPCQPOperations::Params bezier_mpc_p = {p.piecewise_bezier_params, p.mpc_params};
                this->piecewise_mpc_operations_ptr_ = std::make_unique<PiecewiseBezierMPCQPOperations>(bezier_mpc_p, model_ptr);

                // mpc params
                this->h_ = p.mpc_params.h_;
                this->k_hor_ = p.mpc_params.k_hor_;
                this->mpc_tuning_ = p.mpc_params.tuning_;
                // control input predict
                this->U_basis_ = this->piecewise_mpc_operations_ptr_->U_basis(); // [3K, num_piece*dim*num_control_pts]
            }

    template <typename T, unsigned int DIM>
    typename ConnectivityMPCCBFQPOperations<T, DIM>::LinearConstraint
    ConnectivityMPCCBFQPOperations<T, DIM>::safetyCBFConstraint(const State &current_state,
                                                               const Vector &other_pos,
                                                               T slack_value) {
        Vector state(2*DIM);
        state << current_state.pos_, current_state.vel_;
        
        // TODO: Use connectivity_cbf_ptr_ for proper constraint generation
        // For now, implement a simple safety constraint
        Vector a = Vector::Zero(DIM);
        VectorDIM relative_pos = current_state.pos_ - other_pos.segment(0, DIM);
        T distance = relative_pos.norm();
        if (distance > 0.01) { // Avoid division by zero
            a.segment(0, DIM) = relative_pos / distance;
        }
        
        T b = 1.0; // Safety bound (distance >= 1.0)
        Vector A0 = Vector::Zero(DIM*this->k_hor_);
        A0.segment(0, DIM) = a;
        math::Row<T> A_control_pts = -1.0 * A0.transpose() * this->U_basis_;
        
        return LinearConstraint(A_control_pts, std::numeric_limits<T>::lowest(), b + slack_value);
    }

    template <typename T, unsigned int DIM>
    typename ConnectivityMPCCBFQPOperations<T, DIM>::LinearConstraint
    ConnectivityMPCCBFQPOperations<T, DIM>::connectivityConstraint(const State &current_state,
                                                                  const Vector &other_pos,
                                                                  T slack_value) {
        Vector state(2*DIM);
        state << current_state.pos_, current_state.vel_;
        
        // TODO: Use connectivity_cbf_ptr_ for proper constraint generation
        // For now, implement a simple connectivity constraint
        Vector a = Vector::Zero(DIM);
        VectorDIM relative_pos = current_state.pos_ - other_pos.segment(0, DIM);
        T distance = relative_pos.norm();
        if (distance > 0.01) { // Avoid division by zero
            a.segment(0, DIM) = -relative_pos / distance; // Negative for upper bound constraint
        }
        
        T b = 10.0; // Connectivity bound (distance <= 10.0)
        Vector A0 = Vector::Zero(DIM*this->k_hor_);
        A0.segment(0, DIM) = a;
        math::Row<T> A_control_pts = -1.0 * A0.transpose() * this->U_basis_;
        
        return LinearConstraint(A_control_pts, std::numeric_limits<T>::lowest(), b + slack_value);
    }

    template <typename T, unsigned int DIM>
    std::vector<typename ConnectivityMPCCBFQPOperations<T, DIM>::LinearConstraint>
    ConnectivityMPCCBFQPOperations<T, DIM>::predSafetyCBFConstraints(const std::vector<State> &pred_states,
                                                                    const Vector &other_pos) {
        std::vector<LinearConstraint> linear_constraints;
        for (size_t k = 0; k < pred_states.size(); ++k) {
            const State& pred_state = pred_states.at(k);
            Vector state(2*DIM);
            state << pred_state.pos_, pred_state.vel_;
            
            // TODO: Use connectivity_cbf_ptr_ for proper constraint generation
            Vector ak = Vector::Zero(DIM);
            VectorDIM relative_pos = pred_state.pos_ - other_pos.segment(0, DIM);
            T distance = relative_pos.norm();
            if (distance > 0.01) {
                ak.segment(0, DIM) = relative_pos / distance;
            }
            T bk = 1.0; // Safety bound
            
            Vector Ak = Vector::Zero(DIM*this->k_hor_);
            Ak.segment(k*DIM, DIM) = ak;
            math::Row<T> Ak_control_pts = -1.0 * Ak.transpose() * this->U_basis_;
            linear_constraints.push_back(LinearConstraint(Ak_control_pts, std::numeric_limits<T>::lowest(), bk));
        }
        return linear_constraints;
    }

    template <typename T, unsigned int DIM>
    std::vector<typename ConnectivityMPCCBFQPOperations<T, DIM>::LinearConstraint>
    ConnectivityMPCCBFQPOperations<T, DIM>::predConnectivityConstraints(const std::vector<State> &pred_states,
                                                                        const Vector &other_pos) {
        std::vector<LinearConstraint> linear_constraints;
        for (size_t k = 0; k < pred_states.size(); ++k) {
            const State& pred_state = pred_states.at(k);
            Vector state(2*DIM);
            state << pred_state.pos_, pred_state.vel_;
            
            // TODO: Use connectivity_cbf_ptr_ for proper constraint generation
            Vector ak = Vector::Zero(DIM);
            VectorDIM relative_pos = pred_state.pos_ - other_pos.segment(0, DIM);
            T distance = relative_pos.norm();
            if (distance > 0.01) {
                ak.segment(0, DIM) = -relative_pos / distance; // Negative for upper bound
            }
            T bk = 10.0; // Connectivity bound
            
            Vector Ak = Vector::Zero(DIM*this->k_hor_);
            Ak.segment(k*DIM, DIM) = ak;
            math::Row<T> Ak_control_pts = -1.0 * Ak.transpose() * this->U_basis_;
            linear_constraints.push_back(LinearConstraint(Ak_control_pts, std::numeric_limits<T>::lowest(), bk));
        }
        return linear_constraints;
    }

    template <typename T, unsigned int DIM>
    std::unique_ptr<typename ConnectivityMPCCBFQPOperations<T, DIM>::PiecewiseBezierMPCQPOperations>
    ConnectivityMPCCBFQPOperations<T, DIM>::piecewise_mpc_operations_ptr() {
        return std::move(this->piecewise_mpc_operations_ptr_);
    }

    // Explicit template instantiation
    template class ConnectivityMPCCBFQPOperations<double, 3U>;

} // mpc_cbf
