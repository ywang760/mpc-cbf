//
// Created by lishuo on 9/20/24.
//

#include <mpc_cbf/optimization/FovMPCCBFQPOperations.h>

namespace mpc_cbf {
    template <typename T, unsigned int DIM>
    FovMPCCBFQPOperations<T, DIM>::FovMPCCBFQPOperations(
            Params &p,
            std::shared_ptr<DoubleIntegrator> model_ptr,
            std::shared_ptr<FovCBF> fov_cbf_ptr)
            : MPCCBFQPOperationsBase<T, DIM>(model_ptr), fov_cbf_ptr_(fov_cbf_ptr) {
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
    typename FovMPCCBFQPOperations<T, DIM>::LinearConstraint
    FovMPCCBFQPOperations<T, DIM>::safetyCBFConstraint(const State &current_state,
                                                                   const Vector &other_pos,
                                                                   T slack_value) {
        Vector state(2*DIM);
        state << current_state.pos_, current_state.vel_;
        // cbf constraint A,b
        Vector a = fov_cbf_ptr_->getSafetyConstraints(state, other_pos); // [3, 1]
        T b = fov_cbf_ptr_->getSafetyBound(state, other_pos);
        Vector A0 = Vector::Zero(DIM*this->k_hor_); // [3K, 1]
        A0.segment(0, DIM) = a;
        Row A_control_pts = -1.0 * A0.transpose() * this->U_basis_; // [1, num_piece*dim*num_control_pts]
        return LinearConstraint(A_control_pts, std::numeric_limits<T>::lowest(), b+slack_value);
    }

    template <typename T, unsigned int DIM>
    std::vector<typename FovMPCCBFQPOperations<T, DIM>::LinearConstraint>
    FovMPCCBFQPOperations<T, DIM>::fovLBConstraint(const State &current_state,
                                                               const Vector &other_pos,
                                                               T slack_value) {
        Vector state(2*DIM);
        state << current_state.pos_, current_state.vel_;
        // cbf constraint A,b
        Vector a = fov_cbf_ptr_->getLBConstraints(state, other_pos); // [3, 1]
        T b = fov_cbf_ptr_->getLBBound(state, other_pos);

        Vector A0 = Vector::Zero(DIM*this->k_hor_); // [3K, 1]
        A0.segment(0, DIM) = a;
        Row A_control_pts = -1.0 * A0.transpose() * this->U_basis_; // [1, num_piece*dim*num_control_pts]

        std::vector<LinearConstraint> linear_constraints;
        linear_constraints.push_back(LinearConstraint(A_control_pts, std::numeric_limits<T>::lowest(), b+slack_value));

        Vector A1 = Vector::Zero(DIM*this->k_hor_); // [3K, 1]
        A1.segment(DIM, DIM) = a;
        Row A1_control_pts = -1.0 * A1.transpose() * this->U_basis_; // [1, num_piece*dim*num_control_pts]
        return linear_constraints;
    }

    template <typename T, unsigned int DIM>
    std::vector<typename FovMPCCBFQPOperations<T, DIM>::LinearConstraint>
    FovMPCCBFQPOperations<T, DIM>::fovRBConstraint(const State &current_state,
                                                               const Vector &other_pos,
                                                               T slack_value) {
        Vector state(2*DIM);
        state << current_state.pos_, current_state.vel_;
        // cbf constraint A,b
        Vector a = fov_cbf_ptr_->getRBConstraints(state, other_pos); // [3, 1]
        T b = fov_cbf_ptr_->getRBBound(state, other_pos);

        Vector A0 = Vector::Zero(DIM*this->k_hor_); // [3K, 1]
        A0.segment(0, DIM) = a;
        Row A_control_pts = -1.0 * A0.transpose() * this->U_basis_; // [1, num_piece*dim*num_control_pts]
        std::vector<LinearConstraint> linear_constraints;

        linear_constraints.push_back(LinearConstraint(A_control_pts, std::numeric_limits<T>::lowest(), b+slack_value));

        Vector A1 = Vector::Zero(DIM*this->k_hor_); // [3K, 1]
        A1.segment(DIM, DIM) = a;
        Row A1_control_pts = -1.0 * A1.transpose() * this->U_basis_; // [1, num_piece*dim*num_control_pts]
        return linear_constraints;
    }

    template <typename T, unsigned int DIM>
    typename FovMPCCBFQPOperations<T, DIM>::LinearConstraint
    FovMPCCBFQPOperations<T, DIM>::rangeCBFConstraint(const State &current_state,
                                                                   const Vector &other_pos,
                                                                   T slack_value) {
        Vector state(2*DIM);
        state << current_state.pos_, current_state.vel_;
        // cbf constraint A,b
        Vector a = fov_cbf_ptr_->getRangeConstraints(state, other_pos); // [3, 1]
        T b = fov_cbf_ptr_->getRangeBound(state, other_pos);
        Vector A0 = Vector::Zero(DIM*this->k_hor_); // [3K, 1]
        A0.segment(0, DIM) = a;
        Row A_control_pts = -1.0 * A0.transpose() * this->U_basis_; // [1, num_piece*dim*num_control_pts]
        return LinearConstraint(A_control_pts, std::numeric_limits<T>::lowest(), b+slack_value);
    }

    template <typename T, unsigned int DIM>
    std::vector<typename FovMPCCBFQPOperations<T, DIM>::LinearConstraint>
    FovMPCCBFQPOperations<T, DIM>::predSafetyCBFConstraints(const std::vector<State> &pred_states,
                                                                        const Vector &other_pos) {
        std::vector<LinearConstraint> linear_constraints;
        for (size_t k = 0; k < pred_states.size(); ++k) {
            const State& pred_state = pred_states.at(k);
            Vector state(2*DIM);
            state << pred_state.pos_, pred_state.vel_;
            Vector ak = fov_cbf_ptr_->getSafetyConstraints(state, other_pos);
            T bk = fov_cbf_ptr_->getSafetyBound(state, other_pos);

            Vector Ak = Vector::Zero(DIM*this->k_hor_); // [3K, 1]
            Ak.segment(k*DIM, DIM) = ak;
            Row Ak_control_pts = -1.0 * Ak.transpose() * this->U_basis_; // [1, num_piece*dim*num_control_pts]
            linear_constraints.push_back(LinearConstraint(Ak_control_pts, std::numeric_limits<T>::lowest(), bk));
        }
        return linear_constraints;
    }

    template <typename T, unsigned int DIM>
    std::vector<typename FovMPCCBFQPOperations<T, DIM>::LinearConstraint>
    FovMPCCBFQPOperations<T, DIM>::predFovLBConstraints(const std::vector<State> &pred_states,
                                                                    const Vector &other_pos) {
        std::vector<LinearConstraint> linear_constraints;
        for (size_t k = 0; k < pred_states.size(); ++k) {
            const State& pred_state = pred_states.at(k);
            Vector state(2*DIM);
            state << pred_state.pos_, pred_state.vel_;
            Vector ak = fov_cbf_ptr_->getLBConstraints(state, other_pos);
            T bk = fov_cbf_ptr_->getLBBound(state, other_pos);

            Vector Ak = Vector::Zero(DIM*this->k_hor_); // [3K, 1]
            Ak.segment(k*DIM, DIM) = ak;
            Row Ak_control_pts = -1.0 * Ak.transpose() * this->U_basis_; // [1, num_piece*dim*num_control_pts]
            linear_constraints.push_back(LinearConstraint(Ak_control_pts, std::numeric_limits<T>::lowest(), bk));
        }
        return linear_constraints;
    }

    template <typename T, unsigned int DIM>
    std::vector<typename FovMPCCBFQPOperations<T, DIM>::LinearConstraint>
    FovMPCCBFQPOperations<T, DIM>::predFovRBConstraints(const std::vector<State> &pred_states,
                                                                    const Vector &other_pos) {
        std::vector<LinearConstraint> linear_constraints;
        for (size_t k = 0; k < pred_states.size(); ++k) {
            const State& pred_state = pred_states.at(k);
            Vector state(2*DIM);
            state << pred_state.pos_, pred_state.vel_;
            // cbf constraint A,b
            Vector ak = fov_cbf_ptr_->getRBConstraints(state, other_pos);
            T bk = fov_cbf_ptr_->getRBBound(state, other_pos);

            Vector Ak = Vector::Zero(DIM*this->k_hor_); // [3K, 1]
            Ak.segment(k*DIM, DIM) = ak;
            Row Ak_control_pts = -1.0 * Ak.transpose() * this->U_basis_; // [1, num_piece*dim*num_control_pts]
            linear_constraints.push_back(LinearConstraint(Ak_control_pts, std::numeric_limits<T>::lowest(), bk));
        }
        return linear_constraints;
    }

    template <typename T, unsigned int DIM>
    std::vector<typename FovMPCCBFQPOperations<T, DIM>::LinearConstraint>
    FovMPCCBFQPOperations<T, DIM>::predRangeCBFConstraints(const std::vector<State> &pred_states,
                                                                        const Vector &other_pos) {
        std::vector<LinearConstraint> linear_constraints;
        for (size_t k = 0; k < pred_states.size(); ++k) {
            const State& pred_state = pred_states.at(k);
            Vector state(2*DIM);
            state << pred_state.pos_, pred_state.vel_;
            Vector ak = fov_cbf_ptr_->getRangeConstraints(state, other_pos);
            T bk = fov_cbf_ptr_->getRangeBound(state, other_pos);

            Vector Ak = Vector::Zero(DIM*this->k_hor_); // [3K, 1]
            Ak.segment(k*DIM, DIM) = ak;
            Row Ak_control_pts = -1.0 * Ak.transpose() * this->U_basis_; // [1, num_piece*dim*num_control_pts]
            linear_constraints.push_back(LinearConstraint(Ak_control_pts, std::numeric_limits<T>::lowest(), bk));
        }
        return linear_constraints;
    }

    template class FovMPCCBFQPOperations<double, 3U>;

} // mpc_cbf