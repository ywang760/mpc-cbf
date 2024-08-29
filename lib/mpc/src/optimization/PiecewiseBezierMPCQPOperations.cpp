//
// Created by lishuo on 8/19/24.
//

#include <mpc/optimization/PiecewiseBezierMPCQPOperations.h>

namespace mpc {
    template <typename T, unsigned int DIM>
    PiecewiseBezierMPCQPOperation<T, DIM>::PiecewiseBezierMPCQPOperation(Params &p,
                                                                         std::shared_ptr<DoubleIntegrator> model_ptr)
    {
        // init the piecewise bezier MPC QP Generator
        for (size_t piece_index = 0; piece_index < p.piecewise_bezier_params.num_pieces_; ++piece_index) {
            typename BezierQPOperations::Params bezier_p({p.piecewise_bezier_params.num_control_points_,
                                                          p.piecewise_bezier_params.piece_max_parameter_});
            std::unique_ptr<BezierQPOperations> bezier_qp_operation = std::make_unique<BezierQPOperations>(bezier_p);
            addPiece(std::move(bezier_qp_operation));
        }

        h_ = p.mpc_params.h_;
        k_hor_ = p.mpc_params.k_hor_;
        mpc_tuning_ = p.mpc_params.tuning_;
        // compute the A0 and Lambda for the horizon
        model_ptr_ = model_ptr;
        //TODO assert the model 1. ts is h 2. dim is DIM

        // model predict
        A0_ = model_ptr_->get_A0(k_hor_); // A0.pos_: [3K, 6], A0.vel_: [3K, 6],
        Lambda_ = model_ptr_->get_lambda(k_hor_); // Lambda_.pos_: [3K, 3K], Lambda_.vel_: [3K, 3K]
        // control sequence U control point coefficient
        Vector h_samples = Vector::LinSpaced(k_hor_, 0, (k_hor_-1)*h_);
        U_basis_ = evalSamplingBasisMatrix(h_samples, 2); // [3K, num_piece*dim*num_control_pts] here since the control is acc, derivative degree is 2

    }

    template <typename T, unsigned int DIM>
    typename PiecewiseBezierMPCQPOperation<T, DIM>::Matrix
    PiecewiseBezierMPCQPOperation<T, DIM>::evalSamplingBasisMatrix(Vector &h_samples, uint64_t derivative_degree) {
        int hor = h_samples.size();
        Matrix result = Matrix::Zero(DIM*hor, numDecisionVariables());
        for (size_t k = 0; k < h_samples.size(); ++k) {
            T h_sample = h_samples[k];
            const PieceIndexAndParameter& piece_idx_and_parameter = getPieceIndexAndParameter(h_sample);
            size_t piece_idx = piece_idx_and_parameter.piece_idx();
            for (size_t d = 0; d < DIM; ++d) {
                Row piece_dim_basis = piece_operations_ptrs_.at(piece_idx)->evalBasisRow(
                        d, piece_idx_and_parameter.parameter(), derivative_degree);
                size_t piece_num_decision_vars = piece_dim_basis.size();
                Eigen::Ref<Matrix> result_row = result.block(k*DIM+d, 0, 1, numDecisionVariables());
                result_row.block(0, piece_idx*piece_num_decision_vars, 1, piece_num_decision_vars) = piece_dim_basis;
            }
        }
        return result;
    }

    template <typename T, unsigned int DIM>
    typename PiecewiseBezierMPCQPOperation<T, DIM>::CostAddition
    PiecewiseBezierMPCQPOperation<T, DIM>::positionErrorPenaltyCost(const State &current_state, const VectorDIM &target) {
        Matrix quadratic_term(numDecisionVariables(), numDecisionVariables());
        quadratic_term.setZero();
        Vector linear_term(numDecisionVariables());
        linear_term.setZero();

        // weight matrix
        Matrix Q_pe = Matrix::Zero(DIM*k_hor_, DIM*k_hor_);
        Eigen::Ref<Matrix> block_Q_pe = Q_pe.block(
                DIM * (k_hor_ - mpc_tuning_.spd_f_), DIM * (k_hor_ - mpc_tuning_.spd_f_),
                DIM * mpc_tuning_.spd_f_, DIM * mpc_tuning_.spd_f_);

        block_Q_pe = mpc_tuning_.w_pos_err_ * Matrix::Identity(DIM * mpc_tuning_.spd_f_, DIM * mpc_tuning_.spd_f_);


        // compute quadratic term
        Matrix Phi_pred = Lambda_.pos_ * U_basis_; // [3K, num_piece*dim*num_control_pts]
        quadratic_term = Phi_pred.transpose() * Q_pe * Phi_pred;
        // compute the linear term
        Vector x0 = Vector::Zero(2*DIM);
        x0 << current_state.pos_, current_state.vel_; // [6,1]
        Matrix linear_term_coef = (A0_.pos_ * x0).transpose() * Q_pe;
        linear_term_coef += -1.0 * target.replicate(k_hor_, 1).transpose() * Q_pe;
        linear_term = (linear_term_coef * Phi_pred).transpose();

        return CostAddition(quadratic_term, linear_term, 0);
    }

    template <typename T, unsigned int DIM>
    typename PiecewiseBezierMPCQPOperation<T, DIM>::CostAddition
    PiecewiseBezierMPCQPOperation<T, DIM>::controlEffortPenaltyCost() {
        Matrix quadratic_term(numDecisionVariables(), numDecisionVariables());
        quadratic_term.setZero();
        Vector linear_term(numDecisionVariables());
        linear_term.setZero();

        // weight matrix
        Matrix Q_ue = mpc_tuning_.w_u_eff_ * Matrix::Identity(DIM*k_hor_, DIM*k_hor_); // [3K, 3K]
        // compute the control effort penalty cost
        quadratic_term = U_basis_.transpose() * Q_ue * U_basis_; // [num_pieces*dim*num_control_pts, num_pieces*dim*num_control_pts]

        return CostAddition(quadratic_term, linear_term, 0);
    }

    template <typename T, unsigned int DIM>
    std::size_t PiecewiseBezierMPCQPOperation<T, DIM>::numDecisionVariables() const {
        size_t num_pieces = piece_operations_ptrs_.size();
        return num_pieces * piece_operations_ptrs_[0]->numDecisionVariables();
    }

    template <typename T, unsigned int DIM>
    size_t PiecewiseBezierMPCQPOperation<T, DIM>::numPieces() const {
        return cumulative_max_parameters_.size();
    }

    template <typename T, unsigned int DIM>
    T PiecewiseBezierMPCQPOperation<T, DIM>::max_parameter() const {
        if (cumulative_max_parameters_.empty()) {
            throw std::runtime_error("there is no piece in the curve");
        }
        return cumulative_max_parameters_.back();
    }

    template <typename T, unsigned int DIM>
    const std::vector<std::unique_ptr<typename PiecewiseBezierMPCQPOperation<T, DIM>::BezierQPOperations>> &
    PiecewiseBezierMPCQPOperation<T, DIM>::piece_operations_ptrs() const {
        return piece_operations_ptrs_;
    }

    template <typename T, unsigned int DIM>
    const std::vector<T> &PiecewiseBezierMPCQPOperation<T, DIM>::cumulative_max_parameters() const {
        return cumulative_max_parameters_;
    }

    template <typename T, unsigned int DIM>
    void PiecewiseBezierMPCQPOperation<T, DIM>::addPiece(
            std::unique_ptr<BezierQPOperations> &&piece_operations_ptr) {
        // augment the cumulative max parameters
        if (cumulative_max_parameters_.empty()) {
            cumulative_max_parameters_.push_back(
                    piece_operations_ptr->max_parameter());
        } else {
            // add the duration of new piece to the cumulative max parameters
            cumulative_max_parameters_.push_back(
                    cumulative_max_parameters_.back() +
                    piece_operations_ptr->max_parameter());
        }
        // push in the piece operations pointer
        piece_operations_ptrs_.emplace_back(std::move(piece_operations_ptr));

    }


    template <typename T, unsigned int DIM>
    PiecewiseBezierMPCQPOperation<T, DIM>::PieceIndexAndParameter::
    PieceIndexAndParameter(std::size_t piece_idx, T parameter)
            : piece_idx_(piece_idx), parameter_(parameter) {}

    template <typename T, unsigned int DIM>
    std::size_t
    PiecewiseBezierMPCQPOperation<T, DIM>::PieceIndexAndParameter::piece_idx() const {
        return piece_idx_;
    }

    template <typename T, unsigned int DIM>
    T
    PiecewiseBezierMPCQPOperation<T, DIM>::PieceIndexAndParameter::parameter() const {
        return parameter_;
    }

    template <typename T, unsigned int DIM>
    typename PiecewiseBezierMPCQPOperation<T, DIM>::PieceIndexAndParameter
    PiecewiseBezierMPCQPOperation<T, DIM>::getPieceIndexAndParameter(T parameter) const {
        if (parameter < 0 || parameter > cumulative_max_parameters_.back()) {
            throw std::runtime_error("SingleParameterPiecewiseCurveQPGenerator::pieceIndexAndParameter: "
                                     "parameter is out of range [0,max parameter]. max_parameter: " +
                                     std::to_string(cumulative_max_parameters_.back()) +
                                     ", parameter: " + std::to_string(parameter));
        }

        typename std::vector<T>::const_iterator cumulative_parameter_iterator =
                std::lower_bound(cumulative_max_parameters_.begin(),
                                 cumulative_max_parameters_.end(), parameter);

        std::size_t piece_idx = std::distance(cumulative_max_parameters_.begin(),
                                              cumulative_parameter_iterator);

        if (piece_idx >= piece_operations_ptrs_.size()) {
            throw std::runtime_error("SingleParameterPiecewiseCurveQPGenerator::pieceIndexAndParameter: "
                                     "piece_idx is out of range [0, num_pieces). num_pieces: " +
                                     std::to_string(piece_operations_ptrs_.size()) +
                                     ", piece_idx: " + std::to_string(piece_idx) +
                                     ". This should not have happened, please report this bug.");
        }

        const T piece_max_parameter =
                piece_operations_ptrs_[piece_idx]->max_parameter();

        if (piece_idx == 0) {
            return PieceIndexAndParameter(
                    piece_idx, std::clamp(parameter, T(0.0), piece_max_parameter));
        } else {
            return PieceIndexAndParameter(
                    piece_idx,
                    std::clamp(parameter - cumulative_max_parameters_.at(piece_idx - 1),
                               T(0.0), piece_max_parameter));
        }
    }

    template class PiecewiseBezierMPCQPOperation<double, 3U>;
    template class PiecewiseBezierMPCQPOperation<float, 3U>;
    template class PiecewiseBezierMPCQPOperation<double, 2U>;
    template class PiecewiseBezierMPCQPOperation<float, 2U>;
} // mpc