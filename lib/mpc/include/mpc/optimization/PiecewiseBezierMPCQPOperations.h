//
// Created by lishuo on 8/19/24.
//

#ifndef MPC_PIECEWISEBEZIERMPCQPOPERATIONS_H
#define MPC_PIECEWISEBEZIERMPCQPOPERATIONS_H

#include <math/Types.h>
#include <model/DoubleIntegrator.h>
#include <splines/optimization/SingleParameterCurveQPOperations.h>
#include <splines/optimization/BezierQPOperations.h>
#include <splines/optimization/SingleParameterPiecewiseCurveQPGenerator.h>

namespace mpc {
    // params for piecewise bezier
    template <typename T, unsigned int DIM>
    struct PiecewiseBezierParams {
        size_t num_pieces_, num_control_points_;
        T piece_max_parameter_;
    };

    // params for MPC
    template <typename T>
    struct PhysicalLimits {
        math::Vector<T> p_max_, p_min_, a_max_, a_min_;
    };

    template <typename T>
    struct TuningParams {
        T w_pos_err_, w_u_eff_;
        int spd_f_;
    };

    template <typename T>
    struct MPCParams {
        T h_, Ts_;
        int k_hor_;
        TuningParams<T> tuning_;
        PhysicalLimits<T> limits_;
    };

    template <typename T, unsigned int DIM>
    class PiecewiseBezierMPCQPOperation {
    public:
        using QPOperation = splines::QPOperations<T>;
        using BezierQPOperations = splines::BezierQPOperations<T, DIM>;
        using CostAddition = typename QPOperation::CostAddition;
        using LinearConstraint = typename QPOperation::LinearConstraint;
        using DoubleIntegrator = model::DoubleIntegrator<T, DIM>;
        using State = model::State<T, DIM>;
        using StatePropagator = model::StatePropagator<T>;
        using VectorDIM = math::VectorDIM<T, DIM>;
        using Matrix = math::Matrix<T>;
        using Vector = math::Vector<T>;
        using Row = math::Row<T>;


        struct Params {
            PiecewiseBezierParams<T, DIM> &piecewise_bezier_params;
            MPCParams<T> &mpc_params;
        };
        PiecewiseBezierMPCQPOperation(Params &p, std::shared_ptr<DoubleIntegrator> model_ptr);
        ~PiecewiseBezierMPCQPOperation()=default;

        // return the total number of decision variables of the piecewise curve
        std::size_t numDecisionVariables() const;
        // return the number of pieces this piecewise curve has
        size_t numPieces() const;
        // returns the maximum parameter of the piecewise curve. return status is
        // not ok if there is no pieces.
        T max_parameter() const;
        const std::vector<std::unique_ptr<BezierQPOperations>>& piece_operations_ptrs() const;
        const std::vector<T>& cumulative_max_parameters() const;


        Matrix evalSamplingBasisMatrix(Vector &h_samples, uint64_t derivative_degree);

        /// objective
        CostAddition positionErrorPenaltyCost(const State &current_state, const VectorDIM &target);
        CostAddition controlEffortPenaltyCost();

        /// constraints


        // Container class for piece index and parameter pair. return type of
        // getPieceIndexAndParameter
        class PieceIndexAndParameter {
        public:
            PieceIndexAndParameter(std::size_t piece_idx, T parameter);
            std::size_t piece_idx() const;
            T parameter() const;

        private:
            std::size_t piece_idx_;
            T parameter_;
        };
        // returns the piece index and parameter within the piece with the piece
        // index that the `parameter` argument corresponds to.
        //
        // e.g. if there are 4 pieces with max parameters [1.0,2.0,3.0,4.0] in the
        // piecewise curve, getPieceIndexAndParameter(1.2) returns
        // PieceIndexAndParameter(1, 0.2), because parameter 1.2 corresponds to
        // piece with index 1 and parameter 0.2. getPieceIndexAndParameter(1.0)
        // returns PieceIndexAndParameter(1, 0.0) and getPieceIndexAndParameter(0.0)
        // returns PieceIndexAndParameter(0, 0.0).
        PieceIndexAndParameter getPieceIndexAndParameter(T parameter) const;

    private:
        // adds operations for the next piece
        void addPiece(std::unique_ptr<BezierQPOperations> &&piece_operations_ptr);

        // mpc params
        T h_;
        int k_hor_;
        TuningParams<T> mpc_tuning_;
        // model predict
        std::shared_ptr<DoubleIntegrator> model_ptr_;
        StatePropagator A0_;
        StatePropagator Lambda_;
        Matrix U_basis_;

        // contains pointers to operations for pieces
        std::vector<std::unique_ptr<BezierQPOperations>> piece_operations_ptrs_;

        // contains cumulative max parameters of pieces.
        // cumulative_max_parameters_[i] gives the total max parameters of pieces
        // with indices [0, ..., i]
        std::vector<T> cumulative_max_parameters_;

    };

} // mpc

#endif //MPC_PIECEWISEBEZIERMPCQPOPERATIONS_H
