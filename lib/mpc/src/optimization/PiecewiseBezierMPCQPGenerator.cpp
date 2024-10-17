//
// Created by lishuo on 8/22/24.
//

#include <mpc/optimization/PiecewiseBezierMPCQPGenerator.h>

namespace mpc {
    template <typename T, unsigned int DIM>
    void PiecewiseBezierMPCQPGenerator<T, DIM>::addPiecewise(
            std::unique_ptr<PiecewiseBezierMPCQPOperations> &&piecewise_operations_ptr) {
        piecewise_operations_ptr_ = std::move(piecewise_operations_ptr);

        const std::vector<std::unique_ptr<BezierQPOperations>> &all_piece_operations_ptrs =
                piecewise_operations_ptr_->piece_operations_ptrs();

        // allocate the variable
        for (size_t piece_index = 0; piece_index < piecewise_operations_ptr_->numPieces(); ++piece_index) {
            size_t num_decision_variables_of_piece = all_piece_operations_ptrs.at(piece_index)->numDecisionVariables();
            // add variable map for the piece
            variable_map_.emplace_back(num_decision_variables_of_piece);
            std::vector<qpcpp::Variable<T>*>& new_piece_variables = variable_map_.back();
            // initiate decision variables one by one to the model instance
            for (std::size_t decision_variable_idx = 0;
                 decision_variable_idx < num_decision_variables_of_piece;
                 ++decision_variable_idx) {
                qpcpp::Variable<T>* variable_ptr = problem_.addVariable();
                new_piece_variables.at(decision_variable_idx) = variable_ptr;
                variables_.push_back(variable_ptr);
            }
        }

        assert(variables_.size() == variable_map_.size()*variable_map_[0].size());
    }

    template <typename T, unsigned int DIM>
    size_t PiecewiseBezierMPCQPGenerator<T, DIM>::numPieces() const {
        return piecewise_operations_ptr_->piece_operations_ptrs().size();
    }

    template <typename T, unsigned int DIM>
    const std::unique_ptr<typename PiecewiseBezierMPCQPGenerator<T, DIM>::PiecewiseBezierMPCQPOperations> &
    PiecewiseBezierMPCQPGenerator<T, DIM>::piecewise_operations_ptr() const {
        return piecewise_operations_ptr_;
    }

    template <typename T, unsigned int DIM>
    qpcpp::Problem<T>& PiecewiseBezierMPCQPGenerator<T, DIM>::problem() {
        return problem_;
    }

    template <typename T, unsigned int DIM>
    typename PiecewiseBezierMPCQPGenerator<T, DIM>::SingleParameterPiecewiseCurve
    PiecewiseBezierMPCQPGenerator<T, DIM>::generateCurveFromSolution() const {
        SingleParameterPiecewiseCurve piecewise_curve;

        for (std::size_t piece_idx = 0; piece_idx < piecewise_operations_ptr_->piece_operations_ptrs().size();
             ++piece_idx) {
            std::size_t num_decision_variables_of_piece =
                    piecewise_operations_ptr_->piece_operations_ptrs().at(piece_idx)->numDecisionVariables();
            Vector solution(num_decision_variables_of_piece);
            // get the decision variable value.
            for (std::size_t decision_variable_idx = 0;
                 decision_variable_idx < num_decision_variables_of_piece;
                 ++decision_variable_idx) {
                solution(decision_variable_idx) = variable_map_.at(piece_idx).at(decision_variable_idx)
                        ->solution_value();
            }

            std::unique_ptr<SingleParameterCurve> curve =
                    piecewise_operations_ptr_->piece_operations_ptrs().at(piece_idx)->generateCurveFromSolution(
                            solution);

            piecewise_curve.addPiece(std::move(curve));
        }

        return piecewise_curve;
    }

    template <typename T, unsigned int DIM>
    typename PiecewiseBezierMPCQPGenerator<T, DIM>::Vector
    PiecewiseBezierMPCQPGenerator<T, DIM>::getVariablesValue() const {
        size_t num_decision_variables_of_piecewise = piecewise_operations_ptr_->numDecisionVariables();
        Vector solution(num_decision_variables_of_piecewise);
        assert(solution.size() == variables_.size());
        for (size_t decision_variable_idx = 0;
             decision_variable_idx < num_decision_variables_of_piecewise;
             ++decision_variable_idx) {
            solution(decision_variable_idx) = variables_.at(decision_variable_idx)->solution_value();
        }
        return solution;
    }

    template <typename T, unsigned int DIM>
    void PiecewiseBezierMPCQPGenerator<T, DIM>::addPositionErrorPenaltyCost(const State &current_state, const Vector &ref_positions) {
        const CostAddition cost_addition = piecewise_operations_ptr_->positionErrorPenaltyCost(current_state, ref_positions);
        addCostAdditionForPiecewise(cost_addition);
    }

    template <typename T, unsigned int DIM>
    void PiecewiseBezierMPCQPGenerator<T, DIM>::addEvalPositionErrorPenaltyCost(const Vector &ref_positions) {
        Vector h_samples = piecewise_operations_ptr_->h_samples();
        int spd_f = piecewise_operations_ptr_->mpc_tuning().spd_f_;
        T w_pos_err = piecewise_operations_ptr_->mpc_tuning().w_pos_err_;
        for (size_t k = h_samples.size() - spd_f; k < h_samples.size(); ++k) {
            T h_sample = h_samples[k];
            const Vector ref_pos = ref_positions.segment(DIM*k, DIM);
            const PieceIndexAndParameter& piece_idx_and_parameter = piecewise_operations_ptr_->getPieceIndexAndParameter(h_sample);
            size_t piece_idx = piece_idx_and_parameter.piece_idx();
            CostAddition cost_addition = piecewise_operations_ptr_->piece_operations_ptrs().at(piece_idx)
                    ->evalCost(piece_idx_and_parameter.parameter(), 0, ref_pos, w_pos_err);
            addCostAdditionForPiece(piece_idx, cost_addition);
        }
    }

    template <typename T, unsigned int DIM>
    void PiecewiseBezierMPCQPGenerator<T, DIM>::addControlEffortPenaltyCost() {
        const CostAddition cost_addition = piecewise_operations_ptr_->controlEffortPenaltyCost();
        addCostAdditionForPiecewise(cost_addition);
    }

    template <typename T, unsigned int DIM>
    void
    PiecewiseBezierMPCQPGenerator<T, DIM>::addIntegratedSquaredDerivativeCost(uint64_t derivative_degree, T lambda) {
        for (std::size_t piece_idx = 0; piece_idx < piecewise_operations_ptr_->piece_operations_ptrs().size(); ++piece_idx) {
            const CostAddition &cost_addition = piecewise_operations_ptr_->piece_operations_ptrs().at(piece_idx)
                    ->integratedSquaredDerivativeCost(derivative_degree, lambda);
            addCostAdditionForPiece(piece_idx, cost_addition);
        }
    }

    template <typename T, unsigned int DIM>
    void PiecewiseBezierMPCQPGenerator<T, DIM>::addEvalConstraint(T parameter, uint64_t derivative_degree,
                                                                  const VectorDIM &target) {
        const PieceIndexAndParameter& piece_idx_and_parameter =
                piecewise_operations_ptr_->getPieceIndexAndParameter(parameter);
        const std::vector<LinearConstraint> linear_constraints =
                piecewise_operations_ptr_->piece_operations_ptrs().at(piece_idx_and_parameter.piece_idx())
                        ->evalConstraint(piece_idx_and_parameter.parameter(), derivative_degree, target);
        for (const LinearConstraint& linear_constraint : linear_constraints) {
            addLinearConstraintForPiece(piece_idx_and_parameter.piece_idx(), linear_constraint);
        }
    }

    template <typename T, unsigned int DIM>
    void PiecewiseBezierMPCQPGenerator<T, DIM>::addContinuityConstraint(std::size_t piece_idx,
                                                                        uint64_t derivative_degree) {
        using Row = math::Row<T>;
        // piece_idx out_of_range error
        if (piece_idx + 1 >= piecewise_operations_ptr_->cumulative_max_parameters().size()) {
            throw std::runtime_error("PiecewiseBezierMPCQPGenerator::addContinuityConstraint:"
                                     " piece index out of range. piece_idx: " + std::to_string(piece_idx) +
                                     ", piece_idx + 1: " + std::to_string(piece_idx + 1) +
                                     ", number of pieces: " + std::to_string(piecewise_operations_ptr_->cumulative_max_parameters().size()));
        }

        // get the first piece eval basis (still need to multiply with control point of first piece)
        for (unsigned int dimension = 0; dimension < DIM; ++dimension) {
            const Row& first_piece_row =
                    piecewise_operations_ptr_->piece_operations_ptrs().at(piece_idx)->evalBasisRow(
                            /*dimension=*/dimension,
                            /*parameter=*/piecewise_operations_ptr_->piece_operations_ptrs().at(piece_idx)->max_parameter(),
                            /*derivative_degree=*/derivative_degree);
            const Row& second_piece_row =
                    piecewise_operations_ptr_->piece_operations_ptrs().at(piece_idx+1)->evalBasisRow(
                            /*dimension=*/dimension,
                            /*parameter=*/0,
                            /*derivative_degree=*/derivative_degree);

            // add the linear constraint
            qpcpp::LinearConstraint<T>* constraint_ptr = problem_.addLinearConstraint(/*min=*/0, /*max=*/0);

            for (std::size_t variable_idx = 0;
                 variable_idx < variable_map_.at(piece_idx).size();
                 ++variable_idx) {
                constraint_ptr->setCoefficient(variable_map_.at(piece_idx).at(variable_idx),
                                               first_piece_row(variable_idx));
            }

            for (std::size_t variable_idx = 0;
                 variable_idx < variable_map_.at(piece_idx+1).size();
                 ++variable_idx) {
                constraint_ptr->setCoefficient(variable_map_.at(piece_idx+1).at(variable_idx),
                                               -second_piece_row(variable_idx));
            }
        }
    }

    template <typename T, unsigned int DIM>
    void
    PiecewiseBezierMPCQPGenerator<T, DIM>::addHyperplaneConstraintForPiece(std::size_t piece_idx,
                                                                           const Hyperplane& hyperplane,
                                                                           T epsilon) {
        if (piece_idx >= piecewise_operations_ptr_->cumulative_max_parameters().size()) {
            throw std::invalid_argument("PiecewiseBezierMPCQPGenerator::"
                                        "addHyperplaneConstraintForPiece: piece index out of range. "
                                        "piece_idx: " + std::to_string(piece_idx) +
                                        ", number of pieces: " + std::to_string(piecewise_operations_ptr_->cumulative_max_parameters().size()));
        }

        const std::vector<LinearConstraint> linear_constraints =
                piecewise_operations_ptr_->piece_operations_ptrs().at(piece_idx)->hyperplaneConstraintAll(hyperplane, epsilon);

        for (const LinearConstraint& linear_constraint : linear_constraints) {
            addLinearConstraintForPiece(piece_idx, linear_constraint);
        }

    }

    template <typename T, unsigned int DIM>
    void PiecewiseBezierMPCQPGenerator<T, DIM>::addBoundingBoxConstraintAll(
            const AlignedBox &bounding_box, uint64_t derivative_degree) {
        for (std::size_t piece_idx = 0; piece_idx < piecewise_operations_ptr_->piece_operations_ptrs().size();
             ++piece_idx) {
            const std::vector<LinearConstraint>& piece_constraints =
                    piecewise_operations_ptr_->piece_operations_ptrs().at(piece_idx)->boundingBoxConstraintAll(bounding_box, derivative_degree);

            for (const LinearConstraint& piece_constraint : piece_constraints) {
                addLinearConstraintForPiece(piece_idx, piece_constraint);
            }
        }
    }

    template <typename T, unsigned int DIM>
    void
    PiecewiseBezierMPCQPGenerator<T, DIM>::addHyperplaneConstraintAt(
            T parameter, const Hyperplane& hyperplane, T epsilon) {
        const PieceIndexAndParameter& piece_idx_and_parameter =
                piecewise_operations_ptr_->getPieceIndexAndParameter(parameter);


        const std::vector<LinearConstraint>& linear_constraints =
                piecewise_operations_ptr_->piece_operations_ptrs().at(piece_idx_and_parameter.piece_idx())
                ->hyperplaneConstraintAt(piece_idx_and_parameter.parameter(),
                                         hyperplane, epsilon);

        for (const LinearConstraint& linear_constraint : linear_constraints) {
            addLinearConstraintForPiece(piece_idx_and_parameter.piece_idx(), linear_constraint);
        }

    }

    template <typename T, unsigned int DIM>
    void PiecewiseBezierMPCQPGenerator<T, DIM>::addCostAdditionForPiecewise(const CostAddition &cost_addition) {
        constexpr T epsilon = std::numeric_limits<T>::epsilon() * T(100.0);

        size_t num_decision_variables = variables_.size();

        // number decision variable number error
        if (num_decision_variables != cost_addition.linear_term().rows() ||
            num_decision_variables != cost_addition.quadratic_term().rows() ||
            num_decision_variables != cost_addition.quadratic_term().cols()) {
            throw std::runtime_error("PiecewiseBezierMPCQPGenerator::addCostAdditionForPiecewise:"
                                     " number of decision variables of the piecewise does not match the "
                                     "CostAddition structure");
        }

        // add the constant term
        if (cost_addition.constant() != 0) {
            problem_.cost_function()->add_constant(cost_addition.constant());
        }

        for (size_t i = 0; i < num_decision_variables; ++i) {
            const qpcpp::Variable<T>* var1_ptr = variables_.at(i);
            // add the linear term
            if (!math::isApproximatelyEqual<T>(cost_addition.linear_term()(i), T(0.0), epsilon)) {
                problem_.cost_function()->addLinearTerm(var1_ptr, cost_addition.linear_term()(i));
            }
            // add quadratic term
            for (std::size_t j = 0; j < num_decision_variables; ++j) {
                const qpcpp::Variable<T>* var2_ptr = variables_.at(j);
                // add quadratic term
                if (!math::isApproximatelyEqual<T>(cost_addition.quadratic_term()(i,j), T(0.0), epsilon)) {
                    problem_.cost_function()->addQuadraticTerm(var1_ptr, var2_ptr, cost_addition.quadratic_term()(i,j));
                }
            }
        }
    }

    template <typename T, unsigned int DIM>
    void PiecewiseBezierMPCQPGenerator<T, DIM>::addLinearConstraintForPiecewise(
            const LinearConstraint &linear_constraint) {
        size_t num_decision_variables = variables_.size();

        // linear constraint structure error.
        if (num_decision_variables != linear_constraint.coefficients().cols()) {
            throw std::runtime_error("PiecewiseBezierMPCQPGenerator::"
                                     "addLinearConstraintForPiecewise:"
                                     " number of decision variables of the piecewise does not match the "
                                     "LinearConstraint structure");
        }

        // initialize linear constraint with lower bound and upper bound
        qpcpp::LinearConstraint<T>* qpcpp_linear_constraint = problem_.addLinearConstraint(
                linear_constraint.lower_bound(),
                linear_constraint.upper_bound());
        // set coefficients for this linear constraint
        for (std::size_t decision_variable_idx = 0;
             decision_variable_idx < num_decision_variables;
             ++decision_variable_idx) {
            const qpcpp::Variable<T>* var_ptr = variables_.at(decision_variable_idx);
            qpcpp_linear_constraint->setCoefficient(var_ptr, linear_constraint.coefficients()(decision_variable_idx));
        }
    }

    template <typename T, unsigned int DIM>
    void PiecewiseBezierMPCQPGenerator<T, DIM>::addCostAdditionForPiece(std::size_t piece_idx,
                                                                        const CostAddition &cost_addition) {
        constexpr T epsilon = std::numeric_limits<T>::epsilon() * T(100.0);

        std::size_t num_decision_variables = variable_map_.at(piece_idx).size();

        // number decision variable number error
        if (num_decision_variables != cost_addition.linear_term().rows() ||
            num_decision_variables != cost_addition.quadratic_term().rows() ||
            num_decision_variables != cost_addition.quadratic_term().cols()) {
            throw std::runtime_error("PiecewiseBezierMPCQPGenerator::addCostAdditionForPiece:"
                                     " number of decision variables of the piece does not match the "
                                     "CostAddition structure");
        }
        // piece_idx out_of_range error
        if (piece_idx >= piecewise_operations_ptr_->piece_operations_ptrs().size()) {
            throw std::out_of_range("PiecewiseBezierMPCQPGenerator::addCostAdditionForPiece:"
                                    " piece index out of range. piece_idx: " + std::to_string(piece_idx) +
                                    ", num pieces: " + std::to_string(piecewise_operations_ptr_->piece_operations_ptrs().size()));
        }

        // add the constant term
        if (cost_addition.constant() != 0) {
            problem_.cost_function()->add_constant(cost_addition.constant());
        }

        for (std::size_t i = 0; i < num_decision_variables; ++i) {
            const qpcpp::Variable<T>* var1_ptr = variable_map_.at(piece_idx).at(i);
            // add linear term
            if (!math::isApproximatelyEqual<T>(cost_addition.linear_term()(i), T(0.0), epsilon)) {
                problem_.cost_function()->addLinearTerm(var1_ptr, cost_addition.linear_term()(i));
            }
            // add quadratic term
            for (std::size_t j = 0; j < num_decision_variables; ++j) {
                const qpcpp::Variable<T>* var2_ptr = variable_map_.at(piece_idx).at(j);
                // add quadratic term
                if (!math::isApproximatelyEqual<T>(cost_addition.quadratic_term()(i,j), T(0.0), epsilon)) {
                    problem_.cost_function()->addQuadraticTerm(var1_ptr, var2_ptr, cost_addition.quadratic_term()(i,j));
                }
            }
        }
    }

    template <typename T, unsigned int DIM>
    void PiecewiseBezierMPCQPGenerator<T, DIM>::addLinearConstraintForPiece(std::size_t piece_idx,
                                                                            const LinearConstraint &linear_constraint) {
        std::size_t num_decision_variables = variable_map_.at(piece_idx).size();

        // linear constraint structure error.
        if (num_decision_variables != linear_constraint.coefficients().cols()) {
            throw std::runtime_error("PiecewiseBezierMPCQPGenerator::"
                                     "addLinearConstraintForPiece:"
                                     " number of decision variables of the piece does not match the "
                                     "LinearConstraint structure");
        }

        // initialize linear constraint with lower bound and upper bound
        qpcpp::LinearConstraint<T>* qpcpp_linear_constraint = problem_.addLinearConstraint(
                linear_constraint.lower_bound(),
                linear_constraint.upper_bound());
        // set coefficients for this linear constraint
        for (std::size_t decision_variable_idx = 0;
             decision_variable_idx < num_decision_variables;
             ++decision_variable_idx) {
            const qpcpp::Variable<T>* var_ptr = variable_map_.at(piece_idx).at(decision_variable_idx);
            qpcpp_linear_constraint->setCoefficient(var_ptr, linear_constraint.coefficients()(decision_variable_idx));
        }
    }

    template class PiecewiseBezierMPCQPGenerator<double, 2U>;
    template class PiecewiseBezierMPCQPGenerator<float, 2U>;
    template class PiecewiseBezierMPCQPGenerator<double, 3U>;
    template class PiecewiseBezierMPCQPGenerator<float, 3U>;
} // mpc