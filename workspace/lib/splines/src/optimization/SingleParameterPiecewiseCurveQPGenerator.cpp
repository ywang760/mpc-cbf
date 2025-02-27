//
// Created by lishuo on 2/13/24.
//

#include <splines/optimization/SingleParameterPiecewiseCurveQPGenerator.h>

namespace splines {
    template <typename T, unsigned int DIM>
    void SingleParameterPiecewiseCurveQPGenerator<T, DIM>::addPiece(
            std::unique_ptr<SingleParameterCurveQPOperations> &&piece_operations_ptr) {
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
        // augment the variable map
        std::size_t num_decision_variables_of_new_piece =
                piece_operations_ptrs_.back()->numDecisionVariables();

        // add variable map for the new piece
        variable_map_.emplace_back(num_decision_variables_of_new_piece);
        std::vector<qpcpp::Variable<T>*>& new_piece_variables = variable_map_.back();

        // initiate decision variables one by one to the model instance
        for (std::size_t decision_variable_idx = 0;
        decision_variable_idx < num_decision_variables_of_new_piece;
        ++decision_variable_idx) {
            qpcpp::Variable<T>* variable_ptr = problem_.addVariable();
            new_piece_variables.at(decision_variable_idx) = variable_ptr;
        }
    }

    template <typename T, unsigned int DIM>
    std::size_t SingleParameterPiecewiseCurveQPGenerator<T, DIM>::numPieces() const {
        return cumulative_max_parameters_.size();
    }

    template <typename T, unsigned int DIM>
    T SingleParameterPiecewiseCurveQPGenerator<T, DIM>::max_parameter() const {
        if (cumulative_max_parameters_.empty()) {
            throw std::runtime_error("there is no piece in the curve");
        }
        return cumulative_max_parameters_.back();
    }

    template <typename T, unsigned int DIM>
    typename SingleParameterPiecewiseCurveQPGenerator<T, DIM>::SingleParameterPiecewiseCurve
    SingleParameterPiecewiseCurveQPGenerator<T, DIM>::generateCurveFromSolution() const {
        using Vector = typename SingleParameterCurveQPOperations::Vector;
        using SingleParameterCurve =
                typename SingleParameterCurveQPOperations::SingleParameterCurve;

        SingleParameterPiecewiseCurve piecewise_curve;

        for (std::size_t piece_idx = 0; piece_idx < piece_operations_ptrs_.size();
             ++piece_idx) {
            std::size_t num_decision_variables_of_piece =
                    piece_operations_ptrs_.at(piece_idx)->numDecisionVariables();
            Vector solution(num_decision_variables_of_piece);
            // get the decision variable value.
            for (std::size_t decision_variable_idx = 0;
                 decision_variable_idx < num_decision_variables_of_piece;
                 ++decision_variable_idx) {
                solution(decision_variable_idx) = variable_map_.at(piece_idx).at(decision_variable_idx)
                        ->solution_value();
            }

            std::unique_ptr<SingleParameterCurve> curve =
                    piece_operations_ptrs_.at(piece_idx)->generateCurveFromSolution(
                            solution);

            piecewise_curve.addPiece(std::move(curve));
        }

        return piecewise_curve;
    }

    template <typename T, unsigned int DIM>
    void
    SingleParameterPiecewiseCurveQPGenerator<T, DIM>::addIntegratedSquaredDerivativeCost(uint64_t derivative_degree,
                                                                                         T lambda) {
        for (std::size_t piece_idx = 0; piece_idx < piece_operations_ptrs_.size(); ++piece_idx) {
            const CostAddition &cost_addition = piece_operations_ptrs_.at(piece_idx)
                    ->integratedSquaredDerivativeCost(derivative_degree, lambda);
            addCostAdditionForPiece(piece_idx, cost_addition);
        }
    }

    template <typename T, unsigned int DIM>
    void SingleParameterPiecewiseCurveQPGenerator<T, DIM>::addEvalCost(T parameter, uint64_t derivative_degree,
                                                                       const VectorDIM &target, T lambda) {
        const PieceIndexAndParameter& piece_idx_and_parameter = getPieceIndexAndParameter(parameter);

        const CostAddition& cost_addition =
                piece_operations_ptrs_.at(piece_idx_and_parameter.piece_idx())->evalCost(
                        piece_idx_and_parameter.parameter(), derivative_degree, target, lambda);

        addCostAdditionForPiece(piece_idx_and_parameter.piece_idx(), cost_addition);
    }

    template <typename T, unsigned int DIM>
    void SingleParameterPiecewiseCurveQPGenerator<T, DIM>::addEvalConstraint(T parameter, uint64_t derivative_degree,
                                                                             const VectorDIM &target) {
        const PieceIndexAndParameter& piece_idx_and_parameter = getPieceIndexAndParameter(parameter);
        const std::vector<LinearConstraint> linear_constraints =
                piece_operations_ptrs_.at(piece_idx_and_parameter.piece_idx())
                ->evalConstraint(piece_idx_and_parameter.parameter(), derivative_degree, target);
        for (const LinearConstraint& linear_constraint : linear_constraints) {
            addLinearConstraintForPiece(piece_idx_and_parameter.piece_idx(), linear_constraint);
        }
    }

    template <typename T, unsigned int DIM>
    void SingleParameterPiecewiseCurveQPGenerator<T, DIM>::addCostAdditionForPiece(std::size_t piece_idx,
                                                                                   const CostAddition &cost_addition) {
        constexpr T epsilon = std::numeric_limits<T>::epsilon() * T(100.0);

        std::size_t num_decision_variables = variable_map_.at(piece_idx).size();

        // number decision variable number error
        if (num_decision_variables != cost_addition.linear_term().rows() ||
            num_decision_variables != cost_addition.quadratic_term().rows() ||
            num_decision_variables != cost_addition.quadratic_term().cols()) {
            throw std::runtime_error("SingleParameterPiecewiseCurveQPGenerator::addCostAdditionForPiece:"
                                     " number of decision variables of the piece does not match the "
                                     "CostAddition structure");
        }
        // piece_idx out_of_range error
        if (piece_idx >= piece_operations_ptrs_.size()) {
            throw std::out_of_range("SingleParameterPiecewiseCurveQPGenerator::addCostAdditionForPiece:"
                                    " piece index out of range. piece_idx: " + std::to_string(piece_idx) +
                                    ", num pieces: " + std::to_string(piece_operations_ptrs_.size()));
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
    void SingleParameterPiecewiseCurveQPGenerator<T, DIM>::addLinearConstraintForPiece(
            std::size_t piece_idx,
            const LinearConstraint &linear_constraint) {
        std::size_t num_decision_variables = variable_map_.at(piece_idx).size();

        // linear constraint structure error.
        if (num_decision_variables != linear_constraint.coefficients().cols()) {
            throw std::runtime_error("SingleParameterPiecewiseCurveQPGenerator::"
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

    template <typename T, unsigned int DIM>
    void SingleParameterPiecewiseCurveQPGenerator<T, DIM>::addContinuityConstraint(std::size_t piece_idx,
                                                                                   uint64_t derivative_degree) {
        using Row = math::Row<T>;
        // piece_idx out_of_range error
        if (piece_idx + 1 >= cumulative_max_parameters_.size()) {
            throw std::runtime_error("SingleParameterPiecewiseCurveQPGenerator::addContinuityConstraint:"
                                     " piece index out of range. piece_idx: " + std::to_string(piece_idx) +
                                     ", piece_idx + 1: " + std::to_string(piece_idx + 1) +
                                     ", number of pieces: " + std::to_string(cumulative_max_parameters_.size()));
        }

        // get the first piece eval basis (still need to multiply with control point of first piece)
        for (unsigned int dimension = 0; dimension < DIM; ++dimension) {
            const Row& first_piece_row =
                    piece_operations_ptrs_.at(piece_idx)->evalBasisRow(
                            /*dimension=*/dimension,
                            /*parameter=*/piece_operations_ptrs_.at(piece_idx)->max_parameter(),
                            /*derivative_degree=*/derivative_degree);
            const Row& second_piece_row =
                    piece_operations_ptrs_.at(piece_idx+1)->evalBasisRow(
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
    SingleParameterPiecewiseCurveQPGenerator<T, DIM>::addBoundingBoxConstraint(
            const AlignedBox& bounding_box) {
        for (std::size_t piece_idx = 0; piece_idx < piece_operations_ptrs_.size();
             ++piece_idx) {
            DecisionVariableBounds decision_variable_bounds =
                    piece_operations_ptrs_.at(piece_idx)->boundingBoxConstraint(
                            bounding_box);
            addDecisionVariableBoundsForPiece(piece_idx, decision_variable_bounds);
        }
    }

    template <typename T, unsigned int DIM>
    void
    SingleParameterPiecewiseCurveQPGenerator<T, DIM>::addBoundingBoxConstraintAll(
            const AlignedBox& aligned_box, uint64_t derivative_degree) {
        for (std::size_t piece_idx = 0; piece_idx < piece_operations_ptrs_.size();
             ++piece_idx) {
            const std::vector<LinearConstraint>& piece_constraints =
                    piece_operations_ptrs_.at(piece_idx)->boundingBoxConstraintAll(aligned_box, derivative_degree);

            for (const LinearConstraint& piece_constraint : piece_constraints) {
                addLinearConstraintForPiece(piece_idx, piece_constraint);
            }
        }
    }

    template <typename T, unsigned int DIM>
    void
    SingleParameterPiecewiseCurveQPGenerator<T, DIM>::addPositionConstraintForPiece(std::size_t piece_idx,
                                                                                    const VectorDIM &position) {
        if (piece_idx >= cumulative_max_parameters_.size()) {
            throw std::invalid_argument("SingleParameterPiecewiseCurveQPGenerator::"
                                        "addPositionConstraintForPiece: piece index out of range. "
                                        "piece_idx: " + std::to_string(piece_idx) +
                                        ", number of pieces: " + std::to_string(cumulative_max_parameters_.size()));
        }

        DecisionVariableBounds decision_variable_bounds =
                piece_operations_ptrs_.at(piece_idx)->positionConstraint(position);
        addDecisionVariableBoundsForPiece(piece_idx, decision_variable_bounds);
    }

    template <typename T, unsigned int DIM>
    void
    SingleParameterPiecewiseCurveQPGenerator<T, DIM>::addHyperplaneConstraintForPiece(std::size_t piece_idx,
                                                                                      const Hyperplane& hyperplane,
                                                                                      T epsilon) {
        if (piece_idx >= cumulative_max_parameters_.size()) {
            throw std::invalid_argument("SingleParameterPiecewiseCurveQPGenerator::"
                                        "addHyperplaneConstraintForPiece: piece index out of range. "
                                        "piece_idx: " + std::to_string(piece_idx) +
                                        ", number of pieces: " + std::to_string(cumulative_max_parameters_.size()));
        }

        const std::vector<LinearConstraint> linear_constraints =
                piece_operations_ptrs_.at(piece_idx)->hyperplaneConstraintAll(hyperplane, epsilon);

        for (const LinearConstraint& linear_constraint : linear_constraints) {
            addLinearConstraintForPiece(piece_idx, linear_constraint);
        }

    }

    template <typename T, unsigned int DIM>
    void SingleParameterPiecewiseCurveQPGenerator<T, DIM>::addDecisionVariableBoundsForPiece(
            std::size_t piece_idx,
            const DecisionVariableBounds& decision_variable_bounds) {
        std::size_t num_decision_variables = variable_map_.at(piece_idx).size();
        if (num_decision_variables != decision_variable_bounds.lower_bounds().rows() ||
            num_decision_variables != decision_variable_bounds.upper_bounds().rows()) {
            throw std::invalid_argument("SingleParameterPiecewiseCurveQPGenerator::"
                                        "addDecisionVariableBoundsForPiece:"
                                        " number of decision variables of the piece does not match the "
                                        "DecisionVariablesBounds structure");
        }

        for (std::size_t decision_variable_idx = 0;
             decision_variable_idx < num_decision_variables;
             ++decision_variable_idx) {
            qpcpp::Variable<T>* var_ptr =
                    variable_map_.at(piece_idx).at(decision_variable_idx);

            var_ptr->set_min(std::max(
                    var_ptr->min(),
                    decision_variable_bounds.lower_bounds()(decision_variable_idx)));

            var_ptr->set_max(std::min(
                    var_ptr->max(),
                    decision_variable_bounds.upper_bounds()(decision_variable_idx)));
        }
    }

    template <typename T, unsigned int DIM>
    qpcpp::Problem<T>& SingleParameterPiecewiseCurveQPGenerator<T, DIM>::problem() {
        return problem_;
    }

    template <typename T, unsigned int DIM>
    typename SingleParameterPiecewiseCurveQPGenerator<T, DIM>::PieceIndexAndParameter
    SingleParameterPiecewiseCurveQPGenerator<T, DIM>::getPieceIndexAndParameter(T parameter) const {
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

    template <typename T, unsigned int DIM>
    SingleParameterPiecewiseCurveQPGenerator<T, DIM>::PieceIndexAndParameter::
        PieceIndexAndParameter(std::size_t piece_idx, T parameter)
        : piece_idx_(piece_idx), parameter_(parameter) {}

    template <typename T, unsigned int DIM>
    std::size_t
    SingleParameterPiecewiseCurveQPGenerator<T, DIM>::PieceIndexAndParameter::piece_idx() const {
        return piece_idx_;
    }

    template <typename T, unsigned int DIM>
    T
    SingleParameterPiecewiseCurveQPGenerator<T, DIM>::PieceIndexAndParameter::parameter() const {
        return parameter_;
    }

    template class SingleParameterPiecewiseCurveQPGenerator<float, 2U>;
    template class SingleParameterPiecewiseCurveQPGenerator<double, 2U>;
    template class SingleParameterPiecewiseCurveQPGenerator<float, 3U>;
    template class SingleParameterPiecewiseCurveQPGenerator<double, 3U>;
} // splines