//
// Created by lishuo on 4/11/24.
//

#include <separating_hyperplanes/SVM.h>

namespace separating_hyperplanes {
    template <typename T, unsigned int DIM>
    Hyperplane<T, DIM> hardMarginSVM(
            const std::vector<VectorDIM<T, DIM>>& first_set_of_points,
            const std::vector<VectorDIM<T, DIM>>& second_set_of_points) {
        // Input validation
        if (first_set_of_points.empty()) {
            throw std::invalid_argument("hardMarginSVM: first set of points is empty.");
        }

        if (second_set_of_points.empty()) {
            throw std::invalid_argument("hardMarginSVM: second set of points is empty.");
        }

        // Setup QP problem to find the maximum margin hyperplane
        // The variables are the normal vector components and the offset
        qpcpp::Problem<T> problem;
        std::vector<qpcpp::Variable<T>*> variable_ptrs(DIM + 1);

        // Create variables for hyperplane: DIM components for normal vector and 1 for offset
        for (unsigned int d = 0; d < DIM; ++d) {
            variable_ptrs[d] = problem.addVariable();
        }
        variable_ptrs[DIM] = problem.addVariable(); // Offset term

        // Add constraints for first set of points: w·x + b ≤ -1
        // (points must be on the negative side with margin)
        for (const VectorDIM<T, DIM>& point : first_set_of_points) {
            qpcpp::LinearConstraint<T>* linear_constraint =
                    problem.addLinearConstraint(std::numeric_limits<T>::lowest(), -1);
            for (unsigned int d = 0; d < DIM; d++) {
                linear_constraint->setCoefficient(variable_ptrs[d], point(d));
            }
            linear_constraint->setCoefficient(variable_ptrs[DIM], 1);
        }

        // Add constraints for second set of points: w·x + b ≥ 1
        // (points must be on the positive side with margin)
        for (const VectorDIM<T, DIM>& point : second_set_of_points) {
            qpcpp::LinearConstraint<T>* linear_constraint =
                    problem.addLinearConstraint(1, std::numeric_limits<T>::max());
            for (unsigned int d = 0; d < DIM; d++) {
                linear_constraint->setCoefficient(variable_ptrs[d], point(d));
            }
            linear_constraint->setCoefficient(variable_ptrs[DIM], 1);
        }

        // Set objective function to minimize ||w||²/2 (maximize margin)
        for (unsigned int d = 0; d < DIM; ++d) {
            problem.cost_function()->addQuadraticTerm(variable_ptrs[d],
                                                      variable_ptrs[d], 1);
        }

        // Solve QP problem
        qpcpp::CPLEXSolver<T> solver;
        qpcpp::SolveStatus solve_status = solver.solve(problem);

        // Extract solution if available
        if (solve_status == qpcpp::SolveStatus::OPTIMAL ||
            solve_status == qpcpp::SolveStatus::FEASIBLE) {
            Hyperplane<T, DIM> hyperplane;
            for (unsigned int d = 0; d < DIM; ++d) {
                hyperplane.normal()(d) = variable_ptrs[d]->solution_value();
            }
            hyperplane.offset() = variable_ptrs.back()->solution_value();
            return hyperplane;
        }

        // No solution found
        throw std::invalid_argument("SVM is not feasible. solve_status: " +
                                    qpcpp::SolveStatusToStr(solve_status));
    }

    // Explicit template instantiations
    template Hyperplane<double, 3U> hardMarginSVM<double, 3U>(
            const std::vector<VectorDIM<double, 3U>>& first_set_of_points,
            const std::vector<VectorDIM<double, 3U>>& second_set_of_points);

    template Hyperplane<float, 3U> hardMarginSVM<float, 3U>(
            const std::vector<VectorDIM<float, 3U>>& first_set_of_points,
            const std::vector<VectorDIM<float, 3U>>& second_set_of_points);

    template Hyperplane<double, 2U> hardMarginSVM<double, 2U>(
            const std::vector<VectorDIM<double, 2U>>& first_set_of_points,
            const std::vector<VectorDIM<double, 2U>>& second_set_of_points);

    template Hyperplane<float, 2U> hardMarginSVM<float, 2U>(
            const std::vector<VectorDIM<float, 2U>>& first_set_of_points,
            const std::vector<VectorDIM<float, 2U>>& second_set_of_points);

} // separating_hyperplanes