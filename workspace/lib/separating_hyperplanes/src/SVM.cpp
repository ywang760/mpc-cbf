//
// Created by lishuo on 4/11/24.
//

#include <separating_hyperplanes/SVM.h>

namespace separating_hyperplanes {
    template <typename T, unsigned int DIM>
    Hyperplane<T, DIM> hardMarginSVM(
            const std::vector<VectorDIM<T, DIM>>& first_set_of_points,
            const std::vector<VectorDIM<T, DIM>>& second_set_of_points) {
        if (first_set_of_points.empty()) {
            throw std::invalid_argument("hardMarginSVM: first set of points is empty.");
        }

        if (second_set_of_points.empty()) {
            throw std::invalid_argument("hardMarginSVM: second set of points is empty.");
        }

        qpcpp::Problem<T> problem;
        std::vector<qpcpp::Variable<T>*> variable_ptrs(DIM + 1);

        for (unsigned int d = 0; d < DIM; ++d) {
            variable_ptrs[d] = problem.addVariable();
        }
        variable_ptrs[DIM] = problem.addVariable();

        for (const VectorDIM<T, DIM>& point : first_set_of_points) {
            qpcpp::LinearConstraint<T>* linear_constraint =
                    problem.addLinearConstraint(std::numeric_limits<T>::lowest(), -1);
            for (unsigned int d = 0; d < DIM; d++) {
                linear_constraint->setCoefficient(variable_ptrs[d], point(d));
            }
            linear_constraint->setCoefficient(variable_ptrs[DIM], 1);
        }

        for (const VectorDIM<T, DIM>& point : second_set_of_points) {
            qpcpp::LinearConstraint<T>* linear_constraint =
                    problem.addLinearConstraint(1, std::numeric_limits<T>::max());
            for (unsigned int d = 0; d < DIM; d++) {
                linear_constraint->setCoefficient(variable_ptrs[d], point(d));
            }
            linear_constraint->setCoefficient(variable_ptrs[DIM], 1);
        }

        for (unsigned int d = 0; d < DIM; ++d) {
            problem.cost_function()->addQuadraticTerm(variable_ptrs[d],
                                                      variable_ptrs[d], 1);
        }

        qpcpp::CPLEXSolver<T> solver;
        qpcpp::SolveStatus solve_status = solver.solve(problem);

        if (solve_status == qpcpp::SolveStatus::OPTIMAL ||
            solve_status == qpcpp::SolveStatus::FEASIBLE) {
            Hyperplane<T, DIM> hyperplane;
            for (unsigned int d = 0; d < DIM; ++d) {
                hyperplane.normal()(d) = variable_ptrs[d]->solution_value();
            }
            hyperplane.offset() = variable_ptrs.back()->solution_value();
            return hyperplane;
        }

        throw std::invalid_argument("SVM is not feasible. solve_status: " +
                                    qpcpp::SolveStatusToStr(solve_status));
    }

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