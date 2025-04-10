//
// Created by lishuo on 8/27/24.
//

#ifndef SEPARATINGHYPERPLANES_SVM_H
#define SEPARATINGHYPERPLANES_SVM_H

#include <separating_hyperplanes/Types.h>
#include <qpcpp/Problem.h>
#include <qpcpp/solvers/CPLEX.h>

namespace separating_hyperplanes {
    /**
     * Computes a hard-margin SVM hyperplane separating two sets of points.
     *
     * Unlike Voronoi separating hyperplanes which simply find an equidistant
     * boundary between point pairs, the hard-margin SVM computes the optimal
     * hyperplane with maximum margin between classes by solving a quadratic
     * programming problem. The resulting hyperplane:
     *
     * 1. Maximizes the distance to the closest point from either set
     * 2. Creates a global decision boundary that relies only on support vectors
     * 3. Guarantees all points in first_set_of_points lie on the negative side
     *    and all points in second_set_of_points lie on the positive side
     *
     * This approach requires that the two sets are linearly separable.
     *
     * @param first_set_of_points Points to be classified as negative (w·x + b < 0)
     * @param second_set_of_points Points to be classified as positive (w·x + b > 0)
     * @return Optimal separating hyperplane that maximizes margin between classes
     * @throws std::invalid_argument if sets are empty or linearly inseparable
     */
    template <typename T, unsigned int DIM>
    Hyperplane<T, DIM> hardMarginSVM(
            const std::vector<VectorDIM<T, DIM>>& first_set_of_points,
            const std::vector<VectorDIM<T, DIM>>& second_set_of_points);

} // separating_hyperplanes

#endif //SEPARATINGHYPERPLANES_SVM_H
