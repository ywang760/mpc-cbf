//
// Created by lishuo on 8/27/24.
//

#ifndef SEPARATINGHYPERPLANES_SVM_H
#define SEPARATINGHYPERPLANES_SVM_H

#include <separating_hyperplanes/Types.h>
#include <qpcpp/Problem.h>
#include <qpcpp/solvers/CPLEX.h>

namespace separating_hyperplanes {
    // returns hard-margin svm hyperplane between first set of points and the second
    // set of points. first_set_of_points remain in the negative side of the
    // returned hyperplane. return status is not ok if any set of points is empty or
    // no hard-margin svm exists between two sets
    template <typename T, unsigned int DIM>
    Hyperplane<T, DIM> hardMarginSVM(
            const std::vector<VectorDIM<T, DIM>>& first_set_of_points,
            const std::vector<VectorDIM<T, DIM>>& second_set_of_points);

} // separating_hyperplanes

#endif //SEPARATINGHYPERPLANES_SVM_H
