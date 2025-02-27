//
// Created by lishuo on 8/27/24.
//

#ifndef SEPARATINGHYPERPLANES_VORONOI_H
#define SEPARATINGHYPERPLANES_VORONOI_H

#include <separating_hyperplanes/Types.h>

namespace separating_hyperplanes {
    // returns the voronoi hyperplane between first_hypersphere and
    // second_hypersphere such that first_hypersphere is on the non-positive side of
    // the hyperplane. return status is not ok if hyperspheres are intersecting
    template <typename T, unsigned int DIM>
    Hyperplane<T, DIM> voronoi(
            const VectorDIM<T, DIM>& first_position,
            const VectorDIM<T, DIM>& second_position);
} // separating_hyperplanes

#endif //SEPARATINGHYPERPLANES_VORONOI_H
