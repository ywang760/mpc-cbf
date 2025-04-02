//
// Created by lishuo on 8/27/24.
//

#ifndef SEPARATINGHYPERPLANES_VORONOI_H
#define SEPARATINGHYPERPLANES_VORONOI_H

#include <separating_hyperplanes/Types.h>

namespace separating_hyperplanes {
    /**
     * Computes a Voronoi separating hyperplane between two points.
     * 
     * The Voronoi hyperplane is the perpendicular bisector of the line segment
     * connecting the two input points. This hyperplane:
     * 
     * 1. Passes through the midpoint between the two input points
     * 2. Has a normal vector pointing from first_position to second_position
     * 3. Partitions the space such that any point on one side is closer to 
     *    first_position, and any point on the other side is closer to second_position
     * 
     * Unlike SVM hyperplanes which solve a global optimization problem,
     * Voronoi hyperplanes are constructed geometrically and are based purely on
     * the local relationship between two points.
     * 
     * @param first_position First point (will be on the negative side of hyperplane)
     * @param second_position Second point (will be on the positive side of hyperplane)
     * @return Hyperplane that is equidistant from both input points
     */
    template <typename T, unsigned int DIM>
    Hyperplane<T, DIM> voronoi(
            const VectorDIM<T, DIM>& first_position,
            const VectorDIM<T, DIM>& second_position);

} // separating_hyperplanes

#endif //SEPARATINGHYPERPLANES_VORONOI_H
