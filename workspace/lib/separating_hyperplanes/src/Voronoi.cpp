//
// Created by lishuo on 8/27/24.
//

#include <separating_hyperplanes/Voronoi.h>

namespace separating_hyperplanes {

    template <typename T, unsigned int DIM>
    Hyperplane<T, DIM> voronoi(
            const VectorDIM<T, DIM>& first_position,
            const VectorDIM<T, DIM>& second_position) {
        using VectorDIM = VectorDIM<T, DIM>;

        // Calculate the direction vector from first to second point
        VectorDIM direction_from_first_to_second =
                second_position - first_position;
        direction_from_first_to_second.normalize();

        // Voronoi hyperplanes pass through the midpoint between the points
        // and are perpendicular to the line connecting them
        VectorDIM mid_point = (first_position + second_position) / 2.0;

        // Calculate the offset by projecting the midpoint onto the normal direction
        // Note: This creates a hyperplane equidistant from both points
        T negative_offset = direction_from_first_to_second.dot(mid_point);

        // Return hyperplane in standard form: n·x + d = 0
        // where n is the normal vector and d is the offset
        return Hyperplane<T, DIM>(direction_from_first_to_second, -negative_offset);
    }

    // Explicit template instantiations
    template Hyperplane<double, 3U> voronoi<double, 3U>(
            const VectorDIM<double, 3U>&, const VectorDIM<double, 3U>&);
    template Hyperplane<float, 3U> voronoi<float, 3U>(
            const VectorDIM<float, 3U>&, const VectorDIM<float, 3U>&);
    template Hyperplane<double, 2U> voronoi<double, 2U>(
            const VectorDIM<double, 2U>&, const VectorDIM<double, 2U>&);
    template Hyperplane<float, 2U> voronoi<float, 2U>(
            const VectorDIM<float, 2U>&, const VectorDIM<float, 2U>&);

}  // namespace separating_hyperplanes