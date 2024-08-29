#include <math/collision_shapes/AlignedBoxCollisionShape.h>
#include <math/Helpers.h>

#include <iomanip>

namespace math {

    template <typename T, unsigned int DIM>
    AlignedBoxCollisionShape<T, DIM>::AlignedBoxCollisionShape(
            const typename AlignedBoxCollisionShape<T, DIM>::AlignedBox&
            collision_box_at_zero)
            : CollisionShape<T, DIM>(CollisionShape<T, DIM>::Type::ALIGNED_BOX),
              collision_box_at_zero_(collision_box_at_zero) {}

    template <typename T, unsigned int DIM>
    std::vector<typename AlignedBoxCollisionShape<T, DIM>::VectorDIM>
    AlignedBoxCollisionShape<T, DIM>::convexHullPoints(
            const typename AlignedBoxCollisionShape<T, DIM>::VectorDIM& position)
    const {

        AlignedBox box(collision_box_at_zero_.min() + position,
                       collision_box_at_zero_.max() + position);
        return math::cornerPoints<T, DIM>(box);
    }

    template <typename T, unsigned int DIM>
    typename AlignedBoxCollisionShape<T, DIM>::AlignedBox
    AlignedBoxCollisionShape<T, DIM>::boundingBox(
            const typename AlignedBoxCollisionShape<T, DIM>::VectorDIM& position)
    const {

        return AlignedBox(collision_box_at_zero_.min() + position,
                          collision_box_at_zero_.max() + position);
    }

    template <typename T, unsigned int DIM>
    std::shared_ptr<CollisionShape<T, DIM>>
    AlignedBoxCollisionShape<T, DIM>::inflate(T inflation_amount) const {
        const VectorDIM inflation_vector = VectorDIM::Constant(inflation_amount);
        const AlignedBox inflated_box_at_zero =
                AlignedBox(collision_box_at_zero_.min() - inflation_vector,
                           collision_box_at_zero_.max() + inflation_vector);

        return std::make_shared<AlignedBoxCollisionShape<T, DIM>>(
                inflated_box_at_zero);
    }

    template class AlignedBoxCollisionShape<double, 2U>;
    template class AlignedBoxCollisionShape<float, 2U>;
    template class AlignedBoxCollisionShape<double, 3U>;
    template class AlignedBoxCollisionShape<float, 3U>;
}  // namespace math