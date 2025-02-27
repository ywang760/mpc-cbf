#ifndef MATH_COLLISIONSHAPES_ALIGNED_BOX_H
#define MATH_COLLISIONSHAPES_ALIGNED_BOX_H

#include <math/collision_shapes/CollisionShape.h>

namespace math {

    template <typename T, unsigned int DIM>
    class AlignedBoxCollisionShape : public CollisionShape<T, DIM> {
    public:
        using Base = CollisionShape<T, DIM>;
        using VectorDIM = typename Base::VectorDIM;
        using AlignedBox = typename Base::AlignedBox;

        // Constructs AlignedBoxCollisionShape where the collision shape is
        // collision_box_at_zero at zero.
        AlignedBoxCollisionShape(const AlignedBox& collision_box_at_zero);

        virtual std::vector<VectorDIM> convexHullPoints(
                const VectorDIM& position) const override;

        virtual AlignedBox boundingBox(const VectorDIM& position) const override;

        const AlignedBox& collision_box_at_zero() const {
            return collision_box_at_zero_;
        }

        /**
         * @brief Inflate the aligned box collision shape by inflation_amound in
         * primary directions
         *
         * @param inflation_amount Amount by which the collision shape is inflated
         * @return std::shared_ptr<CollisionShape<T,DIM>> Inflated collision shape
         */
        virtual std::shared_ptr<CollisionShape<T, DIM>> inflate(
                T inflation_amount) const override;

        bool operator==(const AlignedBoxCollisionShape& other) const {
            return collision_box_at_zero_.min() ==
                   other.collision_box_at_zero_.min() &&
                   collision_box_at_zero_.max() ==
                   other.collision_box_at_zero_.max();
        }

    private:
        // collision box of the robot when center of mass is at 0.
        AlignedBox collision_box_at_zero_;
    };
}  // namespace math

#endif  // MATH_COLLISIONSHAPES_ALIGNED_BOX_H