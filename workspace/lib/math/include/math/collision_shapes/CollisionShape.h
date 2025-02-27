#ifndef MATH_COLLISION_SHAPES_COLLISION_SHAPE_H
#define MATH_COLLISION_SHAPES_COLLISION_SHAPE_H

#include <math/Types.h>
#include <memory>

namespace math {

    template <typename T, unsigned int DIM>
    class CollisionShape {
    public:
        using VectorDIM = math::VectorDIM<T, DIM>;
        using AlignedBox = math::AlignedBox<T, DIM>;

        /**
         * @brief Types of collision shapes
         *
         */
        enum class Type { ALIGNED_BOX };

        /**
         * @brief Convert type to its string representation
         *
         * @param type Type to be converted
         * @return const std::string& String representation of the type
         */
        static const std::string& typeToString(Type type);

        /**
         * @brief Convert string representation of type to its enum value
         *
         * @param type_string String representation of the type
         * @return Type Type enum value
         */
        static Type stringToType(const std::string& type_string);

        /**
         * @brief Construct a new Collision Shape object
         *
         * @param type Type of the collision shape
         */
        CollisionShape(Type type) : type_(type) {}

        /**
         * @brief Get the type of the collision shape
         *
         * @return Type Type of the collision shape
         */
        Type type() const { return type_; }

        /**
         * @brief Destroy the Collision Shape object
         *
         */
        virtual ~CollisionShape() = default;

        // Returns a set of points, convex hull of which contraints the collision
        // shape at the given position
        virtual std::vector<VectorDIM> convexHullPoints(
                const VectorDIM& position) const = 0;

        // Returns an axis aligned bounding box of the collision shape at position
        virtual AlignedBox boundingBox(const VectorDIM& position) const = 0;


        /**
         * @brief Inflate the collision shape by the given amount
         *
         * @param inflation_amount Amount by which the collision shape is inflated
         * @return std::shared_ptr<CollisionShape<T,DIM>> Inflated collision shape
         */
        virtual std::shared_ptr<CollisionShape<T, DIM>> inflate(
                T inflation_amount) const = 0;

    private:
        Type type_;
    };
}  // namespace math

#endif  // MATH_COLLISION_SHAPES_COLLISION_SHAPE_H