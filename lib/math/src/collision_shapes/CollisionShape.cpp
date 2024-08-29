#include <math/collision_shapes/CollisionShape.h>

namespace math {
    template <typename T, unsigned int DIM>
    typename CollisionShape<T, DIM>::Type CollisionShape<T, DIM>::stringToType(
            const std::string& type_string) {
        static const std::unordered_map<std::string, Type> type_map = {
                {"aligned_box", Type::ALIGNED_BOX},
        };

        return type_map.at(type_string);
    }

    template <typename T, unsigned int DIM>
    const std::string& CollisionShape<T, DIM>::typeToString(Type type) {
        static const std::unordered_map<Type, std::string> type_map = {
                {Type::ALIGNED_BOX, "aligned_box"},
        };

        return type_map.at(type);
    }

    template class CollisionShape<double, 2U>;
    template class CollisionShape<float, 2U>;
    template class CollisionShape<double, 3U>;
    template class CollisionShape<float, 3U>;

}  // namespace math