#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "math/collision_shapes/CollisionShape.h"
#include "math/Helpers.h"
#include <unordered_map>
#include <memory>
#include <stdexcept>

// Mock implementation of CollisionShape for testing
template <typename T, unsigned int DIM>
class MockCollisionShape : public math::CollisionShape<T, DIM> {
public:
    using VectorDIM = typename math::CollisionShape<T, DIM>::VectorDIM;
    using AlignedBox = typename math::CollisionShape<T, DIM>::AlignedBox;
    
    MockCollisionShape() : math::CollisionShape<T, DIM>(math::CollisionShape<T, DIM>::Type::ALIGNED_BOX) {}
    
    std::vector<VectorDIM> convexHullPoints(const VectorDIM& position) const override {
        std::vector<VectorDIM> points;
        // Create a simple box around the position
        AlignedBox box(position - VectorDIM::Ones(), position + VectorDIM::Ones());
        return math::cornerPoints<T, DIM>(box);
    }
    
    AlignedBox boundingBox(const VectorDIM& position) const override {
        return AlignedBox(position - VectorDIM::Ones(), position + VectorDIM::Ones());
    }
    
    std::shared_ptr<math::CollisionShape<T, DIM>> inflate(T inflation_amount) const override {
        auto inflated = std::make_shared<MockCollisionShape<T, DIM>>(*this);
        return inflated;
    }
};

class CollisionShapeTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Setup code before each test if needed
    }

    void TearDown() override {
        // Cleanup code after each test if needed
    }
};

TEST_F(CollisionShapeTest, TestCollisionShapeType) {
    // Test type to string conversions
    EXPECT_EQ((math::CollisionShape<double, 2>::typeToString(
        math::CollisionShape<double, 2>::Type::ALIGNED_BOX)), "aligned_box");
    
    // Test string to type conversions
    EXPECT_EQ((math::CollisionShape<double, 2>::stringToType("aligned_box")),
        (math::CollisionShape<double, 2>::Type::ALIGNED_BOX));
    
    // Test for exception on invalid string
    EXPECT_THROW((math::CollisionShape<double, 2>::stringToType("invalid_type")), 
        std::out_of_range);
}

TEST_F(CollisionShapeTest, TestMockCollisionShape2D) {
    // Test with 2D collision shape
    MockCollisionShape<double, 2> shape;
    
    // Test type
    EXPECT_EQ((shape.type()), (math::CollisionShape<double, 2>::Type::ALIGNED_BOX));
    
    // Test convex hull points
    math::VectorDIM<double, 2> position(1.0, 2.0);
    auto points = shape.convexHullPoints(position);
    
    // Expect 4 corner points for a 2D box
    EXPECT_EQ(points.size(), 4);
    
    // Check bounding box
    auto box = shape.boundingBox(position);
    EXPECT_TRUE((math::isApproximatelyEqual<double, 2>(box.min(), position - math::VectorDIM<double, 2>::Ones(), 1e-10)));
    EXPECT_TRUE((math::isApproximatelyEqual<double, 2>(box.max(), position + math::VectorDIM<double, 2>::Ones(), 1e-10)));
    
    // Test inflate
    auto inflated = shape.inflate(0.5);
    EXPECT_NE(inflated, nullptr);
    EXPECT_EQ((inflated->type()), (math::CollisionShape<double, 2>::Type::ALIGNED_BOX));
}

TEST_F(CollisionShapeTest, TestMockCollisionShape3D) {
    // Test with 3D collision shape
    MockCollisionShape<double, 3> shape;
    
    // Test type
    EXPECT_EQ((shape.type()), (math::CollisionShape<double, 3>::Type::ALIGNED_BOX));
    
    // Test convex hull points
    math::VectorDIM<double, 3> position(1.0, 2.0, 3.0);
    auto points = shape.convexHullPoints(position);
    
    // Expect 8 corner points for a 3D box
    EXPECT_EQ(points.size(), 8);
    
    // Check bounding box
    auto box = shape.boundingBox(position);
    EXPECT_TRUE((math::isApproximatelyEqual<double, 3>(box.min(), position - math::VectorDIM<double, 3>::Ones(), 1e-10)));
    EXPECT_TRUE((math::isApproximatelyEqual<double, 3>(box.max(), position + math::VectorDIM<double, 3>::Ones(), 1e-10)));
    
    // Test inflate
    auto inflated = shape.inflate(0.5);
    EXPECT_NE(inflated, nullptr);
    EXPECT_EQ((inflated->type()), (math::CollisionShape<double, 3>::Type::ALIGNED_BOX));
}

TEST_F(CollisionShapeTest, TestDifferentTypes) {
    // Test with float instead of double
    MockCollisionShape<float, 2> shape_float_2d;
    MockCollisionShape<float, 3> shape_float_3d;
    
    // Basic type checks
    EXPECT_EQ((shape_float_2d.type()), (math::CollisionShape<float, 2>::Type::ALIGNED_BOX));
    EXPECT_EQ((shape_float_3d.type()), (math::CollisionShape<float, 3>::Type::ALIGNED_BOX));
    
    // Test with different positions
    math::VectorDIM<float, 2> position_2d(1.0f, 2.0f);
    math::VectorDIM<float, 3> position_3d(1.0f, 2.0f, 3.0f);
    
    auto points_2d = shape_float_2d.convexHullPoints(position_2d);
    auto points_3d = shape_float_3d.convexHullPoints(position_3d);
    
    EXPECT_EQ(points_2d.size(), 4);
    EXPECT_EQ(points_3d.size(), 8);
}
