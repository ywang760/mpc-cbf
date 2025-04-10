#include <gtest/gtest.h>
#include "math/collision_shapes/AlignedBoxCollisionShape.h"
#include "math/Helpers.h"

class AlignedBoxCollisionShapeTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Setup common test objects
        box2d_min = math::VectorDIM<double, 2>(-1.0, -2.0);
        box2d_max = math::VectorDIM<double, 2>(1.0, 2.0);
        box2d = math::AlignedBox<double, 2>(box2d_min, box2d_max);
        
        box3d_min = math::VectorDIM<double, 3>(-1.5, -2.5, -3.5);
        box3d_max = math::VectorDIM<double, 3>(1.5, 2.5, 3.5);
        box3d = math::AlignedBox<double, 3>(box3d_min, box3d_max);
    }

    // Common test objects
    math::VectorDIM<double, 2> box2d_min, box2d_max;
    math::AlignedBox<double, 2> box2d;
    
    math::VectorDIM<double, 3> box3d_min, box3d_max;
    math::AlignedBox<double, 3> box3d;
};

// Test construction and access to collision_box_at_zero
TEST_F(AlignedBoxCollisionShapeTest, Construction) {
    // 2D test
    math::AlignedBoxCollisionShape<double, 2> shape_2d(box2d);
    EXPECT_TRUE(math::isApproximatelyEqual(box2d_min, shape_2d.collision_box_at_zero().min(), 1e-10));
    EXPECT_TRUE(math::isApproximatelyEqual(box2d_max, shape_2d.collision_box_at_zero().max(), 1e-10));
    
    // 3D test
    math::AlignedBoxCollisionShape<double, 3> shape_3d(box3d);
    EXPECT_TRUE(math::isApproximatelyEqual(box3d_min, shape_3d.collision_box_at_zero().min(), 1e-10));
    EXPECT_TRUE(math::isApproximatelyEqual(box3d_max, shape_3d.collision_box_at_zero().max(), 1e-10));
}

// Test convexHullPoints method
TEST_F(AlignedBoxCollisionShapeTest, ConvexHullPoints) {
    // 2D test
    math::AlignedBoxCollisionShape<double, 2> shape_2d(box2d);
    math::VectorDIM<double, 2> position_2d(2.0, 3.0);
    auto hull_points_2d = shape_2d.convexHullPoints(position_2d);
    
    // Should have 4 corners in 2D
    EXPECT_EQ(4, hull_points_2d.size());
    
    // Check that the points are correctly positioned
    math::AlignedBox<double, 2> expected_box_2d(box2d_min + position_2d, box2d_max + position_2d);
    auto expected_corners_2d = math::cornerPoints<double, 2>(expected_box_2d);
    
    for (size_t i = 0; i < hull_points_2d.size(); ++i) {
        bool found_match = false;
        for (size_t j = 0; j < expected_corners_2d.size(); ++j) {
            if (math::isApproximatelyEqual(hull_points_2d[i], expected_corners_2d[j], 1e-10)) {
                found_match = true;
                break;
            }
        }
        EXPECT_TRUE(found_match) << "Point " << i << " not found in expected corners";
    }
    
    // 3D test
    math::AlignedBoxCollisionShape<double, 3> shape_3d(box3d);
    math::VectorDIM<double, 3> position_3d(1.0, 2.0, 3.0);
    auto hull_points_3d = shape_3d.convexHullPoints(position_3d);
    
    // Should have 8 corners in 3D
    EXPECT_EQ(8, hull_points_3d.size());
    
    // Check that the points are correctly positioned
    math::AlignedBox<double, 3> expected_box_3d(box3d_min + position_3d, box3d_max + position_3d);
    auto expected_corners_3d = math::cornerPoints<double, 3>(expected_box_3d);
    
    for (size_t i = 0; i < hull_points_3d.size(); ++i) {
        bool found_match = false;
        for (size_t j = 0; j < expected_corners_3d.size(); ++j) {
            if (math::isApproximatelyEqual(hull_points_3d[i], expected_corners_3d[j], 1e-10)) {
                found_match = true;
                break;
            }
        }
        EXPECT_TRUE(found_match) << "Point " << i << " not found in expected corners";
    }
}

// Test boundingBox method
TEST_F(AlignedBoxCollisionShapeTest, BoundingBox) {
    // 2D test
    math::AlignedBoxCollisionShape<double, 2> shape_2d(box2d);
    math::VectorDIM<double, 2> position_2d(2.0, 3.0);
    auto bbox_2d = shape_2d.boundingBox(position_2d);
    
    math::VectorDIM<double, 2> expected_min_2d = box2d_min + position_2d;
    math::VectorDIM<double, 2> expected_max_2d = box2d_max + position_2d;
    
    EXPECT_TRUE(math::isApproximatelyEqual(expected_min_2d, bbox_2d.min(), 1e-10));
    EXPECT_TRUE(math::isApproximatelyEqual(expected_max_2d, bbox_2d.max(), 1e-10));
    
    // 3D test
    math::AlignedBoxCollisionShape<double, 3> shape_3d(box3d);
    math::VectorDIM<double, 3> position_3d(1.0, 2.0, 3.0);
    auto bbox_3d = shape_3d.boundingBox(position_3d);
    
    math::VectorDIM<double, 3> expected_min_3d = box3d_min + position_3d;
    math::VectorDIM<double, 3> expected_max_3d = box3d_max + position_3d;
    
    EXPECT_TRUE(math::isApproximatelyEqual(expected_min_3d, bbox_3d.min(), 1e-10));
    EXPECT_TRUE(math::isApproximatelyEqual(expected_max_3d, bbox_3d.max(), 1e-10));
}

// Test inflate method
TEST_F(AlignedBoxCollisionShapeTest, Inflate) {
    // 2D test
    math::AlignedBoxCollisionShape<double, 2> shape_2d(box2d);
    double inflation_amount = 0.5;
    auto inflated_shape_2d = shape_2d.inflate(inflation_amount);
    
    // Verify type
    EXPECT_EQ((math::CollisionShape<double, 2>::Type::ALIGNED_BOX), (inflated_shape_2d->type()));
    
    // Cast back to AlignedBoxCollisionShape to check properties
    auto* inflated_box_2d = dynamic_cast<math::AlignedBoxCollisionShape<double, 2>*>(inflated_shape_2d.get());
    ASSERT_NE(nullptr, inflated_box_2d);
    
    math::VectorDIM<double, 2> expected_min_2d = box2d_min - math::VectorDIM<double, 2>::Constant(inflation_amount);
    math::VectorDIM<double, 2> expected_max_2d = box2d_max + math::VectorDIM<double, 2>::Constant(inflation_amount);
    
    EXPECT_TRUE(math::isApproximatelyEqual(expected_min_2d, inflated_box_2d->collision_box_at_zero().min(), 1e-10));
    EXPECT_TRUE(math::isApproximatelyEqual(expected_max_2d, inflated_box_2d->collision_box_at_zero().max(), 1e-10));
    
    // 3D test
    math::AlignedBoxCollisionShape<double, 3> shape_3d(box3d);
    auto inflated_shape_3d = shape_3d.inflate(inflation_amount);
    
    // Verify type
    EXPECT_EQ((math::CollisionShape<double, 3>::Type::ALIGNED_BOX), (inflated_shape_3d->type()));
    
    // Cast back to AlignedBoxCollisionShape to check properties
    auto* inflated_box_3d = dynamic_cast<math::AlignedBoxCollisionShape<double, 3>*>(inflated_shape_3d.get());
    ASSERT_NE(nullptr, inflated_box_3d);
    
    math::VectorDIM<double, 3> expected_min_3d = box3d_min - math::VectorDIM<double, 3>::Constant(inflation_amount);
    math::VectorDIM<double, 3> expected_max_3d = box3d_max + math::VectorDIM<double, 3>::Constant(inflation_amount);
    
    EXPECT_TRUE(math::isApproximatelyEqual(expected_min_3d, inflated_box_3d->collision_box_at_zero().min(), 1e-10));
    EXPECT_TRUE(math::isApproximatelyEqual(expected_max_3d, inflated_box_3d->collision_box_at_zero().max(), 1e-10));
}

// Test equality operator
TEST_F(AlignedBoxCollisionShapeTest, EqualityOperator) {
    math::AlignedBoxCollisionShape<double, 2> shape1(box2d);
    math::AlignedBoxCollisionShape<double, 2> shape2(box2d);
    math::AlignedBoxCollisionShape<double, 2> shape3(math::AlignedBox<double, 2>(
        math::VectorDIM<double, 2>(-2.0, -3.0), 
        math::VectorDIM<double, 2>(2.0, 3.0)
    ));
    
    EXPECT_TRUE(shape1 == shape2);
    EXPECT_FALSE(shape1 == shape3);
}

// Test with different data types (float)
TEST_F(AlignedBoxCollisionShapeTest, FloatDataType) {
    math::VectorDIM<float, 2> min_f(static_cast<float>(box2d_min(0)), static_cast<float>(box2d_min(1)));
    math::VectorDIM<float, 2> max_f(static_cast<float>(box2d_max(0)), static_cast<float>(box2d_max(1)));
    math::AlignedBox<float, 2> box_f(min_f, max_f);
    
    math::AlignedBoxCollisionShape<float, 2> shape_f(box_f);
    
    // Test basic properties
    EXPECT_TRUE(math::isApproximatelyEqual(min_f, shape_f.collision_box_at_zero().min(), 1e-5f));
    EXPECT_TRUE(math::isApproximatelyEqual(max_f, shape_f.collision_box_at_zero().max(), 1e-5f));
    
    // Test with a position
    math::VectorDIM<float, 2> position_f(1.0f, 2.0f);
    auto bbox_f = shape_f.boundingBox(position_f);
    
    EXPECT_TRUE(math::isApproximatelyEqual(min_f + position_f, bbox_f.min(), 1e-5f));
    EXPECT_TRUE(math::isApproximatelyEqual(max_f + position_f, bbox_f.max(), 1e-5f));
}
