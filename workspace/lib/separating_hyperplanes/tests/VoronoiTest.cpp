#include <gtest/gtest.h>
#include <separating_hyperplanes/Voronoi.h>

namespace separating_hyperplanes {
namespace test {

using VectorDIM2d = VectorDIM<double, 2U>;
using Hyperplane2d = Hyperplane<double, 2U>;

TEST(VoronoiTest, ComputeVoronoiHyperplane2D) {
  // Test case with 2D points
  VectorDIM2d p1(1.0, 1.0);
  VectorDIM2d p2(4.0, 5.0);
  
  // Compute the Voronoi hyperplane
  Hyperplane2d hyperplane = separating_hyperplanes::voronoi<double, 2U>(p1, p2);
  
  // The normal should point from p1 to p2
  VectorDIM2d expected_direction = p2 - p1;
  expected_direction.normalize();
  
  // Verify hyperplane normal direction
  for (int i = 0; i < 2; i++) {
    EXPECT_NEAR(expected_direction(i), hyperplane.normal()(i), 1e-10);
  }
  
  // Hyperplane should pass through the midpoint
  VectorDIM2d midpoint = (p1 + p2) / 2.0;
  double evaluation = hyperplane.normal().dot(midpoint) + hyperplane.offset();
  EXPECT_NEAR(0.0, evaluation, 1e-10);
  
  // p1 should be on the negative side
  EXPECT_LT(hyperplane.normal().dot(p1) + hyperplane.offset(), 0.0);
  
  // p2 should be on the positive side
  EXPECT_GT(hyperplane.normal().dot(p2) + hyperplane.offset(), 0.0);
  
  // Distance from p1 to hyperplane should equal distance from p2 to hyperplane
  double p1_to_hyperplane = std::abs(hyperplane.normal().dot(p1) + hyperplane.offset()) / 
                            hyperplane.normal().norm();
  double p2_to_hyperplane = std::abs(hyperplane.normal().dot(p2) + hyperplane.offset()) / 
                            hyperplane.normal().norm();
  EXPECT_NEAR(p1_to_hyperplane, p2_to_hyperplane, 1e-10);
}

TEST(VoronoiTest, EquidistanceProperty) {
  // Create random points 
  VectorDIM2d p1(2.5, -1.0);
  VectorDIM2d p2(-3.0, 4.0);
  
  // Compute Voronoi hyperplane
  Hyperplane2d hyperplane = separating_hyperplanes::voronoi<double, 2U>(p1, p2);
  
  // Generate test points
  const int num_test_points = 10;
  for (int i = 0; i < num_test_points; i++) {
    // Generate random point on the hyperplane
    double t = -5.0 + 10.0 * i / num_test_points; // Parameter along hyperplane
    VectorDIM2d normal = hyperplane.normal();
    VectorDIM2d perpendicular(-normal(1), normal(0)); // Perpendicular direction
    
    // Point on the hyperplane: solve normal·x + offset = 0
    VectorDIM2d point_on_plane = (perpendicular * t) - (normal * hyperplane.offset() / normal.squaredNorm());
    
    // Verify point is on hyperplane
    EXPECT_NEAR(0.0, hyperplane.normal().dot(point_on_plane) + hyperplane.offset(), 1e-10);
    
    // Distance to p1 and p2 should be equal
    double dist_to_p1 = (point_on_plane - p1).norm();
    double dist_to_p2 = (point_on_plane - p2).norm();
    EXPECT_NEAR(dist_to_p1, dist_to_p2, 1e-10);
  }
}

}  // namespace test
}  // namespace separating_hyperplanes
