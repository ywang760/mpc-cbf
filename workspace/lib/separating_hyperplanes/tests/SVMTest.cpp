#include <gtest/gtest.h>
#include <separating_hyperplanes/SVM.h>
#include <separating_hyperplanes/Voronoi.h>
#include <random>

namespace separating_hyperplanes {
namespace test {

using VectorDIM2d = VectorDIM<double, 2U>;
using Hyperplane2d = Hyperplane<double, 2U>;

TEST(SVMTest, SimpleLinearlySeparable2D) {
  // Create two linearly separable sets of points
  std::vector<VectorDIM2d> first_set;
  std::vector<VectorDIM2d> second_set;
  
  // First set: points in the bottom-left quadrant
  first_set.push_back(VectorDIM2d(-1.0, -1.0));
  first_set.push_back(VectorDIM2d(-2.0, -1.0));
  first_set.push_back(VectorDIM2d(-1.0, -2.0));
  first_set.push_back(VectorDIM2d(-1.5, -1.5));
  
  // Second set: points in the top-right quadrant
  second_set.push_back(VectorDIM2d(1.0, 1.0));
  second_set.push_back(VectorDIM2d(2.0, 1.0));
  second_set.push_back(VectorDIM2d(1.0, 2.0));
  second_set.push_back(VectorDIM2d(1.5, 1.5));
  
  // Calculate SVM hyperplane
  Hyperplane2d hyperplane = hardMarginSVM<double, 2U>(first_set, second_set);
  
  // Verify that all points in the first set are on the negative side
  for (const auto& point : first_set) {
    double evaluation = hyperplane.normal().dot(point) + hyperplane.offset();
    EXPECT_LT(evaluation, 0.0);
  }
  
  // Verify that all points in the second set are on the positive side
  for (const auto& point : second_set) {
    double evaluation = hyperplane.normal().dot(point) + hyperplane.offset();
    EXPECT_GT(evaluation, 0.0);
  }
  
  // Verify that the margin is correct (distance to closest point from either set is maximized)
  double min_distance = std::numeric_limits<double>::max();
  
  for (const auto& point : first_set) {
    double distance = std::abs(hyperplane.normal().dot(point) + hyperplane.offset()) / 
                     hyperplane.normal().norm();
    min_distance = std::min(min_distance, distance);
  }
  
  for (const auto& point : second_set) {
    double distance = std::abs(hyperplane.normal().dot(point) + hyperplane.offset()) / 
                     hyperplane.normal().norm();
    min_distance = std::min(min_distance, distance);
  }
  
  // For a properly normalized SVM, margin should be close to 1/||w||
  double expected_margin = 1.0 / hyperplane.normal().norm();
  EXPECT_NEAR(min_distance, expected_margin, 1e-10);
}

}  // namespace test
}  // namespace separating_hyperplanes