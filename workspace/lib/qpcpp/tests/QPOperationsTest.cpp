#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <qpcpp/QPOperations.h>

namespace qpcpp {
namespace {

class QPOperationsTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Common setup for each test
    }

    void TearDown() override {
        // Common cleanup after each test
    }
};

TEST_F(QPOperationsTest, CostAdditionInitialization) {
    using Matrix = math::Matrix<double>;
    using Vector = math::Vector<double>;
    
    Matrix Q = Matrix::Identity(2, 2);
    Vector c = Vector::Ones(2);
    double constant = 5.0;
    
    QPOperations<double>::CostAddition cost(Q, c, constant);
    
    EXPECT_EQ(cost.quadratic_term(), Q);
    EXPECT_EQ(cost.linear_term(), c);
    EXPECT_EQ(cost.constant(), constant);
}

TEST_F(QPOperationsTest, LinearConstraintInitialization) {
    using Row = math::Row<double>;
    
    Row coefficients = Row::Ones(3);
    double lower_bound = -1.0;
    double upper_bound = 2.0;
    
    QPOperations<double>::LinearConstraint constraint(coefficients, lower_bound, upper_bound);
    
    EXPECT_EQ(constraint.coefficients(), coefficients);
    EXPECT_EQ(constraint.lower_bound(), lower_bound);
    EXPECT_EQ(constraint.upper_bound(), upper_bound);
}

TEST_F(QPOperationsTest, DecisionVariableBoundsInitialization) {
    using Vector = math::Vector<double>;
    
    Vector lower_bounds = Vector::Constant(3, -1.0);
    Vector upper_bounds = Vector::Constant(3, 1.0);
    
    QPOperations<double>::DecisionVariableBounds bounds(lower_bounds, upper_bounds);
    
    EXPECT_EQ(bounds.lower_bounds(), lower_bounds);
    EXPECT_EQ(bounds.upper_bounds(), upper_bounds);
}

TEST_F(QPOperationsTest, CostAdditionWithFloats) {
    using Matrix = math::Matrix<float>;
    using Vector = math::Vector<float>;
    
    Matrix Q = Matrix::Identity(2, 2);
    Vector c = Vector::Ones(2);
    float constant = 5.0f;
    
    QPOperations<float>::CostAddition cost(Q, c, constant);
    
    EXPECT_EQ(cost.quadratic_term(), Q);
    EXPECT_EQ(cost.linear_term(), c);
    EXPECT_EQ(cost.constant(), constant);
}

TEST_F(QPOperationsTest, LinearConstraintWithFloats) {
    using Row = math::Row<float>;
    
    Row coefficients = Row::Ones(3);
    float lower_bound = -1.0f;
    float upper_bound = 2.0f;
    
    QPOperations<float>::LinearConstraint constraint(coefficients, lower_bound, upper_bound);
    
    EXPECT_EQ(constraint.coefficients(), coefficients);
    EXPECT_EQ(constraint.lower_bound(), lower_bound);
    EXPECT_EQ(constraint.upper_bound(), upper_bound);
}

TEST_F(QPOperationsTest, DecisionVariableBoundsWithFloats) {
    using Vector = math::Vector<float>;
    
    Vector lower_bounds = Vector::Constant(3, -1.0f);
    Vector upper_bounds = Vector::Constant(3, 1.0f);
    
    QPOperations<float>::DecisionVariableBounds bounds(lower_bounds, upper_bounds);
    
    EXPECT_EQ(bounds.lower_bounds(), lower_bounds);
    EXPECT_EQ(bounds.upper_bounds(), upper_bounds);
}

}  // namespace
}  // namespace qpcpp
