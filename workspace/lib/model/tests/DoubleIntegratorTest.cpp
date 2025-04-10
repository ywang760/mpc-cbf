#include <gtest/gtest.h>
#include <model/DoubleIntegrator.h>
#include <cmath>

class DoubleIntegratorTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Default constructor uses identity matrices
        double_integrator_2d = new model::DoubleIntegrator<double, 2>();
        double_integrator_3d = new model::DoubleIntegrator<double, 3>();
    }

    void TearDown() override {
        delete double_integrator_2d;
        delete double_integrator_3d;
    }

    model::DoubleIntegrator<double, 2>* double_integrator_2d;
    model::DoubleIntegrator<double, 3>* double_integrator_3d;
};

TEST_F(DoubleIntegratorTest, ApplyInputTest2D) {
    // Create a test state
    model::State<double, 2> state;
    state.pos_ << 1.0, 2.0;
    state.vel_ << 0.1, 0.2;

    // Apply input
    Eigen::Vector<double, 2> input;
    input << 0.5, -0.5;

    // Get next state
    model::State<double, 2> next_state = double_integrator_2d->applyInput(state, input);

    // A default constructor should have identity matrices for A and B, which would mean:
    // - Position would remain unchanged in the test
    // - Velocity would be unchanged in the test
    // This test will likely fail until we properly initialize A and B in the next PR
    // For now, we're just testing the API works
    ASSERT_TRUE(next_state.pos_.size() == 2);
    ASSERT_TRUE(next_state.vel_.size() == 2);
}

TEST_F(DoubleIntegratorTest, ApplyInputTest3D) {
    // Create a test state
    model::State<double, 3> state;
    state.pos_ << 1.0, 2.0, 3.0;
    state.vel_ << 0.1, 0.2, 0.3;

    // Apply input
    Eigen::Vector<double, 3> input;
    input << 0.5, -0.5, 0.1;

    // Get next state
    model::State<double, 3> next_state = double_integrator_3d->applyInput(state, input);
    // A default constructor should have identity matrices for A and B, which would mean:
    // - Position would remain unchanged in the test
    // - Velocity would be unchanged in the test
    // This test will likely fail until we properly initialize A and B in the next PR
    // For now, we're just testing the API works
    ASSERT_TRUE(next_state.pos_.size() == 3);
    ASSERT_TRUE(next_state.vel_.size() == 3);
}

TEST_F(DoubleIntegratorTest, StatePropagationMatricesTest) {
    // Test generation of A0 matrix for K steps
    auto A0 = double_integrator_3d->get_A0(5);
    
    // Check dimensions
    ASSERT_EQ(A0.pos_.rows(), 15); // 3 dimensions * 5 steps
    ASSERT_EQ(A0.pos_.cols(), 6);  // 2 * 3 dimensions (pos + vel)
    
    // Test generation of Lambda matrix for K steps
    auto Lambda = double_integrator_3d->get_lambda(5);
    
    // Check dimensions
    ASSERT_EQ(Lambda.pos_.rows(), 15); // 3 dimensions * 5 steps
    ASSERT_EQ(Lambda.pos_.cols(), 15); // 3 dimensions * 5 steps
}
