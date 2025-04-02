#include <gtest/gtest.h>
#include <model/DoubleIntegratorXYYaw.h>
#include <cmath>

class DoubleIntegratorXYYawTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create a model with a time step of 0.1 seconds
        model = new model::DoubleIntegratorXYYaw<double>(0.1);
    }

    void TearDown() override {
        delete model;
    }

    model::DoubleIntegratorXYYaw<double>* model;
};

TEST_F(DoubleIntegratorXYYawTest, InitializationTest) {
    // Create an initial state
    model::State<double, 3> state;
    state.pos_ << 1.0, 2.0, 0.5;  // x, y, yaw
    state.vel_ << 0.1, 0.2, 0.3;  // vx, vy, vyaw

    // Create a control input (accelerations)
    Eigen::Vector<double, 3> input;
    input << 0.5, 0.6, 0.1;  // ax, ay, ayaw

    // Apply the input to get the next state
    model::State<double, 3> next_state = model->applyInput(state, input);

    // Test expected values based on the double integrator dynamics
    // For position: p_next = p + v*dt + 0.5*a*dt^2
    // For velocity: v_next = v + a*dt

    // Test X component
    EXPECT_NEAR(next_state.pos_[0], 1.0 + 0.1 * 0.1 + 0.5 * 0.5 * 0.01, 1e-10);
    EXPECT_NEAR(next_state.vel_[0], 0.1 + 0.5 * 0.1, 1e-10);

    // Test Y component
    EXPECT_NEAR(next_state.pos_[1], 2.0 + 0.2 * 0.1 + 0.5 * 0.6 * 0.01, 1e-10);
    EXPECT_NEAR(next_state.vel_[1], 0.2 + 0.6 * 0.1, 1e-10);

    // Test Yaw component
    EXPECT_NEAR(next_state.pos_[2], 0.5 + 0.3 * 0.1 + 0.5 * 0.1 * 0.01, 1e-10);
    EXPECT_NEAR(next_state.vel_[2], 0.3 + 0.1 * 0.1, 1e-10);
}

TEST_F(DoubleIntegratorXYYawTest, PredictionMatricesTest) {
    // Test the state propagator matrices
    int horizon = 10;
    auto A0 = model->get_A0(horizon);
    auto Lambda = model->get_lambda(horizon);
    
    // Check dimensions
    ASSERT_EQ(A0.pos_.rows(), 3 * horizon);
    ASSERT_EQ(A0.pos_.cols(), 6);
    ASSERT_EQ(A0.vel_.rows(), 3 * horizon);
    ASSERT_EQ(A0.vel_.cols(), 6);
    
    ASSERT_EQ(Lambda.pos_.rows(), 3 * horizon);
    ASSERT_EQ(Lambda.pos_.cols(), 3 * horizon);
    ASSERT_EQ(Lambda.vel_.rows(), 3 * horizon);
    ASSERT_EQ(Lambda.vel_.cols(), 3 * horizon);
    
    // Could add more detailed tests for the content of these matrices
    // depending on specific expected behavior
}
