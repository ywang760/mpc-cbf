#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <random>
#include <utility>
#include <iostream>
#include <spdlog/spdlog.h>
#include <cbf/detail/ConnectivityCBF.h>
#include <common/logging.hpp>

// Convenience: check vector size then its entries
#define EXPECT_SIZE(vec, n) EXPECT_EQ((vec).size(), (n))

// ---------- Helpers ---------------------------------------------------------

// Generate a random double in [a,b]
static double urand(double a, double b)
{
    static std::mt19937_64 rng{std::random_device{}()};
    std::uniform_real_distribution<double> dist(a, b);
    return dist(rng);
}

// Build a random N×6 state matrix where position components lie in a disk
static Eigen::MatrixXd randomRobotStates(int N, double radius)
{
    Eigen::MatrixXd X(N, 6);
    X.setZero();
    for (int i = 0; i < N; ++i)
    {
        // disc sampling
        double theta = urand(0.0, 2.0 * M_PI);
        double r = urand(0.0, radius);
        X(i, 0) = r * std::cos(theta); // x
        X(i, 1) = r * std::sin(theta); // y
        // z = 0, v = 0, etc.  (modify if you use 3-D or non-zero velocities)
    }
    return X;
}

// ---------- Test fixture ----------------------------------------------------

class InitConnCBFTest : public ::testing::Test
{

    using ConnectivityCBF = cbf::ConnectivityCBF;

protected:
    // CBF parameters
    double min_dist = 0.8;                                      // minimum distance for safety (not used)
    double max_dist = 3.0;                                      // maximum distance for connectivity
    Eigen::VectorXd v_min = Eigen::VectorXd::Constant(3, -1.0); // min velocity limits (-1.0 for all 3 dims)
    Eigen::VectorXd v_max = Eigen::VectorXd::Constant(3, 1.0);  // max velocity limits (1.0 for all 3 dims)

    std::shared_ptr<ConnectivityCBF> connectivity_cbf = std::make_shared<ConnectivityCBF>(min_dist, max_dist, v_min, v_max);

    void checkResult(const std::pair<Eigen::VectorXd, double> &res,
                    const Eigen::VectorXd &expected_Ac,
                    double expected_Bc) const
    {
        const Eigen::VectorXd &Ac = res.first;
        double Bc = res.second;

        // Check for Ac
        EXPECT_SIZE(Ac, 3);
        EXPECT_TRUE(Ac.allFinite());
        EXPECT_TRUE(Ac.isApprox(expected_Ac));

        // Check for Bc
        EXPECT_TRUE(std::isfinite(Bc));
        EXPECT_GE(Bc, 0.0); // Bc should be non-negative
        EXPECT_DOUBLE_EQ(Bc, expected_Bc);

        SPDLOG_INFO("Result: Ac = {}, Bc = {}", Ac.transpose(), Bc);
    }
};

// ---------- TEST CASES ------------------------------------------------------

// TEST_F(InitConnCBFTest, TwoRobotRail)
// {
//     /* Robots on x-axis, 1 m apart (within R_s) */
//     Eigen::MatrixXd robot_states(2, 6);
//     robot_states << 0.0, 0.0, 0, 0, 0, 0,
//         1.0, 0.0, 0, 0, 0, 0;

//     int self_idx = 0;
//     Eigen::VectorXd state = robot_states.row(self_idx).transpose();
//     auto res = connectivity_cbf->initConnCBF(robot_states,
//                                              state,
//                                              self_idx);

//     // Expected results
//     Eigen::VectorXd expected_Ac(3);
//     expected_Ac << 0.0, 0.0, 0.0;
//     double expected_Bc = 0.0;
//     // checkResult(res, expected_Ac, expected_Bc);
// }

TEST_F(InitConnCBFTest, Misc)
{
    auto logger = common::initializeLogging();
    /* Robots on x-axis, 1 m apart (within R_s) */
    Eigen::MatrixXd robot_states(3, 6);
    robot_states << 1.0, 2.0, 0, 0, 0, 0,
        1.0, 4.0, 0, 0, 0, 0,
        1.0, 6.0, 0, 0, 0, 0;

    int self_idx = 0;
    Eigen::VectorXd state = robot_states.row(self_idx).transpose();

    const double epsilon = 0.1; // minimum eigenvalue threshold

    const auto robot_positions = robot_states.leftCols(2); // Extract only the position columns (x, y)
    auto [lambda2_val, eigenvec] = connectivity_cbf->getLambda2(robot_positions);
    double h = lambda2_val - epsilon; // barrier function: h = λ₂ - λ₂_min
    auto [Ac_sym, Bc_sym] = connectivity_cbf->initConnCBF(robot_states, 0);
    Eigen::VectorXd Ac = connectivity_cbf->getConnConstraints(state, robot_states, eigenvec);
    double Bc = connectivity_cbf->getConnBound(state, robot_states, eigenvec, h);
}

TEST_F(InitConnCBFTest, Misc2)
{
    auto logger = common::initializeLogging();
    /* Robots on x-axis, 1 m apart (within R_s) */
    Eigen::MatrixXd robot_states(3, 6);
    robot_states << 0.212, 1.592, 0, -0.293, -0.21, 0.0,
        1.01, 4.20, 0, -1.2, 0.12, 0,
        -1.0, -0.02, 0, -0.2, 0.16, 0;

    int self_idx = 0;
    Eigen::VectorXd state = robot_states.row(self_idx).transpose();
    const double epsilon = 0.1;                            // minimum eigenvalue threshold
    const auto robot_positions = robot_states.leftCols(2); // Extract only the position columns (x, y)
    auto [lambda2_val, eigenvec] = connectivity_cbf->getLambda2(robot_positions);
    double h = lambda2_val - epsilon; // barrier function: h = λ₂ - λ₂_min
    auto [Ac_sym, Bc_sym] = connectivity_cbf->initConnCBF(robot_states, 0);
    Eigen::VectorXd Ac = connectivity_cbf->getConnConstraints(state, robot_states, eigenvec);
    double Bc = connectivity_cbf->getConnBound(state, robot_states, eigenvec, h);
}