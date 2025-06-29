#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <random>
#include <utility>
#include <iostream>
#include <spdlog/spdlog.h>
#include <cbf/detail/ConnectivityCBF.h>

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

class InitCBFTest : public ::testing::Test
{

    using ConnectivityCBF = cbf::ConnectivityCBF;

protected:
    // CBF parameters
    double min_dist = 0.5;                                      // minimum distance for connectivity
    double max_dist = 2.0;                                      // maximum distance for connectivity
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

TEST_F(InitCBFTest, TwoRobotRail)
{
    /* Robots on x-axis, 1 m apart (within R_s) */
    Eigen::MatrixXd states(2, 6);
    states << 0.0, 0.0, 0, 0, 0, 0,
        1.0, 0.0, 0, 0, 0, 0;

    int self_idx = 0;
    Eigen::VectorXd x_self = states.row(self_idx).transpose();
    auto res = connectivity_cbf->initConnectivityCBF(states,
                                                    x_self,
                                                    self_idx);

    // Expected results
    Eigen::VectorXd expected_Ac(3);
    expected_Ac << 0.0, 0.0, 0.0;
    double expected_Bc = 0.0;
    checkResult(res, expected_Ac, expected_Bc);
}

// TEST_F(InitCBFTest, RandomSwarmWithinRadius)
// {
//     int N = 1 + static_cast<int>(urand(2.0, 7.9)); // 3 – 8 robots

//     Eigen::MatrixXd states = randomRobotStates(N, /*radius=*/R_s * 0.9);
//     int self_idx = static_cast<int>(urand(0, N - 1));
//     Eigen::VectorXd x_self = states.row(self_idx).transpose();

//     auto res = cbf.initConnectivityCBF(states,
//                                        x_self,
//                                        self_idx);

//     checkResult(res, 6);
// }


// ---------------------------------------------------------------------------

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    spdlog::set_pattern("[%H:%M:%S] [%^%l%$] [%s:%# %!] %v");

    SPDLOG_INFO("Starting Connectivity CBF tests...");
    return RUN_ALL_TESTS();
}