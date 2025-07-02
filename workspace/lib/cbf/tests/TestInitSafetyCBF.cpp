#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <random>
#include <utility>
#include <iostream>
#include <spdlog/spdlog.h>
#include <cbf/detail/ConnectivityCBF.h>
#include <math/Types.h>

// Convenience: check vector size then its entries
#define EXPECT_SIZE(vec, n) EXPECT_EQ((vec).size(), (n))

// ---------- Test fixture ----------------------------------------------------

class InitSafetyCBFTest : public ::testing::Test
{

    using ConnectivityCBF = cbf::ConnectivityCBF;
    

protected:
    // CBF parameters
    double min_dist = 0.8;                                      // minimum distance for connectivity
    double max_dist = 2.0;                                      // maximum distance for connectivity (doens't matter for safety CBF)
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
        EXPECT_NEAR(Bc, expected_Bc, 1e-6); // Allow

        SPDLOG_INFO("Result: Ac = {}, Bc = {}", Ac.transpose(), Bc);
    }
};

// ---------- TEST CASES ------------------------------------------------------

TEST_F(InitSafetyCBFTest, TwoRobotInSafeRegion)
{
    Eigen::VectorXd state(6);
    state << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0; // ego robot pos, vel

    Eigen::VectorXd other_state(2);
    other_state << 1.0, 0.0; // other robot pos

    // Get safety constraints
    auto Ac = connectivity_cbf->getSafetyConstraints(state, other_state);
    auto Bc = connectivity_cbf->getSafetyBound(state, other_state);

    auto res = std::pair<Eigen::VectorXd, double>(Ac, Bc);
    
    // Since distance is greater than min_dist, we expect that Bc to be positive
    // The constraint becomes Ac * u + Bc >= 0
    // Any u will keep the system safe, so the constraint is inactive
    EXPECT_GT(Bc, 0.0);

    Eigen::VectorXd expected_Ac(3);
    expected_Ac << -2.0, 0.0, 0.0;
    double expected_Bc = 0.06347497291775989;
    checkResult(res, expected_Ac, expected_Bc);
}

TEST_F(InitSafetyCBFTest, TwoRobotInSafeRegionWithHugeVelocity)
{
    Eigen::VectorXd state(6);
    state << 0.0, 0.0, 0.0, 100.0, 100.0, 0.0; // ego robot pos, vel

    Eigen::VectorXd other_state(2);
    other_state << 1.0, 0.0; // other robot pos

    // Get safety constraints
    auto Ac = connectivity_cbf->getSafetyConstraints(state, other_state);
    auto Bc = connectivity_cbf->getSafetyBound(state, other_state);

    auto res = std::pair<Eigen::VectorXd, double>(Ac, Bc);
    
    // In this case, since there's huge velocity, we expect that Bc to be negative
    // The CBF is active to ensure safety
    EXPECT_LT(Bc, 0.0);

    Eigen::VectorXd expected_Ac(3);
    expected_Ac << -2.0, 0.0, 0.0;
    double expected_Bc = -39820583.995200224;
    checkResult(res, expected_Ac, expected_Bc);
}

TEST_F(InitSafetyCBFTest, TwoRobotOnSafetyBound)
{
    Eigen::VectorXd state(6);
    state << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0; // ego robot pos, vel

    Eigen::VectorXd other_state(2);
    other_state << 0.8, 0.0; // other robot pos (on safety bound)

    // Get safety constraints
    auto Ac = connectivity_cbf->getSafetyConstraints(state, other_state);
    auto Bc = connectivity_cbf->getSafetyBound(state, other_state);

    auto res = std::pair<Eigen::VectorXd, double>(Ac, Bc);
    
    // Since distance is equal to min_dist, we expect that Bc to be zero.
    EXPECT_EQ(Bc, 0.0);

    Eigen::VectorXd expected_Ac(3);
    expected_Ac << -1.6, 0.0, 0.0;
    double expected_Bc = 0;
    checkResult(res, expected_Ac, expected_Bc);
}

TEST_F(InitSafetyCBFTest, TwoRobotInUnsafeRegion)
{
    Eigen::VectorXd state(6);
    state << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0; // ego robot pos, vel

    Eigen::VectorXd other_state(2);
    other_state << 0.5, 0.0; // other robot pos (unsafe)

    // Get safety constraints
    auto Ac = connectivity_cbf->getSafetyConstraints(state, other_state);
    auto Bc = connectivity_cbf->getSafetyBound(state, other_state);

    auto res = std::pair<Eigen::VectorXd, double>(Ac, Bc);
    
    // Since distance is less than min_dist, we expect that Bc to be negative
    // The constraint Ac * u + Bc >= 0 should be satisfied
    EXPECT_LT(Bc, 0.0);

    Eigen::VectorXd expected_Ac(3);
    expected_Ac << -1.0, 0.0, 0.0;
    double expected_Bc = -0.13045522572422458;
    checkResult(res, expected_Ac, expected_Bc);
}
