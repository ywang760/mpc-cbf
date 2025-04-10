#include <gtest/gtest.h>
#include "math/Combinatorics.h"
#include <limits>
#include <stdexcept>

class CombinatoricsTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Setup code before each test if needed
    }

    void TearDown() override {
        // Cleanup code after each test if needed
    }
};

TEST_F(CombinatoricsTest, TestFactorial) {
    // Test basic factorial values
    EXPECT_EQ(math::fac(0), 1);
    EXPECT_EQ(math::fac(1), 1);
    EXPECT_EQ(math::fac(2), 2);
    EXPECT_EQ(math::fac(3), 6);
    EXPECT_EQ(math::fac(4), 24);
    EXPECT_EQ(math::fac(5), 120);
    EXPECT_EQ(math::fac(10), 3628800);
    EXPECT_EQ(math::fac(20), 2432902008176640000ULL);

    // Test overflow error
    EXPECT_THROW(math::fac(21), std::runtime_error);
}

TEST_F(CombinatoricsTest, TestCombination) {
    // Test basic combination values
    EXPECT_EQ(math::comb(5, 0), 1);
    EXPECT_EQ(math::comb(5, 1), 5);
    EXPECT_EQ(math::comb(5, 2), 10);
    EXPECT_EQ(math::comb(5, 3), 10);
    EXPECT_EQ(math::comb(5, 4), 5);
    EXPECT_EQ(math::comb(5, 5), 1);

    // Test edge cases
    EXPECT_EQ(math::comb(0, 0), 1);
    EXPECT_EQ(math::comb(10, 0), 1);
    EXPECT_EQ(math::comb(10, 10), 1);

    // Test k > n case
    EXPECT_EQ(math::comb(5, 6), 0);

    // Test larger values
    EXPECT_EQ(math::comb(20, 10), 184756);
}

TEST_F(CombinatoricsTest, TestPermutation) {
    // Test basic permutation values
    EXPECT_EQ(math::perm(5, 0), 1);
    EXPECT_EQ(math::perm(5, 1), 5);
    EXPECT_EQ(math::perm(5, 2), 20);
    EXPECT_EQ(math::perm(5, 3), 60);
    EXPECT_EQ(math::perm(5, 4), 120);
    EXPECT_EQ(math::perm(5, 5), 120);

    // Test edge cases
    EXPECT_EQ(math::perm(0, 0), 1);
    EXPECT_EQ(math::perm(10, 0), 1);

    // Test k > n case
    EXPECT_EQ(math::perm(5, 6), 0);

    // Test larger values
    EXPECT_EQ(math::perm(10, 3), 720);
}

TEST_F(CombinatoricsTest, TestPow) {
    // Test integer powers
    EXPECT_DOUBLE_EQ((math::pow<double, int>(2.0, 3)), 8.0);
    EXPECT_DOUBLE_EQ((math::pow<double, int>(2.0, -2)), 0.25);
    EXPECT_FLOAT_EQ((math::pow<float, int>(3.0f, 2)), 9.0f);

    // Test unsigned powers
    EXPECT_DOUBLE_EQ((math::pow<double, uint64_t>(2.0, 10)), 1024.0);
    EXPECT_FLOAT_EQ((math::pow<float, uint64_t>(1.5f, 3)), 3.375f);

    // Test edge cases
    EXPECT_DOUBLE_EQ((math::pow<double, int>(1.0, 1000000)), 1.0);
    EXPECT_DOUBLE_EQ((math::pow<double, int>(0.0, 5)), 0.0);
    EXPECT_DOUBLE_EQ((math::pow<double, int>(5.0, 0)), 1.0);

    // Special case: 0^0 = 1
    EXPECT_DOUBLE_EQ((math::pow<double, int>(0.0, 0)), 1.0);

    // Test error case
    EXPECT_THROW((math::pow<double, int>(0.0, -1)), std::runtime_error);
}
