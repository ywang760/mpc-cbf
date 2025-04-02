#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <qpcpp/Problem.h>

namespace qpcpp {
namespace {

class ProblemTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Common setup for each test
    }

    void TearDown() override {
        // Common cleanup after each test
    }
};

TEST_F(ProblemTest, InitialProblemStateIsEmpty) {
    Problem<double> problem;
    EXPECT_EQ(problem.numVariables(), 0);
    EXPECT_EQ(problem.numLinearConstraints(), 0);
}

TEST_F(ProblemTest, AddingVariablesIncreasesCount) {
    Problem<double> problem;
    problem.addVariable(-1.0, 1.0);
    EXPECT_EQ(problem.numVariables(), 1);
    problem.addVariable(-2.0, 2.0);
    EXPECT_EQ(problem.numVariables(), 2);
}

TEST_F(ProblemTest, VariableMinMaxLimits) {
    Problem<double> problem;
    auto* var = problem.addVariable(-1.5, 2.5);
    EXPECT_EQ(var->min(), -1.5);
    EXPECT_EQ(var->max(), 2.5);
    
    var->set_min(-3.0);
    var->set_max(4.0);
    EXPECT_EQ(var->min(), -3.0);
    EXPECT_EQ(var->max(), 4.0);
}

TEST_F(ProblemTest, HasVariableTest) {
    Problem<double> problem;
    auto* var1 = problem.addVariable();
    auto* var2 = problem.addVariable();
    
    EXPECT_TRUE(problem.hasVariable(var1));
    EXPECT_TRUE(problem.hasVariable(var2));
    
    Problem<double> other_problem;
    auto* other_var = other_problem.addVariable();
    EXPECT_FALSE(problem.hasVariable(other_var));
}

TEST_F(ProblemTest, VariableSolutionValueTest) {
    Problem<double> problem;
    auto* var = problem.addVariable(-1.0, 1.0);
    
    var->set_solution_value(0.5);
    EXPECT_EQ(var->solution_value(), 0.5);
}

TEST_F(ProblemTest, AddingConstraintsIncreasesCount) {
    Problem<double> problem;
    problem.addLinearConstraint(-1.0, 1.0);
    EXPECT_EQ(problem.numLinearConstraints(), 1);
    problem.addLinearConstraint(-2.0, 2.0);
    EXPECT_EQ(problem.numLinearConstraints(), 2);
}

TEST_F(ProblemTest, ConstraintCoefficients) {
    Problem<double> problem;
    auto* var1 = problem.addVariable();
    auto* var2 = problem.addVariable();
    auto* constraint = problem.addLinearConstraint(0.0, 5.0);
    
    constraint->setCoefficient(var1, 2.0);
    constraint->setCoefficient(var2, 3.0);
    
    EXPECT_EQ(constraint->getCoefficient(var1), 2.0);
    EXPECT_EQ(constraint->getCoefficient(var2), 3.0);
}

TEST_F(ProblemTest, ClearingConstraints) {
    Problem<double> problem;
    problem.addLinearConstraint();
    problem.addLinearConstraint();
    EXPECT_EQ(problem.numLinearConstraints(), 2);
    
    problem.clearLinearConstraints();
    EXPECT_EQ(problem.numLinearConstraints(), 0);
}

TEST_F(ProblemTest, CostFunctionLinearTerms) {
    Problem<double> problem;
    auto* var1 = problem.addVariable();
    auto* var2 = problem.addVariable();
    auto* cost = problem.cost_function();
    
    cost->addLinearTerm(var1, 2.5);
    cost->addLinearTerm(var2, 3.5);
    
    EXPECT_EQ(cost->getLinearCoefficient(var1), 2.5);
    EXPECT_EQ(cost->getLinearCoefficient(var2), 3.5);
}

TEST_F(ProblemTest, CostFunctionQuadraticTerms) {
    Problem<double> problem;
    auto* var1 = problem.addVariable();
    auto* var2 = problem.addVariable();
    auto* cost = problem.cost_function();
    
    cost->addQuadraticTerm(var1, var2, 2.5);
    
    EXPECT_EQ(cost->getQuadraticCoefficient(var1, var2), 2.5);
    // Check symmetry
    EXPECT_EQ(cost->getQuadraticCoefficient(var2, var1), 2.5);
}

TEST_F(ProblemTest, ResetProblem) {
    Problem<double> problem;
    auto* var = problem.addVariable();
    problem.addLinearConstraint();
    auto* cost = problem.cost_function();
    cost->addLinearTerm(var, 1.0);
    
    EXPECT_EQ(problem.numLinearConstraints(), 1);
    EXPECT_EQ(cost->getLinearCoefficient(var), 1.0);
    
    problem.resetProblem();
    
    EXPECT_EQ(problem.numLinearConstraints(), 0);
    EXPECT_EQ(cost->getLinearCoefficient(var), 0.0);
}

}  // namespace
}  // namespace qpcpp
