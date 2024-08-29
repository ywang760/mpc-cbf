//
// Created by lishuo on 2/18/24.
//

#ifndef QPCPP_SOLVERS_SOLVER_H
#define QPCPP_SOLVERS_SOLVER_H

#include <qpcpp/Problem.h>

namespace qpcpp {
    // return type of Solver::solve function. return types are inspired from
    // CPLEX.
    enum class SolveStatus {
        OPTIMAL,
        FEASIBLE,
        UNBOUNDED,
        INFEASIBLE,
        ERROR,
        UNKNOWN,
        INFEASIBLEORUNBOUNDED
    };

    std::string SolveStatusToStr(SolveStatus status);

    // abstract base class for qp solvers
    template <typename T>
    class Solver {
    public:
        using Problem = qpcpp::Problem<T>;

        virtual ~Solver() = default;

        // solves the problem and sets the solution to the variables in the
        // problem. solution is only set if return status is ok and SolveStatus
        // \in {OPTIMAL, FEASIBLE}
        virtual SolveStatus solve(Problem& problem) = 0;
    };
} // qpcpp

#endif //QPCPP_SOLVERS_SOLVER_H
