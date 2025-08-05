//
// Created by lishuo on 2/18/24.
//

#include <qpcpp/solvers/Solver.h>

namespace qpcpp {
std::string SolveStatusToStr(SolveStatus status) {
    switch (status) {
    case SolveStatus::ERROR: {
        return "ERROR";
    }
    case SolveStatus::FEASIBLE: {
        return "FEASIBLE";
    }
    case SolveStatus::INFEASIBLE: {
        return "INFEASIBLE";
    }
    case SolveStatus::INFEASIBLEORUNBOUNDED: {
        return "INFEASIBLEORUNBOUNDED";
    }
    case SolveStatus::OPTIMAL: {
        return "OPTIMAL";
    }
    case SolveStatus::UNBOUNDED: {
        return "UNBOUNDED";
    }
    case SolveStatus::UNKNOWN: {
        return "UNKNOWN";
    }
    }

    return "";
}

} // namespace qpcpp