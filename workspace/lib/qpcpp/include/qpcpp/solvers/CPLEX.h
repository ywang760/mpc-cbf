//
// Created by lishuo on 4/7/24.
//

#ifndef QPCPP_SOLVERS_CPLEX_H
#define QPCPP_SOLVERS_CPLEX_H

#include <ilcplex/ilocplex.h>
#include <qpcpp/solvers/Solver.h>

namespace qpcpp {
template <typename T>
class CPLEXSolver : public Solver<T> {
  public:
    using Base = Solver<T>;
    using Problem = typename Base::Problem;
    using Variable = qpcpp::Variable<T>;

    CPLEXSolver();

    CPLEXSolver(T optimality_tolerance, T feasibility_tolerance, T time_limit);

    virtual ~CPLEXSolver();
    virtual SolveStatus solve(Problem& problem) override;

  private:
    // Set the solution found by cplex_solver within the context of cplex_env
    // from cplex_variables to Variable*s in cplex_indices_problem_variables.
    // cplex_indices_problem_variables is a map from variable indices in cplex
    // to Variable*s.
    void setSolution(const IloEnv& cplex_env, const IloCplex& cplex_solver,
                     const IloNumVarArray& cplex_variables,
                     std::vector<Variable*>& cplex_indices_problem_variables);

    bool set_parameters_;
    T optimality_tolerance_;
    T feasibility_tolerance_;
    T time_limit_;
};

} // namespace qpcpp

#endif //QPCPP_SOLVERS_CPLEX_H
