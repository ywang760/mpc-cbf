#include <ilcplex/ilocplexi.h>
#include <qpcpp/solvers/CPLEX.h>

#include <iostream>
namespace qpcpp {

template <typename T>
CPLEXSolver<T>::CPLEXSolver() : set_parameters_(false) {}

template <typename T>
CPLEXSolver<T>::CPLEXSolver(T optimality_tolerance, T feasibility_tolerance, T time_limit)
    : set_parameters_(true),
      optimality_tolerance_(optimality_tolerance),
      feasibility_tolerance_(feasibility_tolerance),
      time_limit_(time_limit) {}

template <typename T>
CPLEXSolver<T>::~CPLEXSolver() {}

template <typename T>
void CPLEXSolver<T>::setSolution(const IloEnv& cplex_env, const IloCplex& cplex_solver,
                                 const IloNumVarArray& cplex_variables,
                                 std::vector<Variable*>& cplex_indices_problem_variables) {
    IloNumArray cplex_solution(cplex_env);
    cplex_solver.getValues(cplex_solution, cplex_variables);

    for (std::size_t variable_idx = 0; variable_idx < cplex_indices_problem_variables.size();
         ++variable_idx) {
        Variable* variable_ptr = cplex_indices_problem_variables[variable_idx];
        variable_ptr->set_solution_value(cplex_solution[variable_idx]);
    }
}

template <typename T>
SolveStatus CPLEXSolver<T>::solve(Problem& problem) {
    using Variable = qpcpp::Variable<T>;
    using LinearConstraint = qpcpp::LinearConstraint<T>;
    using CostFunction = qpcpp::CostFunction<T>;

    IloEnv cplex_env;
    cplex_env.setOut(cplex_env.getNullStream());

    IloModel cplex_model(cplex_env);
    typename std::forward_list<Variable>& problem_variables =
        problem.mutable_variables(); // variables are mutable because they will
    // be updated by the algorithm once the
    // solution is found
    const typename std::forward_list<LinearConstraint>& problem_linear_constraints =
        problem.linear_constraints();
    const CostFunction* cost_function = problem.cost_function();

    // create variables
    IloNumVarArray cplex_variables(cplex_env);
    // [cplex var idx] => const Variable* in problem
    std::vector<Variable*> cplex_indices_problem_variables(problem.numVariables());
    typename std::forward_list<Variable>::iterator problem_variable_it = problem_variables.begin();
    for (std::size_t variable_idx = 0; variable_idx < problem.numVariables();
         ++variable_idx, ++problem_variable_it) {
        Variable& problem_variable = *problem_variable_it;

        IloNumVar cplex_variable(cplex_env, problem_variable.min(), problem_variable.max(),
                                 ILOFLOAT);
        cplex_variables.add(cplex_variable);
        cplex_model.add(cplex_variable);
        cplex_indices_problem_variables[variable_idx] = &problem_variable;
    }

    // add linear constraints
    for (const LinearConstraint& problem_linear_constraint : problem_linear_constraints) {
        IloExpr cplex_expr(cplex_env);
        for (std::size_t variable_idx = 0; variable_idx < problem.numVariables(); ++variable_idx) {
            const Variable* problem_variable_ptr = cplex_indices_problem_variables[variable_idx];
            T coefficient = problem_linear_constraint.getCoefficient(problem_variable_ptr);
            cplex_expr += cplex_variables[variable_idx] * (coefficient);
        }

        IloRange cplex_range(cplex_env, problem_linear_constraint.min(), cplex_expr,
                             problem_linear_constraint.max());
        cplex_model.add(cplex_range);
    }

    // set cost functions
    IloExpr cplex_quadratic_cost(cplex_env);
    IloExpr cplex_linear_cost(cplex_env);
    for (std::size_t first_variable_idx = 0; first_variable_idx < problem.numVariables();
         ++first_variable_idx) {
        const Variable* first_variable_ptr = cplex_indices_problem_variables[first_variable_idx];

        T coefficient = cost_function->getLinearCoefficient(first_variable_ptr);
        cplex_linear_cost += cplex_variables[first_variable_idx] * (coefficient);

        for (std::size_t second_variable_idx = first_variable_idx;
             second_variable_idx < problem.numVariables(); ++second_variable_idx) {
            const Variable* second_variable_ptr =
                cplex_indices_problem_variables[second_variable_idx];
            T coefficient =
                cost_function->getQuadraticCoefficient(first_variable_ptr, second_variable_ptr);

            cplex_quadratic_cost += cplex_variables[first_variable_idx] *
                                    cplex_variables[second_variable_idx] * (coefficient);
        }
    }

    IloObjective cplex_objective(
        cplex_env, cplex_quadratic_cost + cplex_linear_cost + cost_function->constant(),
        IloObjective::Minimize);
    cplex_model.add(cplex_objective);

    IloCplex cplex_solver(cplex_model);
    cplex_solver.setOut(cplex_env.getNullStream());
    cplex_solver.setWarning(cplex_env.getNullStream());
    cplex_solver.setError(cplex_env.getNullStream());

    cplex_solver.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Primal);
    cplex_solver.setParam(IloCplex::Param::OptimalityTarget, CPX_OPTIMALITYTARGET_OPTIMALCONVEX);
    // cplex_solver.setParam(IloCplex::Param::Parallel,
    //                       CPX_PARALLEL_DETERMINISTIC);
    cplex_solver.setParam(IloCplex::Param::Threads, 1);

    if (set_parameters_) {
        cplex_solver.setParam(IloCplex::Param::Simplex::Tolerances::Feasibility,
                              feasibility_tolerance_);
        cplex_solver.setParam(IloCplex::Param::Simplex::Tolerances::Optimality,
                              optimality_tolerance_);
        cplex_solver.setParam(IloCplex::Param::TimeLimit, time_limit_);
    }

    // std::cout << "cplex before cvx solve" << std::endl;
    try {
        cplex_solver.solve();
        // std::cout << "cplex after cvx solve" << std::endl;
    } catch (const IloCplex::Exception& exp) {
        // std::cout << "cplex cvx solve failed" << std::endl;
        cplex_solver.setParam(IloCplex::Param::OptimalityTarget, CPX_OPTIMALITYTARGET_FIRSTORDER);
        cplex_solver.setParam(IloCplex::Param::RootAlgorithm, IloCplex::AutoAlg);
        cplex_solver.solve();
        // std::cout << "cplex after first order solve" << std::endl;
    }

    IloAlgorithm::Status cplex_solve_status = cplex_solver.getStatus();
    SolveStatus return_status = SolveStatus::UNKNOWN;

    switch (cplex_solve_status) {
    case IloAlgorithm::Status::Unknown: {
        return_status = SolveStatus::UNKNOWN;
        break;
    }
    case IloAlgorithm::Status::Feasible: {
        setSolution(cplex_env, cplex_solver, cplex_variables, cplex_indices_problem_variables);
        return_status = SolveStatus::FEASIBLE;
        break;
    }
    case IloAlgorithm::Status::Optimal: {
        setSolution(cplex_env, cplex_solver, cplex_variables, cplex_indices_problem_variables);
        return_status = SolveStatus::OPTIMAL;
        break;
    }
    case IloAlgorithm::Status::Infeasible: {
        return_status = SolveStatus::INFEASIBLE;
        break;
    }
    case IloAlgorithm::Status::Unbounded: {
        return_status = SolveStatus::UNBOUNDED;
        break;
    }
    case IloAlgorithm::Status::InfeasibleOrUnbounded: {
        return_status = SolveStatus::INFEASIBLEORUNBOUNDED;
        break;
    }
    case IloAlgorithm::Status::Error: {
        return_status = SolveStatus::ERROR;
        break;
    }
    }

    cplex_env.end();
    return return_status;
}

template class CPLEXSolver<float>;
template class CPLEXSolver<double>;
} // namespace qpcpp