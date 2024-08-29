#ifndef QPCPP_PROBLEM_H
#define QPCPP_PROBLEM_H

#include <vector>
#include <forward_list>
#include <unordered_set>
#include <unordered_map>

namespace qpcpp {

    template <typename T>
    class Problem;

    template <typename T>
    class Variable {
    public:
        using Problem = qpcpp::Problem<T>;

        const Problem& context_problem() const;
        T min() const;
        T max() const;

        void set_min(T new_min);
        void set_max(T new_max);

        void set_solution_value(T solution_value);
        T solution_value() const;

    private:
        // variables can't be created from outside. only Problem instances can
        // create them.
        Variable(const Problem& context_problem, T min, T max);

        const Problem& context_problem_;
        T min_;
        T max_;
        T solution_value_;

        friend class qpcpp::Problem<T>;
    };

    template <typename T>
    class LinearConstraint {
    public:
        using Problem = qpcpp::Problem<T>;
        using Variable = qpcpp::Variable<T>;

        void setCoefficient(const Variable* variable, T coefficient);
        T getCoefficient(const Variable* variable) const;

        T min() const;
        T max() const;

    private:
        // linear constraints can only be created by Problem instances
        LinearConstraint(const Problem& context_problem, T min, T max);

        const Problem& context_problem_;
        T min_;
        T max_;

        std::unordered_map<const Variable*, T> coefficients_;

        friend class qpcpp::Problem<T>;
    };

    template <typename T>
    class CostFunction {
    public:
        using Problem = qpcpp::Problem<T>;
        using Variable = qpcpp::Variable<T>;

        // adds coefficient * first_variable * second_variable term to the problem.
        // status is not OK if first_variable or second_variable does not exist in
        // context problem
        void addQuadraticTerm(const Variable* first_variable,
                                      const Variable* second_variable,
                                      T coefficient);

        // get the coefficient of the term first_variable * second_variable. Return
        // status is not OK if first_variable or seocnd_variable does not exist on
        // context problem. Note that first_variable * second_variable is NOT a
        // different term than second_variable * first_variable and they would
        // return the same coefficient from this function.
        T getQuadraticCoefficient(
                const Variable* first_variable, const Variable* second_variable) const;

        // adds coefficient * variable term to the problem. status is not OK if
        // variable does not exist in the context problem
        void addLinearTerm(const Variable* variable, T coefficient);

        // get the cofficient of the term "variable". Return status is not OK if
        // variable does not exist in the context problem
        T getLinearCoefficient(const Variable* variable) const;

        // set the constant of the the cost function
        void add_constant(T constant);

        // get the constant of the cost function
        T constant() const;

        // sets the cost function to zero
        void setZero();

    private:
        // cost functions can only be created by problem instances
        CostFunction(const Problem& context_problem);

        const Problem& context_problem_;

        // [var-a-ptr][var-b-ptr] => coefficient. var-a-ptr <= var-b-ptr is enforced
        std::unordered_map<const Variable*, std::unordered_map<const Variable*, T>>
                quadratic_coefficients_;

        // [var-ptr] => coefficient.
        std::unordered_map<const Variable*, T> linear_coefficients_;

        // constant term of the cost
        T constant_;

        friend class qpcpp::Problem<T>;
    };

    template <typename T>
    class CPLEXSolver;

// Minimize cost_function_ subject to linear_constraints_ over
// variables_
    template <typename T>
    class Problem {
    public:
        using Variable = qpcpp::Variable<T>;
        using LinearConstraint = qpcpp::LinearConstraint<T>;
        using CostFunction = qpcpp::CostFunction<T>;

        Problem();

        // no copy construction and no assignment because pointers to fields are
        // important
        Problem(const Problem<T>& rhs) = delete;
        Problem<T>& operator=(const Problem<T>& rhs) = delete;

        std::size_t numVariables() const;
        Variable* addVariable(T min = std::numeric_limits<T>::lowest(),
                              T max = std::numeric_limits<T>::max());
        bool hasVariable(const Variable* variable_ptr) const;

        std::size_t numLinearConstraints() const;
        LinearConstraint* addLinearConstraint(
                T min = std::numeric_limits<T>::lowest(),
                T max = std::numeric_limits<T>::max());

        const std::forward_list<Variable>& variables() const;
        const std::forward_list<LinearConstraint>& linear_constraints() const;
        CostFunction* cost_function();

        void clearLinearConstraints();
        void setCostFunctionToZero();

    private:
        // We need pointer stability since pointers
        // to Variables are used as identifiers. std::forward_list
        // and std::list provides that while std::vector and
        // std::deque does not. we use std::forward_list because
        // singly linked list is enough for our purposes
        std::forward_list<Variable> variables_;

        // Since std::forward_list does not provide a size() function, we keep track
        // of number of variables using this member
        std::size_t num_variables_;

        // used for fast hasVariable check. contains all
        // variable iterators from the variables_ list
        std::unordered_set<const Variable*> variable_pointers_;

        // We need pointer stability since pointers to
        // LinearConstraints are used outside. std::forward_list
        // and std::list provides that while std::vector and
        // std::deque does not. we use std::forward_list because
        // singly linked list is enough for our purposes
        std::forward_list<LinearConstraint> linear_constraints_;

        // Since std::forward_list does not provide a size() function, we keep track
        // of number of constraints using this member
        std::size_t num_linear_constraints_;

        // quadratic cost function
        CostFunction cost_function_;

        // getter for mutable variables (only Problem class itself or friends can
        // use)
        std::forward_list<Variable>& mutable_variables();

        friend class CPLEXSolver<T>;
    };

}  // namespace qpcpp

#endif