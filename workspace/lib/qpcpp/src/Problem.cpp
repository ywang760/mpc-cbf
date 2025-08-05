#include <qpcpp/Problem.h>

namespace qpcpp {

template <typename T>
Variable<T>::Variable(const Problem& context_problem, T min, T max)
    : context_problem_(context_problem), min_(min), max_(max) {}

template <typename T>
const Problem<T>& Variable<T>::context_problem() const {
    return context_problem_;
}

template <typename T>
T Variable<T>::min() const {
    return min_;
}

template <typename T>
void Variable<T>::set_max(T new_max) {
    max_ = new_max;
}

template <typename T>
void Variable<T>::set_min(T new_min) {
    min_ = new_min;
}

template <typename T>
T Variable<T>::max() const {
    return max_;
}

template <typename T>
void Variable<T>::set_solution_value(T solution_value) {
    solution_value_ = solution_value;
}

template <typename T>
T Variable<T>::solution_value() const {
    return solution_value_;
}

template <typename T>
LinearConstraint<T>::LinearConstraint(const Problem& context_problem, T min, T max)
    : context_problem_(context_problem), min_(min), max_(max) {}

template <typename T>
void LinearConstraint<T>::setCoefficient(const Variable* variable, T coefficient) {
    if (!context_problem_.hasVariable(variable)) {
        throw std::runtime_error(
            "LinearConstraint::setCoefficient: variable does not exist in the convex problem");
    }

    coefficients_[variable] = coefficient;
}

template <typename T>
T LinearConstraint<T>::getCoefficient(const Variable* variable) const {
    if (!context_problem_.hasVariable(variable)) {
        throw std::runtime_error(
            "LinearConstraint::getCoefficient: variable does not exist in the convex problem");
    }

    if (typename std::unordered_map<const Variable*, T>::const_iterator it =
            coefficients_.find(variable);
        it == coefficients_.end()) {
        return T(0);
    } else {
        return it->second;
    }
}

template <typename T>
T LinearConstraint<T>::min() const {
    return min_;
}

template <typename T>
T LinearConstraint<T>::max() const {
    return max_;
}

template <typename T>
CostFunction<T>::CostFunction(const Problem& context_problem)
    : context_problem_(context_problem), constant_(0) {}

template <typename T>
void CostFunction<T>::addQuadraticTerm(const Variable* first_variable,
                                       const Variable* second_variable, T coefficient) {
    if (!context_problem_.hasVariable(first_variable)) {
        throw std::runtime_error(
            "CostFunction::addQuadraticTerm: first_variable does not exist in the context problem");
    }

    if (!context_problem_.hasVariable(second_variable)) {
        throw std::runtime_error("CostFunction::addQuadraticTerm: second_variable does not exist "
                                 "in the context problem");
    }

    if (first_variable > second_variable) {
        std::swap(first_variable, second_variable);
    }

    if (quadratic_coefficients_.find(first_variable) == quadratic_coefficients_.end()) {
        quadratic_coefficients_.emplace(first_variable, std::unordered_map<const Variable*, T>());
    }

    std::unordered_map<const Variable*, T>& first_variable_coefficients =
        quadratic_coefficients_.at(first_variable);

    if (first_variable_coefficients.find(second_variable) == first_variable_coefficients.end()) {
        first_variable_coefficients.emplace(second_variable, T(0));
    }

    quadratic_coefficients_[first_variable][second_variable] += coefficient;
}

template <typename T>
T CostFunction<T>::getQuadraticCoefficient(const Variable* first_variable,
                                           const Variable* second_variable) const {
    if (!context_problem_.hasVariable(first_variable)) {
        throw std::runtime_error("CostFunction::getQuadraticCoefficient: first_variable does not "
                                 "exist in the context problem");
    }

    if (!context_problem_.hasVariable(second_variable)) {
        throw std::runtime_error("CostFunction::getQuadraticCoefficient: second_variable does not "
                                 "exist in the context problem");
    }

    if (first_variable > second_variable) {
        std::swap(first_variable, second_variable);
    }

    if (typename std::unordered_map<const Variable*,
                                    std::unordered_map<const Variable*, T>>::const_iterator
            first_variable_coefficients_it = quadratic_coefficients_.find(first_variable);
        first_variable_coefficients_it != quadratic_coefficients_.end()) {
        const std::unordered_map<const Variable*, T>& first_variable_coefficients =
            first_variable_coefficients_it->second;

        if (typename std::unordered_map<const Variable*, T>::const_iterator coefficient_it =
                first_variable_coefficients.find(second_variable);
            coefficient_it != first_variable_coefficients.end()) {
            return coefficient_it->second;
        }
    }

    return T(0);
}

template <typename T>
void CostFunction<T>::addLinearTerm(const Variable* variable, T coefficient) {
    if (!context_problem_.hasVariable(variable)) {
        throw std::runtime_error(
            "CostFunction::addLinearTerm: variable does not exist in the context problem");
    }

    if (linear_coefficients_.find(variable) == linear_coefficients_.end()) {
        linear_coefficients_.emplace(variable, T(0));
    }

    linear_coefficients_[variable] += coefficient;
}

template <typename T>
T CostFunction<T>::getLinearCoefficient(const Variable* variable) const {
    if (!context_problem_.hasVariable(variable)) {
        throw std::runtime_error(
            "CostFunction::getLinearCoefficient: variable does not exist in the context problem");
    }

    if (linear_coefficients_.find(variable) == linear_coefficients_.end()) {
        return T(0);
    }

    return linear_coefficients_.at(variable);
}

template <typename T>
void CostFunction<T>::add_constant(T constant) {
    constant_ += constant;
}

template <typename T>
T CostFunction<T>::constant() const {
    return constant_;
}

template <typename T>
void CostFunction<T>::setZero() {
    quadratic_coefficients_.clear();
    linear_coefficients_.clear();
    constant_ = 0;
}

template <typename T>
Problem<T>::Problem() : num_variables_(0), num_linear_constraints_(0), cost_function_(*this) {}

template <typename T>
std::size_t Problem<T>::numVariables() const {
    return num_variables_;
}

template <typename T>
Variable<T>* Problem<T>::addVariable(T min, T max) {
    variables_.push_front(Variable(*this, min, max));
    Variable& new_variable = variables_.front();
    variable_pointers_.insert(&new_variable);
    ++num_variables_;

    return &new_variable;
}

template <typename T>
bool Problem<T>::hasVariable(const Variable* variable_ptr) const {
    return variable_pointers_.find(variable_ptr) != variable_pointers_.end();
}

template <typename T>
std::size_t Problem<T>::numLinearConstraints() const {
    return num_linear_constraints_;
}

template <typename T>
LinearConstraint<T>* Problem<T>::addLinearConstraint(T min, T max) {
    linear_constraints_.push_front(LinearConstraint(*this, min, max));
    ++num_linear_constraints_;

    LinearConstraint& new_linear_constraint = linear_constraints_.front();
    return &new_linear_constraint;
}

template <typename T>
const std::forward_list<Variable<T>>& Problem<T>::variables() const {
    return variables_;
}

template <typename T>
const std::forward_list<LinearConstraint<T>>& Problem<T>::linear_constraints() const {
    return linear_constraints_;
}

template <typename T>
CostFunction<T>* Problem<T>::cost_function() {
    return &cost_function_;
}

template <typename T>
std::forward_list<Variable<T>>& Problem<T>::mutable_variables() {
    return variables_;
}

template <typename T>
void Problem<T>::clearLinearConstraints() {
    linear_constraints_.clear();
    num_linear_constraints_ = 0;
}

template <typename T>
void Problem<T>::setCostFunctionToZero() {
    cost_function_.setZero();
}

template <typename T>
void Problem<T>::resetProblem() {
    setCostFunctionToZero();
    clearLinearConstraints();
}

template class Variable<double>;
template class Variable<float>;
template class LinearConstraint<double>;
template class LinearConstraint<float>;
template class CostFunction<double>;
template class CostFunction<float>;
template class Problem<double>;
template class Problem<float>;

} // namespace qpcpp