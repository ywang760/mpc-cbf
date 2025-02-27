//
// Created by lishuo on 8/20/24.
//
#include <iostream>
#include <memory>
#include <math/Types.h>
#include <cmath>

template <typename T, unsigned int DIM>
class c0 {
public:
    int x_ = 0;

};

template <typename T>
class c1 : public c0<T, 3U> {
public:
    using c0<T, 3U>::x_;
    c1(int x, int y) {
        x_ = x;
        y_ = y;
    };
    int y_;
};

template <typename T, unsigned int DIM>
class d0 {
public:
    d0(std::shared_ptr<c0<T, DIM>> c_ptr) {c_ptr_ = c_ptr;};

    std::shared_ptr<c0<T, DIM>> c_ptr_;
};

double dmod(double x, double y) {
    return x - (int)(x/y) * y;
}

int main() {
    c1<double> p = {1, 2};
    std::shared_ptr<c1<double>> var = std::make_shared<c1<double>>(p);
    d0<double, 3U> d_var(std::move(var));
    std::cout << "\n";

    int k = 10;
    double h = 0.1;
    using Vector = math::Vector<double>;
    Vector x_lin_space = Vector::LinSpaced(k, 0, (k-1)*h);
    std::cout << "line spaced: " << x_lin_space << "\n";

    int x = 10;
    int* x_ptr = &x;
    int* x_ptr_1 = x_ptr;
    std::vector<int*> x_ptrs;
    x_ptrs.push_back(x_ptr);
    x_ptrs.push_back(x_ptr_1);

    std::cout << "x: " << x_ptr << "\n";
    std::cout << "x: " << x_ptr_1 << "\n";
    std::cout << "x: " << *x_ptr << "\n";
    std::cout << "x: " << *x_ptr_1 << "\n";

    for (size_t i = 0; i < x_ptrs.size(); ++i) {
        std::cout << "ptrs: " << x_ptrs.at(i) << "\n";
        std::cout << "x " << *x_ptrs.at(i) << "\n";
    }

    double a = 5.20;
    double b = 0.40;
    double test_mod = a / b - int(a / b);
    double test_fmod = std::fmod(a,b);
    double test_dmod = dmod(a,b);
    std::cout << "mod result: " << test_mod << "\n";
    std::cout << "fmod result: " << test_fmod << "\n";
    std::cout << "dmod result: " << test_dmod << "\n";

    return 0;
}