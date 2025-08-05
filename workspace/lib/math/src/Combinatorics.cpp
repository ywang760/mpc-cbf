//
// Created by lishuo on 2/3/24.
//

#include <cmath>
#include <math/Combinatorics.h>

namespace math {

uint64_t fac(uint64_t n) {
    if (n > 20) {
        throw std::runtime_error("warning: fac overflow");
    }
    uint64_t result = 1;
    for (uint64_t i = 2; i <= n; i++) {
        result *= i;
    }
    return result;
}

uint64_t comb(uint64_t n, uint64_t k) {
    if (k > n) return 0;
    k = std::min(k, n - k);
    uint64_t top = 1;
    uint64_t bottom = 1;
    for (uint64_t i = 0; i < k; i++) {
        bottom *= (i + 1);
        top *= (n - i);
    }
    return top / bottom;
}

uint64_t perm(uint64_t n, uint64_t k) {
    if (k > n) return 0;
    uint64_t result = 1;
    for (uint64_t i = n - k + 1; i <= n; ++i) {
        result *= i;
    }
    return result;
}

template <typename T, typename U>
T pow(T base, U exp) {
    if (base == 0 && exp == 0) {
        return 1;
    }
    if (base == 0 && exp < 0) {
        throw std::runtime_error("pow exp is negative");
    }
    return std::pow(base, exp);
}

template double pow<double, uint64_t>(double, uint64_t);
template float pow<float, uint64_t>(float, uint64_t);
template double pow<double, int>(double, int);
template float pow<float, int>(float, int);
template double pow<double, unsigned int>(double, unsigned int);
template float pow<float, unsigned int>(float, unsigned int);

} // namespace math