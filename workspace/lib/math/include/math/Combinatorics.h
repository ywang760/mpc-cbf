//
// Created by lishuo on 2/3/24.
//

#ifndef MATH_COMBINATORICS_H
#define MATH_COMBINATORICS_H

#include <iostream>

namespace math {

uint64_t fac(uint64_t n);

uint64_t comb(uint64_t n, uint64_t k);

uint64_t perm(uint64_t n, uint64_t k);

template <typename T, typename U>
T pow(T base, U exp);

} // namespace math

#endif //MATH_COMBINATORICS_H
