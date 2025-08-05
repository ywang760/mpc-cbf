//
// Created by lishuo on 8/27/24.
//

#ifndef SEPARATINGHYPERPLANES_TYPES_H
#define SEPARATINGHYPERPLANES_TYPES_H

#include <math/Types.h>

namespace separating_hyperplanes {
template <typename T, unsigned int DIM>
using VectorDIM = math::VectorDIM<T, DIM>;

template <typename T, unsigned int DIM>
using Hyperplane = math::Hyperplane<T, DIM>;
} // namespace separating_hyperplanes

#endif //SEPARATINGHYPERPLANES_TYPES_H
