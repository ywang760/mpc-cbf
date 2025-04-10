#include <math/Helpers.h>

namespace math {

    template <typename T, unsigned int DIM>
    std::vector<VectorDIM<T, DIM>> cornerPoints(const AlignedBox<T, DIM>& box) {
        // Since there are 1U << DIM number of corner points, if there are more
        // than or equal to 32 dimensions 1U << 32 overflows.
        static_assert(DIM < 32U);
        std::vector<VectorDIM<T, DIM>> pts(1U << DIM);
        for (unsigned int i = 0; i < (1 << DIM); i++) {
            for (unsigned int d = 0; d < DIM; d++) {
                pts[i](d) = (i & (1 << d)) ? box.min()(d) : box.max()(d);
            }
        }
        return pts;
    }

    template <typename T, unsigned int DIM>
    Hyperplane<T, DIM> shiftHyperplane(const Hyperplane<T, DIM>& hyperplane,
                                       const AlignedBox<T, DIM>& box_at_zero) {
        using Hyperplane = Hyperplane<T, DIM>;
        using VectorDIM = VectorDIM<T, DIM>;

        Hyperplane shifted_hyperplane{hyperplane.normal(),
                                      std::numeric_limits<T>::lowest()};

        const std::vector<VectorDIM> corner_points =
                cornerPoints<T, DIM>(box_at_zero);

        for (const VectorDIM& corner_point : corner_points) {
            shifted_hyperplane.offset() = std::max(
                    shifted_hyperplane.offset(),
                    hyperplane.normal().dot(corner_point) + hyperplane.offset());
        }

        return shifted_hyperplane;
    }

    template <typename T, unsigned int DIM>
    AlignedBox<T, DIM> bufferAlignedBox(const AlignedBox<T, DIM>& aligned_box,
                                        const AlignedBox<T, DIM>& box_at_zero) {
        return AlignedBox<T, DIM>{aligned_box.min() - box_at_zero.min(),
                                  aligned_box.max() - box_at_zero.max()};
    }

    template <typename T, unsigned int DIM>
    Vector<T> stackVectorDIMs(const std::vector<VectorDIM<T, DIM>>& vectordims) {
        Vector<T> result(vectordims.size() * DIM);
        for (std::size_t vectordim_idx = 0; vectordim_idx < vectordims.size();
             ++vectordim_idx) {
            for (unsigned int dimension_idx = 0; dimension_idx < DIM;
                 ++dimension_idx) {
                result(vectordim_idx * DIM + dimension_idx) =
                        vectordims[vectordim_idx](dimension_idx);
            }
        }

        return result;
    }

    template <typename T, unsigned int DIM>
    std::vector<VectorDIM<T, DIM>> unstackVector(const Vector<T>& vector) {
        if (vector.rows() % DIM != 0) {
            throw std::runtime_error("[unstackVector]: vector size is not divisible by dimension.");
        }

        std::size_t vectordim_size = vector.rows() / DIM;
        std::vector<VectorDIM<T, DIM>> vectordims(vectordim_size);
        for (std::size_t vectordim_idx = 0; vectordim_idx < vectordim_size;
             ++vectordim_idx) {
            for (unsigned int dimension_idx = 0; dimension_idx < DIM;
                 ++dimension_idx) {
                vectordims[vectordim_idx](dimension_idx) =
                        vector(vectordim_idx * DIM + dimension_idx);
            }
        }

        return vectordims;
    }

    template <typename T>
    bool isApproximatelyEqual(T lhs, T rhs, T tolerance) {
        static_assert(
                std::is_same<T, double>::value || std::is_same<T, float>::value,
                "[isAlmostEqual]: T should be in {float, double.}");

        return (std::fabs(lhs - rhs) <= tolerance);
    }

    template <typename T, unsigned int DIM>
    bool isApproximatelyEqual(const VectorDIM<T, DIM>& lhs,
                              const VectorDIM<T, DIM>& rhs, T tolerance) {
        for (unsigned int d = 0; d < DIM; ++d) {
            if (!isApproximatelyEqual<T>(lhs(d), rhs(d), tolerance)) {
                return false;
            }
        }
        return true;
    }

    template <typename T, unsigned int R, unsigned int C>
    bool isApproximatelyEqual(const MatrixRC<T, R, C>& lhs,
                              const MatrixRC<T, R, C>& rhs, T tolerance) {
        for (unsigned int r = 0; r < R; ++r) {
            for (unsigned int c = 0; c < C; ++c) {
                if (!isApproximatelyEqual<T>(lhs(r, c), rhs(r, c), tolerance)) {
                    return false;
                }
            }
        }

        return true;
    }

    template <typename Derived1, typename Derived2, typename T>
    bool isApproximatelyEqual(const Eigen::MatrixBase<Derived1> &lhs,
                              const Eigen::MatrixBase<Derived2> &rhs,
                              T tolerance)
    {
            // Check dimensions match
            if (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols())
            {
                    return false;
            }

            // Check each element is approximately equal
            return ((lhs - rhs).array().abs() <= tolerance).all();
    }

    template <typename T, unsigned int DIM>
    bool isInPositiveSide(const Hyperplane<T, DIM>& hyperplane,
                          const AlignedBox<T, DIM>& box) {
        using VectorDIM = math::VectorDIM<T, DIM>;
        const std::vector<VectorDIM> corner_points = cornerPoints<T, DIM>(box);
        for (const VectorDIM& corner_point : corner_points) {
            if (hyperplane.signedDistance(corner_point) <= T(0.0)) {
                return false;
            }
        }
        return true;
    }

    template <typename T, unsigned int DIM>
    bool isInPositiveSide(const Hyperplane<T, DIM>& hyperplane,
                          const std::vector<VectorDIM<T, DIM>>& points) {
        for (const VectorDIM<T, DIM>& point : points) {
            if (hyperplane.signedDistance(point) <= T(0.0)) {
                return false;
            }
        }
        return true;
    }

    template <typename T, unsigned int DIM>
    Hyperplane<T, DIM> snapHyperplane(
            const Hyperplane<T, DIM>& hyperplane,
            const std::vector<VectorDIM<T, DIM>>& points) {
        T closest_distance = std::numeric_limits<T>::max();
        for (const VectorDIM<T, DIM>& point : points) {
            closest_distance =
                    std::min(closest_distance, hyperplane.signedDistance(point));
        }
        return Hyperplane<T, DIM>{hyperplane.normal(),
                                  hyperplane.offset() - closest_distance};
    }


    template <typename T, unsigned int DIM>
    std::vector<Hyperplane<T, DIM>> boundingHyperplanes(
            const AlignedBox<T, DIM>& box) {
        using Hyperplane = Hyperplane<T, DIM>;
        using VectorDIM = VectorDIM<T, DIM>;
        std::vector<Hyperplane> bounding_hyperplanes;

        for (unsigned int d = 0; d < DIM; ++d) {
            VectorDIM direction = VectorDIM::Constant(0.0);
            direction(d) = 1.0;
            bounding_hyperplanes.emplace_back(direction, -box.max()(d));
            bounding_hyperplanes.emplace_back(-direction, box.min()(d));
        }
        return bounding_hyperplanes;
    }

    template <typename T, unsigned int DIM>
    VectorDIM<T, DIM> linearInterpolate(
            const std::pair<T, VectorDIM<T, DIM>>& point1,
            const std::pair<T, VectorDIM<T, DIM>>& point2, T t) {
        T ratio = (t - point1.first) / (point2.first - point1.first);
        return point1.second + ratio * (point2.second - point1.second);
    }

    template <typename T, unsigned int DIM>
    VectorDIM<T, DIM> computeAPerpendicularVector(const VectorDIM<T, DIM>& vector) {
        using VectorDIM = math::VectorDIM<T, DIM>;
        if (math::isApproximatelyEqual<T, DIM>(vector, VectorDIM::Zero())) {
            return VectorDIM::Zero();
        }

        unsigned int non_zero_index = 0;
        while (math::isApproximatelyEqual(vector(non_zero_index), T(0.0))) {
            ++non_zero_index;
        }

        unsigned int another_index =
                (non_zero_index == 0 ? 1 : (non_zero_index - 1));

        VectorDIM perpendicular_vector = VectorDIM::Zero();
        perpendicular_vector(non_zero_index) = vector(another_index);
        perpendicular_vector(another_index) = -vector(non_zero_index);

        return perpendicular_vector;
    }

    template std::vector<VectorDIM<double, 2U>> cornerPoints<double, 2U>(
            const AlignedBox<double, 2U>& box);
    template std::vector<VectorDIM<double, 3U>> cornerPoints<double, 3U>(
            const AlignedBox<double, 3U>& box);
    template std::vector<VectorDIM<float, 2U>> cornerPoints<float, 2U>(
            const AlignedBox<float, 2U>& box);
    template std::vector<VectorDIM<float, 3U>> cornerPoints<float, 3U>(
            const AlignedBox<float, 3U>& box);

    template Hyperplane<double, 3U> shiftHyperplane<double, 3U>(
            const Hyperplane<double, 3U>& hyperplane,
            const AlignedBox<double, 3U>& bounding_box_at_zero);
    template Hyperplane<float, 3U> shiftHyperplane<float, 3U>(
            const Hyperplane<float, 3U>& hyperplane,
            const AlignedBox<float, 3U>& bounding_box_at_zero);
    template Hyperplane<double, 2U> shiftHyperplane<double, 2U>(
            const Hyperplane<double, 2U>& hyperplane,
            const AlignedBox<double, 2U>& bounding_box_at_zero);
    template Hyperplane<float, 2U> shiftHyperplane<float, 2U>(
            const Hyperplane<float, 2U>& hyperplane,
            const AlignedBox<float, 2U>& bounding_box_at_zero);

    template AlignedBox<double, 3U> bufferAlignedBox<double, 3U>(
            const AlignedBox<double, 3U>& aligned_box,
            const AlignedBox<double, 3U>& box_at_zero);
    template AlignedBox<float, 3U> bufferAlignedBox<float, 3U>(
            const AlignedBox<float, 3U>& aligned_box,
            const AlignedBox<float, 3U>& box_at_zero);
    template AlignedBox<double, 2U> bufferAlignedBox<double, 2U>(
            const AlignedBox<double, 2U>& aligned_box,
            const AlignedBox<double, 2U>& box_at_zero);
    template AlignedBox<float, 2U> bufferAlignedBox<float, 2U>(
            const AlignedBox<float, 2U>& aligned_box,
            const AlignedBox<float, 2U>& box_at_zero);

    template Vector<double> stackVectorDIMs<double, 2U>(
            const std::vector<VectorDIM<double, 2U>>& vectordims);
    template Vector<double> stackVectorDIMs<double, 3U>(
            const std::vector<VectorDIM<double, 3U>>& vectordims);
    template Vector<float> stackVectorDIMs<float, 2U>(
            const std::vector<VectorDIM<float, 2U>>& vectordims);
    template Vector<float> stackVectorDIMs<float, 3U>(
            const std::vector<VectorDIM<float, 3U>>& vectordims);

    template std::vector<VectorDIM<double, 3U>>
    unstackVector<double, 3U>(const Vector<double>& vector);
    template std::vector<VectorDIM<double, 2U>>
    unstackVector<double, 2U>(const Vector<double>& vector);
    template std::vector<VectorDIM<float, 3U>>
    unstackVector<float, 3U>(const Vector<float>& vector);
    template std::vector<VectorDIM<float, 2U>>
    unstackVector<float, 2U>(const Vector<float>& vector);

    template bool isApproximatelyEqual<double>(double lhs, double rhs,
                                               double tolerance);
    template bool isApproximatelyEqual<float>(float lhs, float rhs,
                                              float tolerance);

    template bool isApproximatelyEqual<double, 2U>(const VectorDIM<double, 2U>& lhs,
                                                   const VectorDIM<double, 2U>& rhs,
                                                   double tolerance);
    template bool isApproximatelyEqual<float, 2U>(const VectorDIM<float, 2U>& lhs,
                                                  const VectorDIM<float, 2U>& rhs,
                                                  float tolerance);
    template bool isApproximatelyEqual<double, 3U>(const VectorDIM<double, 3U>& lhs,
                                                   const VectorDIM<double, 3U>& rhs,
                                                   double tolerance);
    template bool isApproximatelyEqual<float, 3U>(const VectorDIM<float, 3U>& lhs,
                                                  const VectorDIM<float, 3U>& rhs,
                                                  float tolerance);

    template bool isApproximatelyEqual<double, 2U, 2U>(
            const MatrixRC<double, 2U, 2U>& lhs, const MatrixRC<double, 2U, 2U>& rhs,
            double tolerance);
    template bool isApproximatelyEqual<float, 2U, 2U>(
            const MatrixRC<float, 2U, 2U>& lhs, const MatrixRC<float, 2U, 2U>& rhs,
            float tolerance);
    template bool isApproximatelyEqual<double, 3U, 3U>(
            const MatrixRC<double, 3U, 3U>& lhs, const MatrixRC<double, 3U, 3U>& rhs,
            double tolerance);
    template bool isApproximatelyEqual<float, 3U, 3U>(
            const MatrixRC<float, 3U, 3U>& lhs, const MatrixRC<float, 3U, 3U>& rhs,
            float tolerance);
    template bool isApproximatelyEqual<double, 4U, 4U>(
            const MatrixRC<double, 4U, 4U>& lhs, const MatrixRC<double, 4U, 4U>& rhs,
            double tolerance);
    template bool isApproximatelyEqual<float, 4U, 4U>(
            const MatrixRC<float, 4U, 4U>& lhs, const MatrixRC<float, 4U, 4U>& rhs,
            float tolerance);
    template bool isApproximatelyEqual<double, 6U, 6U>(
            const MatrixRC<double, 6U, 6U>& lhs, const MatrixRC<double, 6U, 6U>& rhs,
            double tolerance);
    template bool isApproximatelyEqual<float, 6U, 6U>(
            const MatrixRC<float, 6U, 6U>& lhs, const MatrixRC<float, 6U, 6U>& rhs,
            float tolerance);

    template bool isInPositiveSide<double, 3U>(
            const Hyperplane<double, 3U>& hyperplane,
            const AlignedBox<double, 3U>& box);
    template bool isInPositiveSide<double, 2U>(
            const Hyperplane<double, 2U>& hyperplane,
            const AlignedBox<double, 2U>& box);
    template bool isInPositiveSide<float, 3U>(
            const Hyperplane<float, 3U>& hyperplane, const AlignedBox<float, 3U>& box);
    template bool isInPositiveSide<float, 2U>(
            const Hyperplane<float, 2U>& hyperplane, const AlignedBox<float, 2U>& box);

    template bool isInPositiveSide<double, 3U>(
            const Hyperplane<double, 3U>& hyperplane,
            const std::vector<VectorDIM<double, 3U>>& points);
    template bool isInPositiveSide<double, 2U>(
            const Hyperplane<double, 2U>& hyperplane,
            const std::vector<VectorDIM<double, 2U>>& points);
    template bool isInPositiveSide<float, 3U>(
            const Hyperplane<float, 3U>& hyperplane,
            const std::vector<VectorDIM<float, 3U>>& points);
    template bool isInPositiveSide<float, 2U>(
            const Hyperplane<float, 2U>& hyperplane,
            const std::vector<VectorDIM<float, 2U>>& points);

    template Hyperplane<double, 3U> snapHyperplane<double, 3U>(
            const Hyperplane<double, 3U>& hyperplane,
            const std::vector<VectorDIM<double, 3U>>& points);
    template Hyperplane<float, 3U> snapHyperplane<float, 3U>(
            const Hyperplane<float, 3U>& hyperplane,
            const std::vector<VectorDIM<float, 3U>>& points);
    template Hyperplane<double, 2U> snapHyperplane<double, 2U>(
            const Hyperplane<double, 2U>& hyperplane,
            const std::vector<VectorDIM<double, 2U>>& points);
    template Hyperplane<float, 2U> snapHyperplane<float, 2U>(
            const Hyperplane<float, 2U>& hyperplane,
            const std::vector<VectorDIM<float, 2U>>& points);

    template std::vector<Hyperplane<double, 2U>> boundingHyperplanes<double, 2U>(
            const AlignedBox<double, 2U>& box);
    template std::vector<Hyperplane<double, 3U>> boundingHyperplanes<double, 3U>(
            const AlignedBox<double, 3U>& box);
    template std::vector<Hyperplane<float, 2U>> boundingHyperplanes<float, 2U>(
            const AlignedBox<float, 2U>& box);
    template std::vector<Hyperplane<float, 3U>> boundingHyperplanes<float, 3U>(
            const AlignedBox<float, 3U>& box);

    template VectorDIM<double, 3U> linearInterpolate<double, 3U>(
            const std::pair<double, VectorDIM<double, 3U>>&,
            const std::pair<double, VectorDIM<double, 3U>>&, double);
    template VectorDIM<float, 3U> linearInterpolate<float, 3U>(
            const std::pair<float, VectorDIM<float, 3U>>&,
            const std::pair<float, VectorDIM<float, 3U>>&, float);
    template VectorDIM<double, 2U> linearInterpolate<double, 2U>(
            const std::pair<double, VectorDIM<double, 2U>>&,
            const std::pair<double, VectorDIM<double, 2U>>&, double);
    template VectorDIM<float, 2U> linearInterpolate<float, 2U>(
            const std::pair<float, VectorDIM<float, 2U>>&,
            const std::pair<float, VectorDIM<float, 2U>>&, float);

    template VectorDIM<double, 3U> computeAPerpendicularVector<double, 3U>(
            const VectorDIM<double, 3U>& v);
    template VectorDIM<float, 3U> computeAPerpendicularVector<float, 3U>(
            const VectorDIM<float, 3U>& v);
    template VectorDIM<double, 2U> computeAPerpendicularVector<double, 2U>(
            const VectorDIM<double, 2U>& v);
    template VectorDIM<float, 2U> computeAPerpendicularVector<float, 2U>(
            const VectorDIM<float, 2U>& v);

    template bool isApproximatelyEqual<Eigen::Matrix<double, 2, 1>, Eigen::Matrix<double, 2, 1>, double>(
        const Eigen::MatrixBase<Eigen::Matrix<double, 2, 1>> &lhs,
        const Eigen::MatrixBase<Eigen::Matrix<double, 2, 1>> &rhs,
        double tolerance);

    template bool isApproximatelyEqual<Eigen::Matrix<float, 2, 1>, Eigen::Matrix<float, 2, 1>, float>(
        const Eigen::MatrixBase<Eigen::Matrix<float, 2, 1>> &lhs,
        const Eigen::MatrixBase<Eigen::Matrix<float, 2, 1>> &rhs,
        float tolerance);

    // For Vector3d and Vector3f
    template bool isApproximatelyEqual<Eigen::Matrix<double, 3, 1>, Eigen::Matrix<double, 3, 1>, double>(
        const Eigen::MatrixBase<Eigen::Matrix<double, 3, 1>> &lhs,
        const Eigen::MatrixBase<Eigen::Matrix<double, 3, 1>> &rhs,
        double tolerance);

    template bool isApproximatelyEqual<Eigen::Matrix<float, 3, 1>, Eigen::Matrix<float, 3, 1>, float>(
        const Eigen::MatrixBase<Eigen::Matrix<float, 3, 1>> &lhs,
        const Eigen::MatrixBase<Eigen::Matrix<float, 3, 1>> &rhs,
        float tolerance);

    // For explicitly comparing with expression types
    template bool isApproximatelyEqual<Eigen::CwiseBinaryOp<Eigen::internal::scalar_sum_op<double, double>,
                                                            const Eigen::Matrix<double, 2, 1>, const Eigen::Matrix<double, 2, 1>>,
                                       Eigen::Matrix<double, 2, 1>, double>(
        const Eigen::MatrixBase<Eigen::CwiseBinaryOp<Eigen::internal::scalar_sum_op<double, double>,
                                                     const Eigen::Matrix<double, 2, 1>, const Eigen::Matrix<double, 2, 1>>> &lhs,
        const Eigen::MatrixBase<Eigen::Matrix<double, 2, 1>> &rhs,
        double tolerance);

    template bool isApproximatelyEqual<Eigen::CwiseBinaryOp<Eigen::internal::scalar_sum_op<float, float>,
                                                            const Eigen::Matrix<float, 2, 1>, const Eigen::Matrix<float, 2, 1>>,
                                       Eigen::Matrix<float, 2, 1>, float>(
        const Eigen::MatrixBase<Eigen::CwiseBinaryOp<Eigen::internal::scalar_sum_op<float, float>,
                                                     const Eigen::Matrix<float, 2, 1>, const Eigen::Matrix<float, 2, 1>>> &lhs,
        const Eigen::MatrixBase<Eigen::Matrix<float, 2, 1>> &rhs,
        float tolerance);

}  // namespace math