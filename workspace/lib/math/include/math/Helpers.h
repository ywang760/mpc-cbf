#ifndef MATH_HELPERS_H
#define MATH_HELPERS_H

#include <math/Types.h>
#include <memory>
#include <vector>
#include <limits>

namespace math {
    /**
     * @brief Returns the corner points of the given box
     *
     * @tparam T Floating point number type
     * @tparam DIM Dimension that the box resides in
     * @param box Box for which corner points are computed
     * @return std::vector<VectorDIM<T, DIM>> Corner points of the box
     */
    template <typename T, unsigned int DIM>
    std::vector<VectorDIM<T, DIM>> cornerPoints(const AlignedBox<T, DIM>& box);

    /**
     * @brief Shifts hyperplane so that when the position of the given box is in the
     * negative side of the resulting hyperplane, box itself will be in the negative
     * side of the given hyperplane.
     *
     * @details Box is assumed to be translated using Minkowski sum between the
     * box_at_zero and {position}.
     *
     * @tparam T Floating point number type
     * @tparam DIM Dimension that the objects reside in
     * @param hyperplane Hyperplane to be shifted
     * @param box_at_zero Box when placed at position 0
     * @return Hyperplane<T, DIM> Shifted hyperplane
     */
    template <typename T, unsigned int DIM>
    Hyperplane<T, DIM> shiftHyperplane(const Hyperplane<T, DIM>& hyperplane,
                                       const AlignedBox<T, DIM>& box_at_zero);

    /**
     * @brief Buffers aligned_box so that when the position of the given box_at_zero
     * is inside of the resulting box, the box_at_zero shifted by its position will
     * be inside of the given aligned_box.
     *
     * @tparam T Floating point number type
     * @tparam DIM Dimension that the objects reside in
     * @param aligned_box AlignedBox to be buffered
     * @param box_at_zero Box when placed at position 0
     * @return AlignedBox<T, DIM> Buffered aligned box
     */
    template <typename T, unsigned int DIM>
    AlignedBox<T, DIM> bufferAlignedBox(const AlignedBox<T, DIM>& aligned_box,
                                        const AlignedBox<T, DIM>& box_at_zero);

    /**
     * @brief Stacks VectorDIMs into a single Vector
     *
     * @tparam T Floating point number type
     * @tparam DIM Dimension of VectorDIMs
     * @param vectordims VectorDIMs that are stacked
     * @return Vector<T> Stacked Vector of VectorDIMs
     */
    template <typename T, unsigned int DIM>
    Vector<T> stackVectorDIMs(const std::vector<VectorDIM<T, DIM>>& vectordims);

    /**
     * @brief Unstack a vector to VectorDIMs such that each DIM number of rows in
     * vector becomes a VectorDIM.
     *
     * @tparam T Floating point number type
     * @tparam DIM Dimension of the VectorDIMs
     * @param vector Vector to be unstacked
     * @return std::vector<VectorDIM<T, DIM>> std::vector of
     * unstacked VectorDIMs. Return status is not OK if vector's number of rows is
     * not divisible by DIM.
     */
    template <typename T, unsigned int DIM>
    std::vector<VectorDIM<T, DIM>> unstackVector(const Vector<T>& vector);

    /**
     * @brief Returns whether two floating point numbers are approximately equal
     *
     * @tparam T Floating point number type
     * @param lhs Left hand side of the comparison
     * @param rhs Right hand side of the comparison
     * @param tolerance Tolerance while checking equality. If difference of lhs and
     * rhs is less than or equal to tolerance, they are equal.
     * @return true If lhs is approximately equal to rhs
     * @return false If lhs is not approximately equal to rhs
     */
    template <typename T>
    bool isApproximatelyEqual(T lhs, T rhs, T tolerance = 1e-10);

    /**
     * @brief Returns whether two VectorDIMs are component-wise approximately equal
     *
     * @tparam T Floating point number type
     * @tparam DIM Dimension of the vectors
     * @param lhs Left hand side of the comparison
     * @param rhs Right hand side of the comparison
     * @param tolerance Tolerance while checking equality. If difference of any
     * coordinate of lhs and rhs is not less than or equal to the tolerance, lhs and
     * rhs are unequal.
     * @return true If lhs is component-wise approximately equal to rhs
     * @return false If lhs is not component-wise approximately equal to rhs
     */
    template <typename T, unsigned int DIM>
    bool isApproximatelyEqual(const VectorDIM<T, DIM>& lhs,
                              const VectorDIM<T, DIM>& rhs, T tolerance = 1e-10);

    /**
     * @brief Generic function to check if two Eigen expressions are approximately equal
     *
     * @tparam Derived1 First Eigen derived type
     * @tparam Derived2 Second Eigen derived type
     * @tparam T Scalar type used for tolerance
     * @param lhs Left-hand side Eigen expression
     * @param rhs Right-hand side Eigen expression
     * @param tolerance Maximum allowed difference between elements
     * @return true if expressions are approximately equal
     * @return false otherwise
     */
    template <typename Derived1, typename Derived2, typename T = typename Derived1::Scalar>
    bool isApproximatelyEqual(const Eigen::MatrixBase<Derived1> &lhs,
                              const Eigen::MatrixBase<Derived2> &rhs,
                              T tolerance = static_cast<T>(1e-10));

    /**
     * @brief Return whether two MatrixRCs are component-wise approximately equal
     *
     * @tparam T Floating point number type
     * @tparam R Number of rows of matrices
     * @tparam C Number of columns of matrices
     * @param lhs Left hand side of the comparison
     * @param rhs Right hand side of the comparison
     * @param tolerance Tolerance while checking equality. If difference of any
     * element of lhs and rhs is not less than or equal to the tolerance, lhs and
     * rhs are unequal.
     * @return true If lhs is component-wise approximately equal to rhs
     * @return false If lhs is not component-wise approximately equal to rhs
     */
    template <typename T, unsigned int R, unsigned int C>
    bool isApproximatelyEqual(const MatrixRC<T, R, C>& lhs,
                              const MatrixRC<T, R, C>& rhs, T tolerance = 1e-10);

    /**
     * @brief Return whether box is completely in the positive side of the
     * hyperplane
     *
     * @tparam T Floating point number type
     * @tparam DIM Dimension the objects are contained in
     * @param hyperplane Hyperplane against which box's placement is checked
     * @param box Box for which placement is checked
     * @return true If box is completely in the positive side of the hyperplane
     * @return false If box is not completely in the positive side of the hyperplane
     * (i.e. either completely in the negative side or intersects with the
     * hyperplane)
     */
    template <typename T, unsigned int DIM>
    bool isInPositiveSide(const Hyperplane<T, DIM>& hyperplane,
                          const AlignedBox<T, DIM>& box);

    /**
     * @brief Return whether all points are in the positive side of the hyperplane
     *
     * @tparam T Floating point number type
     * @tparam DIM Dimension the points are contained in
     * @param hyperplane Hyperplane against which points' placement is checked
     * @param points Points for which placement is checked
     * @return true If all points are in the positive side of the hyperplane
     * @return false If any point is not in the positive side of the hyperplane
     */
    template <typename T, unsigned int DIM>
    bool isInPositiveSide(const Hyperplane<T, DIM>& hyperplane,
                          const std::vector<VectorDIM<T, DIM>>& points);

    /**
     * @brief Snap hyperplane to the convex hull of the points so that all points
     * are on the non-negative side of the hyperplane but any additional positive
     * shift along the normal causes at least one point to be in the negative side
     * of the hyperplane.
     *
     * @details The resulting hyperplane has the same normal as the input
     * hyperplane.
     *
     * @tparam T Floating point number type
     * @tparam DIM Ambient dimension objects reside in
     * @param hyperplane Hyperplane to be snapped
     * @param points Points that hyperplanes are snapped to
     * @return Hyperplane<T, DIM> Snapped hyperplane
     */
    template <typename T, unsigned int DIM>
    Hyperplane<T, DIM> snapHyperplane(const Hyperplane<T, DIM>& hyperplane,
                                      const std::vector<VectorDIM<T, DIM>>& points);


    /**
     * @brief Return the vector of hyperplanes bounding the given aligned box
     *
     * @tparam T Floating point type
     * @tparam DIM Ambient dimension of the objects
     * @param box Aligned box for which bounding hyperplanes are computed
     * @return std::vector<Hyperplane<T, DIM>> Vector of hyperplanes bounding the
     * aligned box
     */
    template <typename T, unsigned int DIM>
    std::vector<Hyperplane<T, DIM>> boundingHyperplanes(
            const AlignedBox<T, DIM>& box);

    /**
     * @brief Linear interpolate on the line defined by two points
     *
     * @tparam T Floating point type
     * @tparam DIM Ambient dimension
     * @param point1 First point and its parameter
     * @param point2 Second point and its parameter
     * @param t Queried parameter
     * @return VectorDIM<T, DIM> Interpolated point at the queried parameter
     */
    template <typename T, unsigned int DIM>
    VectorDIM<T, DIM> linearInterpolate(
            const std::pair<T, VectorDIM<T, DIM>>& point1,
            const std::pair<T, VectorDIM<T, DIM>>& point2, T t);

    /**
     * @brief Compute a perpendicular vector to the given vector
     *
     * @details The returned vector is not necessarily normalized. If vector is
     * zero, return value is also zero
     *
     * @tparam T Floating point type
     * @tparam DIM Ambient dimension
     * @param vector Vector for which a perpendicular vector is computed
     * @return VectorDIM<T, DIM> Perpendicular vector
     */
    template <typename T, unsigned int DIM>
    VectorDIM<T, DIM> computeAPerpendicularVector(const VectorDIM<T, DIM>& vector);
}  // namespace math

#endif //MATH_HELPERS_H