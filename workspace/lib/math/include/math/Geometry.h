#pragma once

#include <math/Types.h>
#include <model/DoubleIntegrator.h>
#include <Eigen/Dense>
#include <cmath>
#include <cassert>

namespace math {

/**
 * Find the closest point on an ellipse to the given robot position
 * @param robot_pos The position of the robot
 * @param target_mean The mean position of the target/ellipse center
 * @param target_cov The covariance matrix defining the ellipse
 * @return The closest point on the ellipse boundary
 */
Vector<double> closestPointOnEllipse(const VectorDIM<double, 3U>& robot_pos,
                                     const Vector<double>& target_mean,
                                     const Matrix<double>& target_cov);

/**
 * Check if a target is within the field of view of a robot
 * @param robot The robot position and orientation [x, y, yaw]
 * @param target The target position [x, y, ...]
 * @param fov The field of view angle (radians)
 * @param range The maximum detection range
 * @return True if the target is within the FOV, false otherwise
 */
bool insideFOV(const Eigen::VectorXd& robot, 
               const Eigen::VectorXd& target, 
               double fov, 
               double range);

/**
* Convert yaw angle to be within [-pi, pi]
* @param yaw The yaw angle to convert
* @return The yaw angle within range [-pi, pi]
*/
double convertYawInRange(double yaw);

/**
 * Rotate a control input from world frame to body frame
 * @param control_wf Control input in world frame
 * @param radian Rotation angle in radians
 * @return Control input in body frame
 */
template<unsigned int DIM>
VectorDIM<double, DIM> rotateControlInputToBodyFrame(const VectorDIM<double, DIM>& control_wf, double radian) {
    Matrix<double> R(3, 3);
    R << cos(radian), sin(radian), 0,
         -sin(radian), cos(radian), 0,
         0, 0, 1;
    return R * control_wf;
}

/**
 * Rotate a control input from body frame to world frame
 * @param control_bf Control input in body frame
 * @param radian Rotation angle in radians
 * @return Control input in world frame
 */
template<unsigned int DIM>
VectorDIM<double, DIM> rotateControlInputToWorldFrame(const VectorDIM<double, DIM>& control_bf, double radian) {
    Matrix<double> R(3, 3);
    R << cos(radian), sin(radian), 0,
         -sin(radian), cos(radian), 0,
         0, 0, 1;
    return R.transpose() * control_bf;
}


/**
 * Convert goal orientation to the closest equivalent orientation
 * to the current state to avoid large rotations
 * @param state The current state
 * @param goal The goal position and orientation
 * @return The goal with the closest equivalent orientation
 */
template<unsigned int DIM = 3U>
VectorDIM<double, DIM> convertToClosestYaw(
    const model::State<double, DIM>& state, 
    const VectorDIM<double, DIM>& goal) {
    
    using Vector = math::Vector<double>;
    
    double current_yaw = state.pos_(2);
    // Generate the candidate desire yaws
    Vector candidate_yaws(3);
    candidate_yaws << goal(2), goal(2) + 2 * M_PI, goal(2) - 2 * M_PI;
    
    Vector candidate_yaws_offset(3);
    candidate_yaws_offset << std::abs(candidate_yaws(0) - current_yaw), 
                             std::abs(candidate_yaws(1) - current_yaw), 
                             std::abs(candidate_yaws(2) - current_yaw);
    
    // Find the index of the minimum element
    int argmin_index = 0;
    double min_value = candidate_yaws_offset(0);
    
    for (int i = 1; i < candidate_yaws_offset.size(); ++i) {
        if (candidate_yaws_offset(i) < min_value) {
            min_value = candidate_yaws_offset(i);
            argmin_index = i;
        }
    }
    
    VectorDIM<double, DIM> converted_goal;
    converted_goal << goal(0), goal(1), candidate_yaws(argmin_index);
    return converted_goal;
}

} // namespace math
