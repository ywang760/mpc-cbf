#include <math/Geometry.h>
#include <Eigen/Eigenvalues>
#include <cmath>

namespace math {

Vector<double> closestPointOnEllipse(const VectorDIM<double, 3U>& robot_pos,
                                     const Vector<double>& target_mean,
                                     const Matrix<double>& target_cov) {
    if (!std::isinf(target_cov(0, 0))) {
        Eigen::EigenSolver<Matrix<double>> es(target_cov.block(0, 0, 3U-1, 3U-1));
        Vector<double> eigenvalues = es.eigenvalues().real();
        Matrix<double> eigenvectors = es.eigenvectors().real();

        // s = 4.605 for 90% confidence interval
        // s = 5.991 for 95% confidence interval
        // s = 9.210 for 99% confidence interval
        double s = 4.605;
        double a = sqrt(s * eigenvalues(0)); // major axis
        double b = sqrt(s * eigenvalues(1)); // minor axis

        // a could be smaller than b, so swap them
        if (a < b) {
            double temp = a;
            a = b;
            b = temp;
        }

        int m = 0; // higher eigenvalue index
        int l = 1; // lower eigenvalue index
        if (eigenvalues(1) > eigenvalues(0)) {
            m = 1;
            l = 0;
        }

        double theta = atan2(eigenvectors(1, m), eigenvectors(0, m)); // angle of the major axis wrt positive x-axis (ccw rotation)
        if (theta < 0.0) {
            theta += M_PI;
        } // angle in [0, 2pi]

        double slope = atan2(-target_mean(1) + robot_pos(1), -target_mean(0) + robot_pos(0));
        double x_n = target_mean(0) + a * cos(slope - theta) * cos(theta) - b * sin(slope - theta) * sin(theta);
        double y_n = target_mean(1) + a * cos(slope - theta) * sin(theta) + b * sin(slope - theta) * cos(theta);

        Vector<double> p_near(2);
        p_near << x_n, y_n;
        return p_near;
    } else {
        Vector<double> p_near(2);
        p_near << 0, 0;
        return p_near;
    }
}

bool insideFOV(const Eigen::VectorXd& robot, const Eigen::VectorXd& target, double fov, double range) {
    double yaw = robot(2);

    Eigen::Matrix3d R;
    R << cos(yaw), sin(yaw), 0.0,
         -sin(yaw), cos(yaw), 0.0,
         0.0, 0.0, 1.0;
    Eigen::VectorXd t_local = R.block<2,2>(0,0) * (target.head(2) - robot.head(2));
    double dist = t_local.norm();
    double angle = abs(atan2(t_local(1), t_local(0)));
    if (angle <= 0.5*fov && dist <= range) {
        return true;
    } else {
        return false;
    }
}

double convertYawInRange(double yaw) {
    assert(yaw < 2 * M_PI && yaw > -2 * M_PI);
    if (yaw > M_PI) {
        return yaw - 2 * M_PI;
    } else if (yaw < -M_PI) {
        return yaw + 2 * M_PI;
    } else {
        return yaw;
    }
}

} // namespace math
