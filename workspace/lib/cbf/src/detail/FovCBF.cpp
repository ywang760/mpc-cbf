#include <cbf/detail/FovCBF.h>
#include <vector>

// Default alpha function for Control Barrier Functions (CBF)
// Implements a linear scaling of the barrier function
// Parameters:
//   b: The barrier function value
//   gamma: Scaling factor to control convergence rate
GiNaC::ex defaultAlpha(GiNaC::ex b, double gamma = 1.0) {
    return gamma * b;
}

// Cubic alpha function for Control Barrier Functions
// Implements a cubic (power of 3) scaling of the barrier function
// Parameters:
//   myh: The barrier function value
//   mygamma: Scaling factor to control convergence rate
GiNaC::ex myAlpha(GiNaC::ex myh, double mygamma) {
    return mygamma * GiNaC::pow(myh, 3);
}

// Fifth-order alpha function for Control Barrier Functions
// Implements a fifth-order (power of 5) scaling of the barrier function
// Parameters:
//   myh: The barrier function value
//   mygamma: Scaling factor to control convergence rate
GiNaC::ex fifthAlpha(GiNaC::ex myh, double mygamma) {
    return mygamma * GiNaC::pow(myh, 5);
}

namespace cbf {
// Constructor for Field of View Control Barrier Function (FovCBF)
// Initializes the CBF with field of view constraints, safety distances, and velocity limits
// Parameters:
//   fov: Field of view angle in radians
//   safety_dist: Minimum safe distance to maintain
//   max_dist: Maximum sensing/control distance
//   vmin: Minimum velocity limits for each control dimension
//   vmax: Maximum velocity limits for each control dimension
FovCBF::FovCBF(double fov, double safety_dist, double max_dist, Eigen::VectorXd vmin,
               Eigen::VectorXd vmax)
    : fov(fov),
      Ds(safety_dist),
      Rs(max_dist),
      vmin(vmin),
      vmax(vmax),
      px("px"),
      py("py"),
      th("th"),
      vx("vx"),
      vy("vy"),
      w("w"),
      xt("xt"),
      yt("yt") {
    // Define dimensions of the state and control spaces
    STATE_VARS = 6;
    CONTROL_VARS = 3;
    gamma = 0.1; // Convergence rate parameter

    // Define symbolic variables for the state and control
    // px, py: position; th: heading; vx, vy: linear velocities; w: angular velocity
    // xt, yt: target position

    // System dynamics matrix (state transition matrix)
    // Represents continuous-time kinematics of the system
    A = {{0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 1},
         {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}};

    // Control input matrix
    // Maps control inputs to state derivatives
    B = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

    // Create state vector and transpose it for proper matrix operations
    state = {{px, py, th, vx, vy, w}};
    state = state.transpose();

    // Create target state vector and transpose it
    x_target = {{xt, yt}};
    x_target = x_target.transpose();

    // Define system dynamics (f + g*u)
    // f = A*x (drift term)
    // g = B (control input matrix)
    f = A.mul(state);
    g = B;

    // Set default alpha function to fifth-order for better performance
    alpha = fifthAlpha;

    // Initialize all Control Barrier Functions (CBFs)

    // Safety CBF: Ensures minimum safe distance is maintained
    auto res = initSafetyCBF();
    Ac_safe = res.first;  // Constraint matrix
    Bc_safe = res.second; // Constraint bound

    // Left border CBF: Constrains the target to be within the left FOV boundary
    auto res2 = initBorder1CBF();
    Ac_lb = res2.first;
    Bc_lb = res2.second;

    // Right border CBF: Constrains the target to be within the right FOV boundary
    auto res3 = initBorder2CBF();
    Ac_rb = res3.first;
    Bc_rb = res3.second;

    // Range CBF: Constrains the target to be within maximum range
    auto res4 = initRangeCBF();
    Ac_range = res4.first;
    Bc_range = res4.second;

    // Velocity constraint CBFs: Enforce maximum velocity limits
    // Setup max velocity constraint for x-direction
    GiNaC::ex bv1 = -state[3] + vmax(0);
    res = initVelCBF(bv1);
    Ac_v1_max = res.first;
    Bc_v1_max = res.second;

    // Setup max velocity constraint for y-direction
    GiNaC::ex bv2 = -state[4] + vmax(1);
    res = initVelCBF(bv2);
    Ac_v2_max = res.first;
    Bc_v2_max = res.second;

    // Setup max velocity constraint for angular velocity
    GiNaC::ex bv3 = -state[5] + vmax(2);
    res = initVelCBF(bv3);
    Ac_v3_max = res.first;
    Bc_v3_max = res.second;

    // Setup min velocity constraint for x-direction
    GiNaC::ex bv4 = state[3] - vmin(0);
    res = initVelCBF(bv4);
    Ac_v1_min = res.first;
    Bc_v1_min = res.second;

    // Setup min velocity constraint for y-direction
    GiNaC::ex bv5 = state[4] - vmin(1);
    res = initVelCBF(bv5);
    Ac_v2_min = res.first;
    Bc_v2_min = res.second;

    // Setup min velocity constraint for angular velocity
    GiNaC::ex bv6 = state[5] - vmin(2);
    res = initVelCBF(bv6);
    Ac_v3_min = res.first;
    Bc_v3_min = res.second;
}

// Destructor for the FovCBF class
FovCBF::~FovCBF() {
    std::cout << "Closing CBF ..." << std::endl;
}

// Initialize the safety Control Barrier Function
// This CBF ensures that the robot maintains a minimum safe distance (Ds) from the target
// Returns:
//   A pair containing the CBF constraint matrix (Ac) and bound (Bc)
std::pair<GiNaC::matrix, GiNaC::ex> FovCBF::initSafetyCBF() {
    // Calculate relative position of target in world frame
    GiNaC::matrix d2 = x_target.sub(GiNaC::matrix{{px}, {py}});

    // Rotation matrix to transform to robot frame
    GiNaC::matrix R = {{GiNaC::cos(th), GiNaC::sin(th)}, {-GiNaC::sin(th), GiNaC::cos(th)}};

    // Transform target position to robot frame
    GiNaC::matrix xt_rel = R.mul(d2);

    // Calculate squared distance to target
    GiNaC::ex norm2 = GiNaC::pow(xt_rel(0, 0), 2) + GiNaC::pow(xt_rel(1, 0), 2);

    // Define barrier function: h(x) = ||x - x_target||^2 - Ds^2
    // h > 0 when distance is greater than safety distance
    GiNaC::ex b1 = norm2 - GiNaC::pow(Ds, 2);

    // Calculate gradient of barrier function with respect to state variables
    GiNaC::matrix grad_b1 = GiNaC::matrix(STATE_VARS, 1);
    grad_b1(0, 0) = GiNaC::diff(b1, px);
    grad_b1(1, 0) = GiNaC::diff(b1, py);
    grad_b1(2, 0) = GiNaC::diff(b1, th);
    grad_b1(3, 0) = GiNaC::diff(b1, vx);
    grad_b1(4, 0) = GiNaC::diff(b1, vy);
    grad_b1(5, 0) = GiNaC::diff(b1, w);

    // Calculate Lie derivative of h along f: L_f h = ∇h · f
    GiNaC::ex lfb1 = 0.0;
    for (int i = 0; i < STATE_VARS; i++) {
        lfb1 = lfb1 + grad_b1(i, 0) * f(i, 0);
    }

    // Calculate gradient of L_f h
    GiNaC::matrix grad2_b1 = GiNaC::matrix(STATE_VARS, 1);
    grad2_b1(0, 0) = GiNaC::diff(lfb1, px);
    grad2_b1(1, 0) = GiNaC::diff(lfb1, py);
    grad2_b1(2, 0) = GiNaC::diff(lfb1, th);
    grad2_b1(3, 0) = GiNaC::diff(lfb1, vx);
    grad2_b1(4, 0) = GiNaC::diff(lfb1, vy);
    grad2_b1(5, 0) = GiNaC::diff(lfb1, w);

    // Calculate second Lie derivative: L_f^2 h = ∇(L_f h) · f
    GiNaC::ex lf2b1 = 0.0;
    for (int i = 0; i < STATE_VARS; i++) {
        lf2b1 = lf2b1 + grad2_b1(i, 0) * f(i, 0);
    }

    // Calculate gradient of alpha(h) where alpha is the class of K functions
    GiNaC::matrix grad_bc = GiNaC::matrix(STATE_VARS, 1);
    GiNaC::ex alpha_b = alpha(b1, gamma);
    grad_bc(0, 0) = GiNaC::diff(alpha_b, px);
    grad_bc(1, 0) = GiNaC::diff(alpha_b, py);
    grad_bc(2, 0) = GiNaC::diff(alpha_b, th);
    grad_bc(3, 0) = GiNaC::diff(alpha_b, vx);
    grad_bc(4, 0) = GiNaC::diff(alpha_b, vy);
    grad_bc(5, 0) = GiNaC::diff(alpha_b, w);

    // Calculate Lie derivative of alpha(h) along f
    GiNaC::ex lfb_c = 0.0;
    for (int i = 0; i < STATE_VARS; i++) {
        lfb_c = lfb_c + grad_bc(i, 0) * f(i, 0);
    }

    // Calculate L_g L_f h = ∇(L_f h) · g for each control input
    // This forms the constraint matrix Ac for the QP solver
    GiNaC::matrix Ac1 = GiNaC::matrix(1, CONTROL_VARS);
    for (int j = 0; j < CONTROL_VARS; j++) {
        GiNaC::ex Ac1j = 0.0;
        for (int i = 0; i < STATE_VARS; i++) {
            Ac1j += grad2_b1(i, 0) * g(i, j);
        }
        Ac1(0, j) = Ac1j;
    }

    // Calculate the constraint bound: L_f^2 h + L_f alpha(h) + alpha(L_f h + alpha(h))
    // This ensures the relative degree 2 CBF condition is satisfied
    GiNaC::ex B1 = lf2b1;
    GiNaC::ex B2 = lfb_c;
    GiNaC::ex B3 = alpha(lfb1 + alpha_b, gamma);
    GiNaC::ex Bc1 = B1 + B2 + B3;

    return std::make_pair(Ac1, Bc1);
}

// Initialize the left border Control Barrier Function
// This CBF ensures that the target remains within the left field of view boundary
// Returns:
//   A pair containing the CBF constraint matrix (Ac) and bound (Bc)
std::pair<GiNaC::matrix, GiNaC::ex> FovCBF::initBorder1CBF() {
    // Calculate relative position of target in world frame
    GiNaC::matrix d2 = x_target.sub(GiNaC::matrix{{px}, {py}});

    // Rotation matrix to transform to robot frame
    GiNaC::matrix R = {{GiNaC::cos(th), GiNaC::sin(th)}, {-GiNaC::sin(th), GiNaC::cos(th)}};

    // Transform target position to robot frame
    GiNaC::matrix xt_rel = R.mul(d2);

    // Define appropriate barrier function based on field of view angle
    GiNaC::ex b2;
    if (fov < M_PI) {
        // For FOV less than 180 degrees, use tangent-based boundary
        b2 = GiNaC::tan(fov / 2) * xt_rel(0, 0) + xt_rel(1, 0);
    } else if (fov == M_PI) {
        // For FOV of exactly 180 degrees, use x-coordinate as boundary
        b2 = xt_rel(0, 0);
    } else {
        // For FOV greater than 180 degrees, handle special cases
        if (py >= 0 || math::isApproximatelyEqual(fov, 2 * M_PI)) {
            // For full 360 degree FOV or when robot is above the target,
            // no constraint is needed (return max value)
            GiNaC::matrix Ac = GiNaC::matrix(1, CONTROL_VARS);
            for (int j = 0; j < CONTROL_VARS; j++) {
                Ac(0, j) = 0;
            }
            GiNaC::ex Bc = std::numeric_limits<double>::max();
            return std::make_pair(Ac, Bc);
        } else {
            // For FOV greater than 180 but less than 360 degrees
            b2 = GiNaC::tan((2 * M_PI - fov) / 2) * xt_rel(0, 0) - xt_rel(1, 0);
        }
    }

    // Calculate gradient of barrier function with respect to state variables
    GiNaC::matrix grad_b2 = GiNaC::matrix(STATE_VARS, 1);
    grad_b2(0, 0) = GiNaC::diff(b2, px);
    grad_b2(1, 0) = GiNaC::diff(b2, py);
    grad_b2(2, 0) = GiNaC::diff(b2, th);
    grad_b2(3, 0) = GiNaC::diff(b2, vx);
    grad_b2(4, 0) = GiNaC::diff(b2, vy);
    grad_b2(5, 0) = GiNaC::diff(b2, w);

    // Calculate Lie derivative of h along f: L_f h = ∇h · f
    GiNaC::ex lfb2 = 0.0;
    for (int i = 0; i < STATE_VARS; i++) {
        lfb2 = lfb2 + grad_b2(i, 0) * f(i, 0);
    }

    // Calculate gradient of L_f h
    GiNaC::matrix grad2_b2 = GiNaC::matrix(STATE_VARS, 1);
    grad2_b2(0, 0) = GiNaC::diff(lfb2, px);
    grad2_b2(1, 0) = GiNaC::diff(lfb2, py);
    grad2_b2(2, 0) = GiNaC::diff(lfb2, th);
    grad2_b2(3, 0) = GiNaC::diff(lfb2, vx);
    grad2_b2(4, 0) = GiNaC::diff(lfb2, vy);
    grad2_b2(5, 0) = GiNaC::diff(lfb2, w);

    // Calculate second Lie derivative: L_f^2 h = ∇(L_f h) · f
    GiNaC::ex lf2b2 = 0.0;
    for (int i = 0; i < STATE_VARS; i++) {
        lf2b2 = lf2b2 + grad2_b2(i, 0) * f(i, 0);
    }

    // Calculate gradient of alpha(h) where alpha is the class of K functions
    GiNaC::matrix grad_bc = GiNaC::matrix(STATE_VARS, 1);
    GiNaC::ex alpha_b2 = alpha(b2, gamma);
    grad_bc(0, 0) = GiNaC::diff(alpha_b2, px);
    grad_bc(1, 0) = GiNaC::diff(alpha_b2, py);
    grad_bc(2, 0) = GiNaC::diff(alpha_b2, th);
    grad_bc(3, 0) = GiNaC::diff(alpha_b2, vx);
    grad_bc(4, 0) = GiNaC::diff(alpha_b2, vy);
    grad_bc(5, 0) = GiNaC::diff(alpha_b2, w);

    // Calculate Lie derivative of alpha(h) along f
    GiNaC::ex lfb_c = 0.0;
    for (int i = 0; i < STATE_VARS; i++) {
        lfb_c = lfb_c + grad_bc(i, 0) * f(i, 0);
    }

    // Calculate L_g L_f h = ∇(L_f h) · g for each control input
    GiNaC::matrix Ac2 = GiNaC::matrix(1, CONTROL_VARS);
    for (int j = 0; j < CONTROL_VARS; j++) {
        GiNaC::ex Ac2j = 0.0;
        for (int i = 0; i < STATE_VARS; i++) {
            Ac2j += grad2_b2(i, 0) * g(i, j);
        }
        Ac2(0, j) = Ac2j;
    }

    // Calculate the constraint bound for the left FOV boundary
    GiNaC::ex B1 = lf2b2;
    GiNaC::ex B2 = lfb_c;
    GiNaC::ex B3 = alpha(lfb2 + alpha_b2, gamma);
    GiNaC::ex Bc2 = B1 + B2 + B3;

    return std::make_pair(Ac2, Bc2);
}

// Initialize the right border Control Barrier Function
// This CBF ensures that the target remains within the right field of view boundary
// Returns:
//   A pair containing the CBF constraint matrix (Ac) and bound (Bc)
std::pair<GiNaC::matrix, GiNaC::ex> FovCBF::initBorder2CBF() {
    // Calculate relative position of target in world frame
    GiNaC::matrix d2 = x_target.sub(GiNaC::matrix{{px}, {py}});

    // Rotation matrix to transform to robot frame
    GiNaC::matrix R = {{GiNaC::cos(th), GiNaC::sin(th)}, {-GiNaC::sin(th), GiNaC::cos(th)}};

    // Transform target position to robot frame
    GiNaC::matrix xt_rel = R.mul(d2);

    // Define appropriate barrier function based on field of view angle
    GiNaC::ex b3;
    if (fov < M_PI) {
        // For FOV less than 180 degrees, use tangent-based boundary
        b3 = GiNaC::tan(fov / 2) * xt_rel(0, 0) - xt_rel(1, 0);
    } else if (fov == M_PI) {
        // For FOV of exactly 180 degrees, use x-coordinate as boundary
        b3 = xt_rel(0, 0);
    } else {
        // For FOV greater than 180 degrees, handle special cases
        if (py < 0 || math::isApproximatelyEqual(fov, 2 * M_PI)) {
            // For full 360 degree FOV or when robot is below the target,
            // no constraint is needed (return max value)
            GiNaC::matrix Ac = GiNaC::matrix(1, CONTROL_VARS);
            for (int j = 0; j < CONTROL_VARS; j++) {
                Ac(0, j) = 0;
            }
            GiNaC::ex Bc = std::numeric_limits<double>::max();
            return std::make_pair(Ac, Bc);
        } else {
            // For FOV greater than 180 but less than 360 degrees
            b3 = GiNaC::tan((2 * M_PI - fov) / 2) * xt_rel(0, 0) + xt_rel(1, 0);
        }
    }

    // Calculate gradient of barrier function with respect to state variables
    GiNaC::matrix grad_b3 = GiNaC::matrix(STATE_VARS, 1);
    grad_b3(0, 0) = GiNaC::diff(b3, px);
    grad_b3(1, 0) = GiNaC::diff(b3, py);
    grad_b3(2, 0) = GiNaC::diff(b3, th);
    grad_b3(3, 0) = GiNaC::diff(b3, vx);
    grad_b3(4, 0) = GiNaC::diff(b3, vy);
    grad_b3(5, 0) = GiNaC::diff(b3, w);

    // Calculate Lie derivative of h along f: L_f h = ∇h · f
    GiNaC::ex lfb3 = 0.0;
    for (int i = 0; i < STATE_VARS; i++) {
        lfb3 = lfb3 + grad_b3(i, 0) * f(i, 0);
    }

    // Calculate gradient of L_f h
    GiNaC::matrix grad2_b3 = GiNaC::matrix(STATE_VARS, 1);
    grad2_b3(0, 0) = GiNaC::diff(lfb3, px);
    grad2_b3(1, 0) = GiNaC::diff(lfb3, py);
    grad2_b3(2, 0) = GiNaC::diff(lfb3, th);
    grad2_b3(3, 0) = GiNaC::diff(lfb3, vx);
    grad2_b3(4, 0) = GiNaC::diff(lfb3, vy);
    grad2_b3(5, 0) = GiNaC::diff(lfb3, w);

    // Calculate second Lie derivative: L_f^2 h = ∇(L_f h) · f
    GiNaC::ex lf2b3 = 0.0;
    for (int i = 0; i < STATE_VARS; i++) {
        lf2b3 = lf2b3 + grad2_b3(i, 0) * f(i, 0);
    }

    // Calculate gradient of alpha(h) where alpha is the class of K functions
    GiNaC::matrix grad_bc = GiNaC::matrix(STATE_VARS, 1);
    GiNaC::ex alpha_b3 = alpha(b3, gamma);
    grad_bc(0, 0) = GiNaC::diff(alpha_b3, px);
    grad_bc(1, 0) = GiNaC::diff(alpha_b3, py);
    grad_bc(2, 0) = GiNaC::diff(alpha_b3, th);
    grad_bc(3, 0) = GiNaC::diff(alpha_b3, vx);
    grad_bc(4, 0) = GiNaC::diff(alpha_b3, vy);
    grad_bc(5, 0) = GiNaC::diff(alpha_b3, w);

    // Calculate Lie derivative of alpha(h) along f
    GiNaC::ex lfb_c = 0.0;
    for (int i = 0; i < STATE_VARS; i++) {
        lfb_c = lfb_c + grad_bc(i, 0) * f(i, 0);
    }

    // Calculate L_g L_f h = ∇(L_f h) · g for each control input
    GiNaC::matrix Ac3 = GiNaC::matrix(1, CONTROL_VARS);
    for (int j = 0; j < CONTROL_VARS; j++) {
        GiNaC::ex Ac3j = 0.0;
        for (int i = 0; i < STATE_VARS; i++) {
            Ac3j += grad2_b3(i, 0) * g(i, j);
        }
        Ac3(0, j) = Ac3j;
    }

    // Calculate the constraint bound for the right FOV boundary
    GiNaC::ex B1 = lf2b3;
    GiNaC::ex B2 = lfb_c;
    GiNaC::ex B3 = alpha(lfb3 + alpha_b3, gamma);
    GiNaC::ex Bc3 = B1 + B2 + B3;

    return std::make_pair(Ac3, Bc3);
}

// Initialize the range Control Barrier Function
// This CBF ensures that the target remains within the maximum sensing/control range
// Returns:
//   A pair containing the CBF constraint matrix (Ac) and bound (Bc)
std::pair<GiNaC::matrix, GiNaC::ex> FovCBF::initRangeCBF() {
    // Calculate relative position of target in world frame
    GiNaC::matrix d2 = x_target.sub(GiNaC::matrix{{px}, {py}});

    // Rotation matrix to transform to robot frame
    GiNaC::matrix R = {{GiNaC::cos(th), GiNaC::sin(th)}, {-GiNaC::sin(th), GiNaC::cos(th)}};

    // Transform target position to robot frame
    GiNaC::matrix xt_rel = R.mul(d2);

    // Calculate squared distance to target
    GiNaC::ex norm2 = GiNaC::pow(xt_rel(0, 0), 2) + GiNaC::pow(xt_rel(1, 0), 2);

    // Define barrier function: h(x) = -||x - x_target||^2 + Rs^2
    // h > 0 when distance is less than maximum range
    GiNaC::ex b4 = -norm2 + GiNaC::pow(Rs, 2);

    // Calculate gradient of barrier function with respect to state variables
    GiNaC::matrix grad_b4 = GiNaC::matrix(STATE_VARS, 1);
    grad_b4(0, 0) = GiNaC::diff(b4, px);
    grad_b4(1, 0) = GiNaC::diff(b4, py);
    grad_b4(2, 0) = GiNaC::diff(b4, th);
    grad_b4(3, 0) = GiNaC::diff(b4, vx);
    grad_b4(4, 0) = GiNaC::diff(b4, vy);
    grad_b4(5, 0) = GiNaC::diff(b4, w);

    // Calculate Lie derivative of h along f: L_f h = ∇h · f
    GiNaC::ex lfb4 = 0.0;
    for (int i = 0; i < STATE_VARS; i++) {
        lfb4 = lfb4 + grad_b4(i, 0) * f(i, 0);
    }

    // Calculate gradient of L_f h
    GiNaC::matrix grad2_b4 = GiNaC::matrix(STATE_VARS, 1);
    grad2_b4(0, 0) = GiNaC::diff(lfb4, px);
    grad2_b4(1, 0) = GiNaC::diff(lfb4, py);
    grad2_b4(2, 0) = GiNaC::diff(lfb4, th);
    grad2_b4(3, 0) = GiNaC::diff(lfb4, vx);
    grad2_b4(4, 0) = GiNaC::diff(lfb4, vy);
    grad2_b4(5, 0) = GiNaC::diff(lfb4, w);

    // Calculate second Lie derivative: L_f^2 h = ∇(L_f h) · f
    GiNaC::ex lf2b4 = 0.0;
    for (int i = 0; i < STATE_VARS; i++) {
        lf2b4 = lf2b4 + grad2_b4(i, 0) * f(i, 0);
    }

    // Calculate gradient of alpha(h) where alpha is the class of K functions
    GiNaC::matrix grad_bc = GiNaC::matrix(STATE_VARS, 1);
    GiNaC::ex alpha_b4 = alpha(b4, gamma);
    grad_bc(0, 0) = GiNaC::diff(alpha_b4, px);
    grad_bc(1, 0) = GiNaC::diff(alpha_b4, py);
    grad_bc(2, 0) = GiNaC::diff(alpha_b4, th);
    grad_bc(3, 0) = GiNaC::diff(alpha_b4, vx);
    grad_bc(4, 0) = GiNaC::diff(alpha_b4, vy);
    grad_bc(5, 0) = GiNaC::diff(alpha_b4, w);

    // Calculate Lie derivative of alpha(h) along f
    GiNaC::ex lfb_c = 0.0;
    for (int i = 0; i < STATE_VARS; i++) {
        lfb_c = lfb_c + grad_bc(i, 0) * f(i, 0);
    }

    // Calculate L_g L_f h = ∇(L_f h) · g for each control input
    GiNaC::matrix Ac4 = GiNaC::matrix(1, CONTROL_VARS);
    for (int j = 0; j < CONTROL_VARS; j++) {
        GiNaC::ex Ac4j = 0.0;
        for (int i = 0; i < STATE_VARS; i++) {
            Ac4j += grad2_b4(i, 0) * g(i, j);
        }
        Ac4(0, j) = Ac4j;
    }

    // Calculate the constraint bound for the range constraint
    GiNaC::ex B1 = lf2b4;
    GiNaC::ex B2 = lfb_c;
    GiNaC::ex B3 = alpha(lfb4 + alpha_b4, gamma);
    GiNaC::ex Bc4 = B1 + B2 + B3;

    return std::make_pair(Ac4, Bc4);
}

// Initialize velocity Control Barrier Functions
// This creates CBFs to enforce velocity limits (min or max)
// Parameters:
//   bv: The barrier function for the velocity constraint
// Returns:
//   A pair containing the CBF constraint matrix (Ac) and bound (Bc)
std::pair<GiNaC::matrix, GiNaC::ex> FovCBF::initVelCBF(GiNaC::ex bv) {
    // Calculate gradient of barrier function with respect to state variables
    GiNaC::matrix grad_bv = GiNaC::matrix(STATE_VARS, 1);
    grad_bv(0, 0) = GiNaC::diff(bv, px);
    grad_bv(1, 0) = GiNaC::diff(bv, py);
    grad_bv(2, 0) = GiNaC::diff(bv, th);
    grad_bv(3, 0) = GiNaC::diff(bv, vx);
    grad_bv(4, 0) = GiNaC::diff(bv, vy);
    grad_bv(5, 0) = GiNaC::diff(bv, w);

    // Calculate Lie derivative of h along f: L_f h = ∇h · f
    GiNaC::ex lfbv = 0.0;
    for (int i = 0; i < STATE_VARS; i++) {
        lfbv = lfbv + grad_bv(i, 0) * f(i, 0);
    }

    // Calculate L_g h = ∇h · g for each control input
    // For velocity constraints, this is simpler as they are relative degree 1
    GiNaC::matrix Ac_v1 = GiNaC::matrix(1, CONTROL_VARS);
    for (int j = 0; j < CONTROL_VARS; j++) {
        GiNaC::ex Ac_v1j = 0.0;
        for (int i = 0; i < STATE_VARS; i++) {
            Ac_v1j += grad_bv(i, 0) * g(i, j);
        }
        Ac_v1(0, j) = Ac_v1j;
    }

    // Calculate the constraint bound: L_f h + α(h)
    // Using default linear alpha for velocity constraints
    GiNaC::ex Bc_v1 = lfbv + defaultAlpha(bv, 1);
    return std::make_pair(Ac_v1, Bc_v1);
}

// Substitute state and target values into a symbolic matrix
// Parameters:
//   a: Symbolic matrix
//   state: Current state vector
//   target_state: Target state vector
// Returns:
//   Matrix with numerical values substituted
GiNaC::ex FovCBF::matrixSubs(GiNaC::matrix a, Eigen::VectorXd state, Eigen::VectorXd target_state) {
    // Substitute each state variable with its numerical value
    GiNaC::ex tmp = GiNaC::subs(a, px == state(0));
    tmp = GiNaC::subs(tmp, py == state(1));
    tmp = GiNaC::subs(tmp, th == state(2));
    tmp = GiNaC::subs(tmp, vx == state(3));
    tmp = GiNaC::subs(tmp, vy == state(4));
    tmp = GiNaC::subs(tmp, w == state(5));
    tmp = GiNaC::subs(tmp, xt == target_state(0));
    tmp = GiNaC::subs(tmp, yt == target_state(1));
    return tmp;
}

// Substitute state and target values into a symbolic expression
// Parameters:
//   a: Symbolic expression
//   state: Current state vector
//   target_state: Target state vector
// Returns:
//   Expression with numerical values substituted
GiNaC::ex FovCBF::valueSubs(GiNaC::ex a, Eigen::VectorXd state, Eigen::VectorXd target_state) {
    // Substitute each state variable with its numerical value
    GiNaC::ex tmp = GiNaC::subs(a, px == state(0));
    tmp = GiNaC::subs(tmp, py == state(1));
    tmp = GiNaC::subs(tmp, th == state(2));
    tmp = GiNaC::subs(tmp, vx == state(3));
    tmp = GiNaC::subs(tmp, vy == state(4));
    tmp = GiNaC::subs(tmp, w == state(5));
    tmp = GiNaC::subs(tmp, xt == target_state(0));
    tmp = GiNaC::subs(tmp, yt == target_state(1));
    return tmp;
}

// Get the safety constraint vector for the current state and target
// Parameters:
//   state: Current state vector
//   target_state: Target state vector
// Returns:
//   Constraint vector for the QP solver
Eigen::VectorXd FovCBF::getSafetyConstraints(Eigen::VectorXd state, Eigen::VectorXd target_state) {
    // Substitute numerical values into symbolic matrix and convert to Eigen vector
    GiNaC::ex matrix_expr = matrixSubs(Ac_safe, state, target_state);
    Eigen::VectorXd Ac;
    Ac.resize(CONTROL_VARS);
    Ac.setZero();
    for (int i = 0; i < CONTROL_VARS; i++) {
        GiNaC::ex val = matrix_expr[i].evalf();
        Ac(i) = GiNaC::ex_to<GiNaC::numeric>(val).to_double();
    }

    return Ac;
}

// Get the left boundary constraint vector for the current state and target
// Parameters:
//   state: Current state vector
//   target_state: Target state vector
// Returns:
//   Constraint vector for the QP solver
Eigen::VectorXd FovCBF::getLBConstraints(Eigen::VectorXd state, Eigen::VectorXd target_state) {
    // Substitute numerical values into symbolic matrix and convert to Eigen vector
    GiNaC::ex matrix_expr = matrixSubs(Ac_lb, state, target_state);
    Eigen::VectorXd Ac;
    Ac.resize(CONTROL_VARS);
    Ac.setZero();
    for (int i = 0; i < CONTROL_VARS; i++) {
        GiNaC::ex val = matrix_expr[i].evalf();
        Ac(i) = GiNaC::ex_to<GiNaC::numeric>(val).to_double();
    }

    return Ac;
}

// Get the right boundary constraint vector for the current state and target
// Parameters:
//   state: Current state vector
//   target_state: Target state vector
// Returns:
//   Constraint vector for the QP solver
Eigen::VectorXd FovCBF::getRBConstraints(Eigen::VectorXd state, Eigen::VectorXd target_state) {
    // Substitute numerical values into symbolic matrix and convert to Eigen vector
    GiNaC::ex matrix_expr = matrixSubs(Ac_rb, state, target_state);
    Eigen::VectorXd Ac;
    Ac.resize(CONTROL_VARS);
    Ac.setZero();
    for (int i = 0; i < CONTROL_VARS; i++) {
        GiNaC::ex val = matrix_expr[i].evalf();
        Ac(i) = GiNaC::ex_to<GiNaC::numeric>(val).to_double();
    }

    return Ac;
}

// Get the range constraint vector for the current state and target
// Parameters:
//   state: Current state vector
//   target_state: Target state vector
// Returns:
//   Constraint vector for the QP solver
Eigen::VectorXd FovCBF::getRangeConstraints(Eigen::VectorXd state, Eigen::VectorXd target_state) {
    // Substitute numerical values into symbolic matrix and convert to Eigen vector
    GiNaC::ex matrix_expr = matrixSubs(Ac_range, state, target_state);
    Eigen::VectorXd Ac;
    Ac.resize(CONTROL_VARS);
    Ac.setZero();
    for (int i = 0; i < CONTROL_VARS; i++) {
        GiNaC::ex val = matrix_expr[i].evalf();
        Ac(i) = GiNaC::ex_to<GiNaC::numeric>(val).to_double();
    }

    return Ac;
}

// Get the maximum velocity constraint matrix for the current state
// Parameters:
//   state: Current state vector
// Returns:
//   Constraint matrix for the QP solver (one row per control variable)
Eigen::MatrixXd FovCBF::getMaxVelContraints(Eigen::VectorXd state) {
    // Create dummy target (not needed for velocity constraints)
    Eigen::Vector2d dummy_target;
    dummy_target.setZero();

    // Get constraint matrices for each velocity component
    std::vector<GiNaC::ex> expressions(CONTROL_VARS);
    expressions[0] = matrixSubs(Ac_v1_max, state, dummy_target);
    expressions[1] = matrixSubs(Ac_v2_max, state, dummy_target);
    expressions[2] = matrixSubs(Ac_v3_max, state, dummy_target);

    // Convert to Eigen matrix
    Eigen::MatrixXd Acs;
    Acs.resize(CONTROL_VARS, CONTROL_VARS);
    Acs.setZero();
    for (int i = 0; i < CONTROL_VARS; ++i) {
        for (int j = 0; j < CONTROL_VARS; ++j) {
            GiNaC::ex val = expressions[i][j].evalf();
            Acs(i, j) = GiNaC::ex_to<GiNaC::numeric>(val).to_double();
        }
    }

    return Acs;
}

// Get the minimum velocity constraint matrix for the current state
// Parameters:
//   state: Current state vector
// Returns:
//   Constraint matrix for the QP solver (one row per control variable)
Eigen::MatrixXd FovCBF::getMinVelContraints(Eigen::VectorXd state) {
    // Create dummy target (not needed for velocity constraints)
    Eigen::Vector2d dummy_target;
    dummy_target.setZero();

    // Get constraint matrices for each velocity component
    std::vector<GiNaC::ex> expressions(CONTROL_VARS);
    expressions[0] = matrixSubs(Ac_v1_min, state, dummy_target);
    expressions[1] = matrixSubs(Ac_v2_min, state, dummy_target);
    expressions[2] = matrixSubs(Ac_v3_min, state, dummy_target);

    // Convert to Eigen matrix
    Eigen::MatrixXd Acs;
    Acs.resize(CONTROL_VARS, CONTROL_VARS);
    Acs.setZero();
    for (int i = 0; i < CONTROL_VARS; ++i) {
        for (int j = 0; j < CONTROL_VARS; ++j) {
            GiNaC::ex val = expressions[i][j].evalf();
            Acs(i, j) = GiNaC::ex_to<GiNaC::numeric>(val).to_double();
        }
    }

    return Acs;
}

// Get the safety constraint bound for the current state and target
// Parameters:
//   state: Current state vector
//   target_state: Target state vector
// Returns:
//   Bound value for the safety constraint
double FovCBF::getSafetyBound(Eigen::VectorXd state, Eigen::VectorXd target_state) {
    // Substitute numerical values and evaluate
    GiNaC::ex expr = valueSubs(Bc_safe, state, target_state);
    double Bc = GiNaC::ex_to<GiNaC::numeric>(expr).to_double();

    return Bc;
}

// Get the range constraint bound for the current state and target
// Parameters:
//   state: Current state vector
//   target_state: Target state vector
// Returns:
//   Bound value for the range constraint
double FovCBF::getRangeBound(Eigen::VectorXd state, Eigen::VectorXd target_state) {
    // Substitute numerical values and evaluate
    GiNaC::ex expr = valueSubs(Bc_range, state, target_state);
    double Bc = GiNaC::ex_to<GiNaC::numeric>(expr).to_double();

    return Bc;
}

// Get the left boundary constraint bound for the current state and target
// Parameters:
//   state: Current state vector
//   target_state: Target state vector
// Returns:
//   Bound value for the left boundary constraint
double FovCBF::getLBBound(Eigen::VectorXd state, Eigen::VectorXd target_state) {
    // Substitute numerical values and evaluate
    GiNaC::ex expr = valueSubs(Bc_lb, state, target_state);
    double Bc = GiNaC::ex_to<GiNaC::numeric>(expr).to_double();

    return Bc;
}

// Get the right boundary constraint bound for the current state and target
// Parameters:
//   state: Current state vector
//   target_state: Target state vector
// Returns:
//   Bound value for the right boundary constraint
double FovCBF::getRBBound(Eigen::VectorXd state, Eigen::VectorXd target_state) {
    // Substitute numerical values and evaluate
    GiNaC::ex expr = valueSubs(Bc_rb, state, target_state);
    double Bc = GiNaC::ex_to<GiNaC::numeric>(expr).to_double();

    return Bc;
}

// Get the maximum velocity constraint bounds for the current state
// Parameters:
//   state: Current state vector
// Returns:
//   Vector of bound values for maximum velocity constraints
Eigen::VectorXd FovCBF::getMaxVelBounds(Eigen::VectorXd state) {
    // Create dummy target (not needed for velocity constraints)
    Eigen::Vector2d dummy_target;
    dummy_target.setZero();

    // Get bounds for each velocity component
    std::vector<GiNaC::ex> expressions(CONTROL_VARS);
    expressions[0] = valueSubs(Bc_v1_max, state, dummy_target);
    expressions[1] = valueSubs(Bc_v2_max, state, dummy_target);
    expressions[2] = valueSubs(Bc_v3_max, state, dummy_target);

    // Convert to Eigen vector
    Eigen::VectorXd Bs;
    Bs.resize(CONTROL_VARS);

    for (int i = 0; i < CONTROL_VARS; ++i) {
        Bs(i) = GiNaC::ex_to<GiNaC::numeric>(expressions[i]).to_double();
    }

    return Bs;
}

// Get the minimum velocity constraint bounds for the current state
// Parameters:
//   state: Current state vector
// Returns:
//   Vector of bound values for minimum velocity constraints
Eigen::VectorXd FovCBF::getMinVelBounds(Eigen::VectorXd state) {
    // Create dummy target (not needed for velocity constraints)
    Eigen::Vector2d dummy_target;
    dummy_target.setZero();

    // Get bounds for each velocity component
    std::vector<GiNaC::ex> expressions(CONTROL_VARS);
    expressions[0] = valueSubs(Bc_v1_min, state, dummy_target);
    expressions[1] = valueSubs(Bc_v2_min, state, dummy_target);
    expressions[2] = valueSubs(Bc_v3_min, state, dummy_target);

    // Convert to Eigen vector
    Eigen::VectorXd Bs;
    Bs.resize(CONTROL_VARS);

    for (int i = 0; i < CONTROL_VARS; ++i) {
        Bs(i) = GiNaC::ex_to<GiNaC::numeric>(expressions[i]).to_double();
    }

    return Bs;
}

// Set a custom alpha function for the CBFs
// Parameters:
//   newAlpha: New alpha function to use
void FovCBF::setAlpha(std::function<GiNaC::ex(GiNaC::ex, double)> newAlpha) {
    alpha = newAlpha;
}

} // namespace cbf
