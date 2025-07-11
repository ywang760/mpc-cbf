#include <vector>
#include <cbf/detail/ConnectivityCBF.h>
#include <cbf/Helpers.hpp>

// Default alpha function for Control Barrier Functions (CBF)
// Implements a linear scaling of the barrier function
// Parameters:
//   b: The barrier function value
//   gamma: Scaling factor to control convergence rate
GiNaC::ex defaultAlpha(GiNaC::ex b, double gamma = 1.0)
{
    return gamma * b;
}

// Cubic alpha function for Control Barrier Functions
// Implements a cubic (power of 3) scaling of the barrier function
// Parameters:
//   myh: The barrier function value
//   mygamma: Scaling factor to control convergence rate
GiNaC::ex myAlpha(GiNaC::ex myh, double mygamma)
{
    return mygamma * GiNaC::pow(myh, 3);
}

// Fifth-order alpha function for Control Barrier Functions
// Implements a fifth-order (power of 5) scaling of the barrier function
// Parameters:
//   myh: The barrier function value
//   mygamma: Scaling factor to control convergence rate
GiNaC::ex fifthAlpha(GiNaC::ex myh, double mygamma)
{
    return mygamma * GiNaC::pow(myh, 5);
}

namespace cbf
{
    auto logger = spdlog::default_logger();

    // Constructor for Connectivity Control Barrier Function
    // Initializes the CBF with minimum/maximum distance constraints and velocity limits
    // Parameters:
    //   min_dist: Minimum allowed distance between agents
    //   max_dist: Maximum allowed distance between agents (connectivity range)
    //   vmin: Minimum velocity limits for each control dimension
    //   vmax: Maximum velocity limits for each control dimension
    ConnectivityCBF::ConnectivityCBF(double min_dist, double max_dist, Eigen::VectorXd vmin, Eigen::VectorXd vmax)
        : dmin(min_dist), dmax(max_dist), vmin(vmin), vmax(vmax),
          px("px"), py("py"), th("th"), vx("vx"), vy("vy"), w("w"), px_n("px_n"), py_n("py_n"), vx_n("vx_n"), vy_n("vy_n"), h("h")
    {
        // Define dimensions of the state and control spaces
        STATE_VARS = 6;
        CONTROL_VARS = 3;
        gamma = 5.0;   // More aggressive convergence rate (increased from 1.0)
        epsilon = 0.1; // lambda2_min for connectivity CBF

        // System dynamics matrix (state transition matrix)
        // Represents continuous-time kinematics of the system
        A = {{0, 0, 0, 1, 0, 0},
             {0, 0, 0, 0, 1, 0},
             {0, 0, 0, 0, 0, 1},
             {0, 0, 0, 0, 0, 0},
             {0, 0, 0, 0, 0, 0},
             {0, 0, 0, 0, 0, 0}};

        // Control input matrix
        // Maps control inputs to state derivatives
        B = {{0, 0, 0},
             {0, 0, 0},
             {0, 0, 0},
             {1, 0, 0},
             {0, 1, 0},
             {0, 0, 1}};

        // Create state vector of the ego agent
        state = {{px, py, th, vx, vy, w}};
        state = state.transpose();

        // Create state vector of the neighbor agent
        p_n = {{px_n, py_n}};
        p_n = p_n.transpose();
        v_n = {{vx_n, vy_n}};
        v_n = v_n.transpose();

        // Define system dynamics (f + g*u)
        // f = A*x (drift term)
        // g = B (control input matrix)
        f = A.mul(state);
        g = B;

        // Set default alpha function to cubic for better performance
        alpha = myAlpha;

        // Initialize all Control Barrier Functions (CBFs)

        // Safety CBF: Ensures minimum safe distance is maintained between agents
        auto res = initSafetyCBF();
        Ac_safe = res.first;  // Constraint matrix
        Bc_safe = res.second; // Constraint bound

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

    // Destructor for the ConnectivityCBF class
    ConnectivityCBF::~ConnectivityCBF()
    {
        logger->info("Closing ConnectivityCBF instance");
    }

    // Initialize the minimum distance Control Barrier Function
    // This CBF ensures that the robot maintains a minimum safe distance (dmin) from other agents
    // Returns:
    //   A pair containing the CBF constraint matrix (Ac) and bound (Bc)
    std::pair<GiNaC::matrix, GiNaC::ex> ConnectivityCBF::initSafetyCBF()
    {

        GiNaC::ex dx = px - px_n;  // x-distance to other agent
        GiNaC::ex dy = py - py_n;  // y-distance to other agent
        GiNaC::ex dvx = vx - vx_n; // x-velocity difference to other agent
        GiNaC::ex dvy = vy - vy_n; // y-velocity difference

        // Calculate barrier function: h(p) = ||p - p_n||^2 - dmin^2
        GiNaC::ex h = dx * dx + dy * dy - GiNaC::pow(dmin, 2);

        // First Lie derivative of the barrier function: L_f h = ∇h · f
        // This could be analytically derived
        GiNaC::ex Lf_h = 2 * (dx * dvx + dy * dvy);

        // Second Lie derivative of the barrier function: L_f^2 h = ∇(L_f h) · f
        // This could also be analytically derived
        GiNaC::ex Lf2_h = 2 * (dvx * dvx + dvy * dvy); // L_f^2 h = 2 * (dvx*dvx + dvy*dvy)

        // Calculate gradient of alpha(h) where alpha is the class of K functions
        GiNaC::matrix grad_alpha = GiNaC::matrix(STATE_VARS, 1);
        GiNaC::ex alpha_b = alpha(h, gamma);
        grad_alpha(0, 0) = GiNaC::diff(alpha_b, px);
        grad_alpha(1, 0) = GiNaC::diff(alpha_b, py);
        grad_alpha(2, 0) = GiNaC::diff(alpha_b, th);
        grad_alpha(3, 0) = GiNaC::diff(alpha_b, vx);
        grad_alpha(4, 0) = GiNaC::diff(alpha_b, vy);
        grad_alpha(5, 0) = GiNaC::diff(alpha_b, w);

        // Calculate Lie derivative of alpha(h) along f
        GiNaC::ex Lf_alpha = 0.0;
        for (int i = 0; i < STATE_VARS; i++)
        {
            Lf_alpha = Lf_alpha + grad_alpha(i, 0) * f(i, 0);
        }

        // Calculate Ac = L_g L_f h = ∇(L_f h) · g
        // This can be analytically derived as well
        // TODO: expand to account for DIM=3
        GiNaC::matrix Ac = GiNaC::matrix(1, CONTROL_VARS);
        Ac(0, 0) = 2 * dx;
        Ac(0, 1) = 2 * dy;

        // Calculate the constraint bound: L_f^2 h + L_f alpha(h) + alpha(L_f h + alpha(h))
        GiNaC::ex psi1 = Lf_h + alpha_b;
        GiNaC::ex Bc = Lf2_h + Lf_alpha + alpha(psi1, gamma);

        return std::make_pair(Ac, Bc);
    }

    // Initialize velocity Control Barrier Functions
    // This creates CBFs to enforce velocity limits (min or max)
    // Parameters:
    //   bv: The barrier function for the velocity constraint
    // Returns:
    //   A pair containing the CBF constraint matrix (Ac) and bound (Bc)
    std::pair<GiNaC::matrix, GiNaC::ex> ConnectivityCBF::initVelCBF(GiNaC::ex bv)
    {
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
        for (int i = 0; i < STATE_VARS; i++)
        {
            lfbv = lfbv + grad_bv(i, 0) * f(i, 0);
        }

        // Calculate L_g h = ∇h · g for each control input
        // For velocity constraints, this is simpler as they are relative degree 1
        GiNaC::matrix Ac_v1 = GiNaC::matrix(1, CONTROL_VARS);
        for (int j = 0; j < CONTROL_VARS; j++)
        {
            GiNaC::ex Ac_v1j = 0.0;
            for (int i = 0; i < STATE_VARS; i++)
            {
                Ac_v1j += grad_bv(i, 0) * g(i, j);
            }
            Ac_v1(0, j) = Ac_v1j;
        }

        // Calculate the constraint bound: L_f h + α(h)
        // Using default linear alpha for velocity constraints
        GiNaC::ex Bc_v1 = lfbv + defaultAlpha(bv, 1);
        return std::make_pair(Ac_v1, Bc_v1);
    }

    // Get the minimum distance constraint vector for the current state and agent
    // Parameters:
    //   state: Current state vector
    //   neighbor_state: Other agent state vector
    // Returns:
    //   Constraint vector for the QP solver
    Eigen::VectorXd ConnectivityCBF::getSafetyConstraints(Eigen::VectorXd state, Eigen::VectorXd neighbor_state)
    {
        // Substitute numerical values into symbolic matrix and convert to Eigen vector
        GiNaC::ex matrix_expr = matrixSubs(Ac_safe, state, neighbor_state, *this);
        Eigen::VectorXd Ac;
        Ac.resize(CONTROL_VARS);
        Ac.setZero();
        for (int i = 0; i < CONTROL_VARS; i++)
        {
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
    Eigen::MatrixXd ConnectivityCBF::getMaxVelContraints(Eigen::VectorXd state)
    {
        // Create dummy agent (not needed for velocity constraints)
        Eigen::VectorXd dummy_agent(6);
        dummy_agent.setZero();

        // Get constraint matrices for each velocity component
        std::vector<GiNaC::ex> expressions(CONTROL_VARS);
        expressions[0] = matrixSubs(Ac_v1_max, state, dummy_agent, *this);
        expressions[1] = matrixSubs(Ac_v2_max, state, dummy_agent, *this);
        expressions[2] = matrixSubs(Ac_v3_max, state, dummy_agent, *this);

        // Convert to Eigen matrix
        Eigen::MatrixXd Acs;
        Acs.resize(CONTROL_VARS, CONTROL_VARS);
        Acs.setZero();
        for (int i = 0; i < CONTROL_VARS; ++i)
        {
            for (int j = 0; j < CONTROL_VARS; ++j)
            {
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
    Eigen::MatrixXd ConnectivityCBF::getMinVelContraints(Eigen::VectorXd state)
    {
        // Create dummy agent (not needed for velocity constraints)
        Eigen::VectorXd dummy_agent(6);
        dummy_agent.setZero();

        // Get constraint matrices for each velocity component
        std::vector<GiNaC::ex> expressions(CONTROL_VARS);
        expressions[0] = matrixSubs(Ac_v1_min, state, dummy_agent, *this);
        expressions[1] = matrixSubs(Ac_v2_min, state, dummy_agent, *this);
        expressions[2] = matrixSubs(Ac_v3_min, state, dummy_agent, *this);

        // Convert to Eigen matrix
        Eigen::MatrixXd Acs;
        Acs.resize(CONTROL_VARS, CONTROL_VARS);
        Acs.setZero();
        for (int i = 0; i < CONTROL_VARS; ++i)
        {
            for (int j = 0; j < CONTROL_VARS; ++j)
            {
                GiNaC::ex val = expressions[i][j].evalf();
                Acs(i, j) = GiNaC::ex_to<GiNaC::numeric>(val).to_double();
            }
        }

        return Acs;
    }

    // Get sigma_value for the lambda2 function
    // By default, sigma is set to the fourth power of dmax divided by log(2)
    double ConnectivityCBF::getSigma() const
    {
        return dmax * dmax * dmax * dmax / std::log(2.0);
    }

    // Given the numerical robot positions, compute the second smallest eigenvalue (lambda2) for the Laplacian matrix
    // Parameters:
    //   robot_positions: Matrix of robot positions
    // Returns:
    //   A pair containing the second smallest eigenvalue (lambda2) and the corresponding eigenvector
    std::pair<double, Eigen::VectorXd> ConnectivityCBF::getLambda2(
        const Eigen::MatrixXd &robot_positions)
    {
        const int N = robot_positions.rows();
        // Construct the Adjacency matrix A in numerical form
        Eigen::MatrixXd A_num = Eigen::MatrixXd::Zero(N, N);
        double Rs_value2 = std::pow(dmax, 2); // Rs^2
        double sigma_value = getSigma();      // Get sigma value for the lambda2 function
        for (int i = 0; i < N; ++i)
        {
            const auto pi = robot_positions.row(i);
            for (int j = 0; j < N; ++j)
            {
                if (i == j)
                    continue;
                double dij2 = (pi - robot_positions.row(j)).squaredNorm();

                if (dij2 <= Rs_value2)
                {
                    double weight = std::exp(std::pow(Rs_value2 - dij2, 2) / sigma_value) - 1;
                    A_num(i, j) = weight;
                }
            }
        }

        // Construct the Diagonal matrix D and Laplacian matrix L in numerical form
        Eigen::MatrixXd D_num = Eigen::MatrixXd::Zero(N, N);
        for (int i = 0; i < N; ++i)
            D_num(i, i) = A_num.row(i).sum();
        Eigen::MatrixXd L_num = D_num - A_num;

        // Numerically solve for the second smallest eigenvalue and corresponding eigenvector
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(L_num);
        Eigen::VectorXd eigenvals = solver.eigenvalues();
        Eigen::MatrixXd eigenvecs = solver.eigenvectors();

        Eigen::VectorXd eigenvec = eigenvecs.col(1); // Second smallest eigenvector
        eigenvec.normalize();                        // Normalize the eigenvector
        return {eigenvals(1), eigenvec};
    }

    void ConnectivityCBF::initSymbolLists(int N)
    {
        px_list.clear();
        py_list.clear();
        eigenvec_list.clear();
        for (int i = 0; i < N; ++i)
        {
            px_list.emplace_back(GiNaC::symbol("px" + std::to_string(i)));
            py_list.emplace_back(GiNaC::symbol("py" + std::to_string(i)));
            eigenvec_list.emplace_back(GiNaC::symbol("eigenvec" + std::to_string(i)));
        }
    }

    // Parameters:
    //   N: Number of agents
    //   Rs: Symbolic representation of the maximum distance for connectivity
    //   sigma: Symbolic representation of the weight function parameter
    GiNaC::matrix ConnectivityCBF::compute_full_grad_h(int N, const GiNaC::ex &Rs, const GiNaC::ex &sigma)
    {
        initSymbolLists(N);
        GiNaC::matrix grad(N, 2); // TODO: for 3d: this needs to be N x 3 (figure out a way to generalize this)
        for (int i = 0; i < N; ++i)
        {
            GiNaC::ex dLdx = 0, dLdy = 0;
            for (int j = 0; j < N; ++j)
            {
                if (i == j)
                    continue;
                GiNaC::ex dx = px_list[i] - px_list[j];
                GiNaC::ex dy = py_list[i] - py_list[j];
                GiNaC::ex dij2 = dx * dx + dy * dy;
                GiNaC::ex diff = GiNaC::pow(Rs, 2) - dij2; // TODO: here check diff is positive?
                GiNaC::ex Aij = GiNaC::exp(GiNaC::pow(diff, 2) / sigma) - 1;

                // Explicit closed-form expressions for the derivatives
                GiNaC::ex dAij_dx = -4 * (Aij + 1) * diff / sigma * dx;
                GiNaC::ex dAij_dy = -4 * (Aij + 1) * diff / sigma * dy;

                // Calculate each component in the gradient (formula (12) in the paper)
                GiNaC::ex vdiff2 = GiNaC::pow(eigenvec_list[i] - eigenvec_list[j], 2);
                dLdx += dAij_dx * vdiff2;
                dLdy += dAij_dy * vdiff2;
            }
            grad(i, 0) = dLdx;
            grad(i, 1) = dLdy;
        }
        return grad;
    }

    std::pair<GiNaC::matrix, GiNaC::ex> ConnectivityCBF::initConnCBF(
        int N, int self_idx)
    {
        // ∇h (gradient of h), shape = N×2
        const double sigma_val = getSigma();
        GiNaC::matrix full_grad_sym = compute_full_grad_h(N, dmax, sigma_val); // Full gradient, shape Nx2

        // Gradient for the current robot, shape 1x2
        GiNaC::matrix grad_h(1, CONTROL_VARS);
        grad_h(0, 0) = full_grad_sym(self_idx, 0); // ∂h/∂p_x  (symbolic)
        grad_h(0, 1) = full_grad_sym(self_idx, 1); // ∂h/∂p_y  (symbolic)
        grad_h(0, 2) = 0;                          // ∂h/∂θ (not used in connectivity CBF)

        // L_f h = ∇h · f
        // Only the terms related to vx and vy are non-zero
        GiNaC::ex lfh_sym = GiNaC::ex(grad_h(0, 0)) * vx + GiNaC::ex(grad_h(0, 1)) * vy;

        // Lf² h = ∇(L_f h) · f
        // Compute the Hessian matrix
        GiNaC::matrix hess_sym(2, 2);
        hess_sym(0, 0) = GiNaC::diff(grad_h(0, 0), px_list[self_idx]);
        hess_sym(0, 1) = GiNaC::diff(grad_h(0, 0), py_list[self_idx]);
        hess_sym(1, 0) = GiNaC::diff(grad_h(0, 1), px_list[self_idx]);
        hess_sym(1, 1) = GiNaC::diff(grad_h(0, 1), py_list[self_idx]);
        GiNaC::ex lf2h_sym = vx * (hess_sym(0, 0) * vx + hess_sym(0, 1) * vy) +
                             vy * (hess_sym(1, 0) * vx + hess_sym(1, 1) * vy);

        setAlpha(defaultAlpha);

        // Bc = L_f² h + L_f(α(h)) + α(L_f h + α(h))
        // Calculate the first Lie derivative of alpha(h)
        // TODO: this version is not quite working

        // GiNaC::matrix grad_alpha = GiNaC::matrix(STATE_VARS, 1);
        // GiNaC::ex alpha_h = alpha(h, gamma);
        // grad_alpha(0, 0) = GiNaC::diff(alpha_h, px);
        // grad_alpha(1, 0) = GiNaC::diff(alpha_h, py);
        // grad_alpha(2, 0) = GiNaC::diff(alpha_h, th);
        // grad_alpha(3, 0) = GiNaC::diff(alpha_h, vx);
        // grad_alpha(4, 0) = GiNaC::diff(alpha_h, vy);
        // grad_alpha(5, 0) = GiNaC::diff(alpha_h, w);

        // // Calculate Lie derivative of alpha(h) along f
        // GiNaC::ex Lf_alpha = 0.0;
        // for (int i = 0; i < STATE_VARS; i++)
        // {
        //     Lf_alpha = Lf_alpha + grad_alpha(i, 0) * f(i, 0);
        // }

        // Correct version:
        GiNaC::ex Bc_sym = lf2h_sym + alpha(lfh_sym, gamma) + alpha(lfh_sym + alpha(h, gamma), gamma);

        Ac_conn = grad_h;
        Bc_conn = Bc_sym;

        return std::make_pair(Ac_conn, Bc_conn);
    }

    Eigen::VectorXd ConnectivityCBF::getConnConstraints(Eigen::VectorXd state, Eigen::MatrixXd robot_states, Eigen::VectorXd eigenvec)
    {
        Eigen::VectorXd Ac(CONTROL_VARS);
        Ac.setZero();
        GiNaC::matrix Ac_eval = matrixSubsFull(Ac_conn, robot_states, eigenvec, state, *this); // 1×2 numeric
        for (int i = 0; i < CONTROL_VARS; i++)
        {
            GiNaC::ex val = Ac_eval(0, i).evalf();
            Ac(i) = GiNaC::ex_to<GiNaC::numeric>(val).to_double();
        }

        logger->debug("Ac = L_g L_f h = ∇(L_f h) · g = ∇h ^ T\n{}", Ac);

        return Ac;
    }

    double ConnectivityCBF::getConnBound(Eigen::VectorXd state, Eigen::MatrixXd robot_states, Eigen::VectorXd eigenvec, double h_val)
    {
        GiNaC::ex Bc_eval = valueSubsFull(Bc_conn, robot_states, eigenvec, state, *this);
        Bc_eval = GiNaC::subs(Bc_eval, h == h_val);
        double Bc = GiNaC::ex_to<GiNaC::numeric>(Bc_eval).to_double();
        logger->debug("Bc = L_f² h + L_f(α(h)) + α(L_f h + α(h)) = {}", Bc);
        return Bc;
    }

    // Get the minimum distance constraint bound for the current state and agent
    // Parameters:
    //   state: Current state vector
    //   neighbor_state: Other agent state vector
    // Returns:
    //   Bound value for the minimum distance constraint
    double ConnectivityCBF::getSafetyBound(Eigen::VectorXd state, Eigen::VectorXd neighbor_state)
    {
        // Substitute numerical values and evaluate
        GiNaC::ex expr = valueSubs(Bc_safe, state, neighbor_state, *this);
        double Bc = GiNaC::ex_to<GiNaC::numeric>(expr).to_double();

        return Bc;
    }

    // Get the maximum velocity constraint bounds for the current state
    // Parameters:
    //   state: Current state vector
    // Returns:
    //   Vector of bound values for maximum velocity constraints
    Eigen::VectorXd ConnectivityCBF::getMaxVelBounds(Eigen::VectorXd state)
    {
        // Create dummy agent (not needed for velocity constraints)
        Eigen::VectorXd dummy_agent(6);
        dummy_agent.setZero();

        // Get bounds for each velocity component
        std::vector<GiNaC::ex> expressions(CONTROL_VARS);
        expressions[0] = valueSubs(Bc_v1_max, state, dummy_agent, *this);
        expressions[1] = valueSubs(Bc_v2_max, state, dummy_agent, *this);
        expressions[2] = valueSubs(Bc_v3_max, state, dummy_agent, *this);

        // Convert to Eigen vector
        Eigen::VectorXd Bs;
        Bs.resize(CONTROL_VARS);

        for (int i = 0; i < CONTROL_VARS; ++i)
        {
            Bs(i) = GiNaC::ex_to<GiNaC::numeric>(expressions[i]).to_double();
        }

        return Bs;
    }

    // Get the minimum velocity constraint bounds for the current state
    // Parameters:
    //   state: Current state vector
    // Returns:
    //   Vector of bound values for minimum velocity constraints
    Eigen::VectorXd ConnectivityCBF::getMinVelBounds(Eigen::VectorXd state)
    {
        // Create dummy agent (not needed for velocity constraints)
        Eigen::VectorXd dummy_agent(6);
        dummy_agent.setZero();

        // Get bounds for each velocity component
        std::vector<GiNaC::ex> expressions(CONTROL_VARS);
        expressions[0] = valueSubs(Bc_v1_min, state, dummy_agent, *this);
        expressions[1] = valueSubs(Bc_v2_min, state, dummy_agent, *this);
        expressions[2] = valueSubs(Bc_v3_min, state, dummy_agent, *this);

        // Convert to Eigen vector
        Eigen::VectorXd Bs;
        Bs.resize(CONTROL_VARS);

        for (int i = 0; i < CONTROL_VARS; ++i)
        {
            Bs(i) = GiNaC::ex_to<GiNaC::numeric>(expressions[i]).to_double();
        }

        return Bs;
    }

    // Set a custom alpha function for the CBFs
    // Parameters:
    //   newAlpha: New alpha function to use
    void ConnectivityCBF::setAlpha(std::function<GiNaC::ex(GiNaC::ex, double)> newAlpha)
    {
        alpha = newAlpha;
    }
}