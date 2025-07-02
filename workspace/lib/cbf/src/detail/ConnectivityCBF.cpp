#include <vector>
#include <cbf/detail/ConnectivityCBF.h>

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
    // Constructor for Connectivity Control Barrier Function
    // Initializes the CBF with minimum/maximum distance constraints and velocity limits
    // Parameters:
    //   min_dist: Minimum allowed distance between agents
    //   max_dist: Maximum allowed distance between agents (connectivity range)
    //   vmin: Minimum velocity limits for each control dimension
    //   vmax: Maximum velocity limits for each control dimension
    ConnectivityCBF::ConnectivityCBF(double min_dist, double max_dist, Eigen::VectorXd vmin, Eigen::VectorXd vmax)
        : dmin(min_dist), dmax(max_dist), vmin(vmin), vmax(vmax),
          px("px"), py("py"), th("th"), vx("vx"), vy("vy"), w("w"), xt("xt"), yt("yt")
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

        // Create state vector and transpose it for proper matrix operations
        state = {{px, py, th, vx, vy, w}};
        state = state.transpose();

        // Create agent state vector and transpose it
        x_agent = {{xt, yt}};
        x_agent = x_agent.transpose();

        // Define system dynamics (f + g*u)
        // f = A*x (drift term)
        // g = B (control input matrix)
        f = A.mul(state);
        g = B;

        // Set default alpha function to cubic for better performance
        alpha = myAlpha;

        // Initialize all Control Barrier Functions (CBFs)

        // Minimum distance CBF: Ensures minimum safe distance is maintained between agents
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
        std::cout << "Closing Connectivity CBF ..." << std::endl;
    }

    // Initialize the minimum distance Control Barrier Function
    // This CBF ensures that the robot maintains a minimum safe distance (dmin) from other agents
    // Returns:
    //   A pair containing the CBF constraint matrix (Ac) and bound (Bc)
    // TODO: current the CBF has no knowledge of the other agents' velocity. SHOULD IT?
    std::pair<GiNaC::matrix, GiNaC::ex> ConnectivityCBF::initSafetyCBF()
    {
        // Calculate relative position to other agent
        GiNaC::matrix d = x_agent.sub(GiNaC::matrix{{px}, {py}});

        // Calculate squared distance to other agent
        GiNaC::ex norm2 = GiNaC::pow(d(0, 0), 2) + GiNaC::pow(d(1, 0), 2);

        // Define barrier function: h(x) = ||x - x_agent||^2 - dmin^2
        // h > 0 when distance is greater than minimum distance
        GiNaC::ex b1 = norm2 - GiNaC::pow(dmin, 2);

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
        for (int i = 0; i < STATE_VARS; i++)
        {
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
        for (int i = 0; i < STATE_VARS; i++)
        {
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
        for (int i = 0; i < STATE_VARS; i++)
        {
            lfb_c = lfb_c + grad_bc(i, 0) * f(i, 0);
        }

        // Calculate L_g L_f h = ∇(L_f h) · g for each control input
        // This forms the constraint matrix Ac for the QP solver
        GiNaC::matrix Ac1 = GiNaC::matrix(1, CONTROL_VARS);
        for (int j = 0; j < CONTROL_VARS; j++)
        {
            GiNaC::ex Ac1j = 0.0;
            for (int i = 0; i < STATE_VARS; i++)
            {
                Ac1j += grad2_b1(i, 0) * g(i, j);
            }
            Ac1(0, j) = Ac1j;
        }

        // Calculate the constraint bound: L_f^2 h + L_f alpha(h) + alpha(L_f h + alpha(h))
        // This ensures the relative degree 2 CBF condition is satisfied
        GiNaC::ex B1 = lf2b1;
        GiNaC::ex B2 = lfb_c;
        GiNaC::ex psi1 = lfb1 + alpha_b;
        GiNaC::ex B3 = alpha(psi1, gamma);
        GiNaC::ex Bc1 = B1 + B2 + B3;

        return std::make_pair(Ac1, Bc1);
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

    // Substitute state and agent values into a symbolic matrix
    // Parameters:
    //   a: Symbolic matrix
    //   state: Current state vector
    //   agent_state: Other agent state vector
    // Returns:
    //   Matrix with numerical values substituted
    GiNaC::ex ConnectivityCBF::matrixSubs(GiNaC::matrix a, Eigen::VectorXd state, Eigen::VectorXd agent_state)
    {
        // Substitute each state variable with its numerical value
        GiNaC::ex tmp = GiNaC::subs(a, px == state(0));
        tmp = GiNaC::subs(tmp, py == state(1));
        tmp = GiNaC::subs(tmp, th == state(2));
        tmp = GiNaC::subs(tmp, vx == state(3));
        tmp = GiNaC::subs(tmp, vy == state(4));
        tmp = GiNaC::subs(tmp, w == state(5));
        tmp = GiNaC::subs(tmp, xt == agent_state(0));
        tmp = GiNaC::subs(tmp, yt == agent_state(1));
        return tmp;
    }

    // Substitute state and agent values into a symbolic expression
    // Parameters:
    //   a: Symbolic expression
    //   state: Current state vector
    //   agent_state: Other agent state vector
    // Returns:
    //   Expression with numerical values substituted
    GiNaC::ex ConnectivityCBF::valueSubs(GiNaC::ex a, Eigen::VectorXd state, Eigen::VectorXd agent_state)
    {
        // Substitute each state variable with its numerical value
        GiNaC::ex tmp = GiNaC::subs(a, px == state(0));
        tmp = GiNaC::subs(tmp, py == state(1));
        tmp = GiNaC::subs(tmp, th == state(2));
        tmp = GiNaC::subs(tmp, vx == state(3));
        tmp = GiNaC::subs(tmp, vy == state(4));
        tmp = GiNaC::subs(tmp, w == state(5));
        tmp = GiNaC::subs(tmp, xt == agent_state(0));
        tmp = GiNaC::subs(tmp, yt == agent_state(1));
        return tmp;
    }

    // Get the minimum distance constraint vector for the current state and agent
    // Parameters:
    //   state: Current state vector
    //   agent_state: Other agent state vector
    // Returns:
    //   Constraint vector for the QP solver
    Eigen::VectorXd ConnectivityCBF::getSafetyConstraints(Eigen::VectorXd state, Eigen::VectorXd agent_state)
    {
        // Substitute numerical values into symbolic matrix and convert to Eigen vector
        GiNaC::ex matrix_expr = matrixSubs(Ac_safe, state, agent_state);
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
        Eigen::Vector2d dummy_agent;
        dummy_agent.setZero();

        // Get constraint matrices for each velocity component
        std::vector<GiNaC::ex> expressions(CONTROL_VARS);
        expressions[0] = matrixSubs(Ac_v1_max, state, dummy_agent);
        expressions[1] = matrixSubs(Ac_v2_max, state, dummy_agent);
        expressions[2] = matrixSubs(Ac_v3_max, state, dummy_agent);

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
        Eigen::Vector2d dummy_agent;
        dummy_agent.setZero();

        // Get constraint matrices for each velocity component
        std::vector<GiNaC::ex> expressions(CONTROL_VARS);
        expressions[0] = matrixSubs(Ac_v1_min, state, dummy_agent);
        expressions[1] = matrixSubs(Ac_v2_min, state, dummy_agent);
        expressions[2] = matrixSubs(Ac_v3_min, state, dummy_agent);

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

    // Given the numerical robot positions, compute the second smallest eigenvalue (lambda2) for the Laplacian matrix
    // Parameters:
    //   robot_positions: Matrix of robot positions
    //   Rs_value: Maximum distance for connectivity
    //   sigma_value: Parameter for the weight function
    // Returns:
    //   A pair containing the second smallest eigenvalue (lambda2) and the corresponding eigenvector
    std::pair<double, Eigen::VectorXd> getLambda2FromL(
        const Eigen::MatrixXd &robot_positions,
        double Rs_value,
        double sigma_value)
    {
        const int N = robot_positions.rows();
        // Construct the Adjacency matrix A in numerical form
        Eigen::MatrixXd A_num = Eigen::MatrixXd::Zero(N, N);
        double Rs_value2 = Rs_value * Rs_value;
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
        SPDLOG_INFO("Diagonal matrix D:\n{}", D_num);
        SPDLOG_INFO("Adjacency matrix A:\n{}", A_num);
        SPDLOG_INFO("Laplacian matrix L:\n{}", L_num);

        // Numerically solve for the second smallest eigenvalue and corresponding eigenvector
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(L_num);
        Eigen::VectorXd eigenvals = solver.eigenvalues();
        Eigen::MatrixXd eigenvecs = solver.eigenvectors();
        return {eigenvals(1), eigenvecs.col(1)};
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

    // Substitute symbolic matrix with numerical values for robot positions and eigenvector
    // Parameters:
    //   expr_matrix: Symbolic matrix to substitute
    //   robot_positions: Matrix of robot positions (N x 2)
    //   eigenvec: Eigenvector containing velocity values for each robot
    //   self_position: Position of the current robot (optional, defaults to zero vector)
    GiNaC::matrix ConnectivityCBF::matrixSubsMatrix(
        const GiNaC::matrix &expr_matrix,
        const Eigen::MatrixXd &robot_positions,
        const Eigen::VectorXd &eigenvec,
        const Eigen::Vector2d &self_position)
    {
        GiNaC::exmap substitutions;
        // 机器人位置和特征值向量替换
        for (int i = 0; i < robot_positions.rows(); ++i)
        {
            substitutions[px_list[i]] = robot_positions(i, 0);
            substitutions[py_list[i]] = robot_positions(i, 1);
            substitutions[eigenvec_list[i]] = eigenvec(i);
            std::cout << "Substituting px" << i << ": " << robot_positions(i, 0) << ", py" << i << ": " << robot_positions(i, 1) << ", eigenvec" << i << ": " << eigenvec(i) << std::endl;
        }
        // 当前机器人的自身位置
        substitutions[px] = self_position(0);
        substitutions[py] = self_position(1);
        // 对整个矩阵统一替换
        GiNaC::matrix result = GiNaC::ex_to<GiNaC::matrix>(expr_matrix.subs(substitutions));
        return result;
    }

    // Symbolically compute the gradient of the barrier function with respect to state variables
    // Parameters:
    //   N: Number of agents
    //   Rs: Symbolic representation of the maximum distance for connectivity
    //   sigma: Symbolic representation of the weight function parameter
    // TODO: this could be optimized: the A matrices are already calulated earlier
    GiNaC::matrix ConnectivityCBF::compute_dh_dx(int N, const GiNaC::ex &Rs, const GiNaC::ex &sigma)
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
                // print eigenvec_list[i], eigenvec_list[j];
                GiNaC::ex vdiff2 = GiNaC::pow(eigenvec_list[i] - eigenvec_list[j], 2);
                dLdx += dAij_dx * vdiff2;
                dLdy += dAij_dy * vdiff2;
            }
            grad(i, 0) = dLdx;
            grad(i, 1) = dLdy;
        }
        return grad;
    }

    // Compute the second derivative (Hessian) of the barrier function with respect to state variables
    // Parameters:
    //   dh_dx_sym: Symbolic gradient of the barrier function
    //   self_idx: Index of the current robot in the gradient matrix
    // Returns:
    //   Hessian matrix (2x2) of the barrier function with respect to px and py in symbolic form
    GiNaC::matrix ConnectivityCBF::compute_d2h_dx2(const GiNaC::matrix &dh_dx_sym, int self_idx)
    {
        GiNaC::matrix hess(2, 2);
        hess(0, 0) = GiNaC::diff(dh_dx_sym(self_idx, 0), px_list[self_idx]);
        hess(0, 1) = GiNaC::diff(dh_dx_sym(self_idx, 0), py_list[self_idx]);
        hess(1, 0) = GiNaC::diff(dh_dx_sym(self_idx, 1), px_list[self_idx]);
        hess(1, 1) = GiNaC::diff(dh_dx_sym(self_idx, 1), py_list[self_idx]);
        return hess;
    }

    //
    Eigen::VectorXd ConnectivityCBF::compute_dLf_h_dx(
        const GiNaC::matrix &dh_dx_sym,
        int self_idx,
        const Eigen::MatrixXd &robot_positions,
        const Eigen::VectorXd &eigenvec,
        const Eigen::VectorXd &x_self)
    {
        // === Step 1: 构造符号 Hessian 矩阵 ===
        GiNaC::matrix hess_sym = compute_d2h_dx2(dh_dx_sym, self_idx);
        // 输出符号 Hessian 表达式（调试用）
        GiNaC::matrix hess_eval = matrixSubsMatrix(hess_sym, robot_positions, eigenvec, x_self.head<2>());
        SPDLOG_INFO("Step 1: ∇²h evaluated:\n{}", hess_eval);
        // === Step 2: 计算 Hessian 项：∇²h · f(x) ===
        double fx = x_self(3);
        double fy = x_self(4);
        Eigen::Vector2d hess_term;
        hess_term(0) = GiNaC::ex_to<GiNaC::numeric>(hess_eval(0, 0) * fx + hess_eval(0, 1) * fy).to_double();
        hess_term(1) = GiNaC::ex_to<GiNaC::numeric>(hess_eval(1, 0) * fx + hess_eval(1, 1) * fy).to_double();
        SPDLOG_INFO("Step 2: Hessian contribution ∇²h·f");
        // logVector("hess_term", hess_term);

        // === Step 3: ∇h 数值化，替换变量获得数值梯度 ===
        GiNaC::matrix dh_eval = matrixSubsMatrix(dh_dx_sym, robot_positions, eigenvec, x_self.head<2>());
        Eigen::VectorXd dh_dx = Eigen::VectorXd::Zero(6);
        // 这里只对 px 和 py 非零，其他导数为零
        dh_dx(0) = GiNaC::ex_to<GiNaC::numeric>(dh_eval(self_idx, 0)).to_double();
        dh_dx(1) = GiNaC::ex_to<GiNaC::numeric>(dh_eval(self_idx, 1)).to_double();
        SPDLOG_INFO("Step 3: ∇h(px, py)\n{}", dh_dx);

        // === Step 4: 构造 Jf^T ∇h 项 ===
        // 由于 f(x) = [vx, vy, w]，其雅可比矩阵 Jf 关于状态向量 x 的非零偏导为：
        // ∂vx/∂x4=1, ∂vy/∂x5=1, ∂w/∂x6=1，对应 transpose 后影响的就是 dh_dx 中 vx, vy 项
        Eigen::VectorXd jac_term = Eigen::VectorXd::Zero(6);
        jac_term(3) = dh_dx(0); // vx 对应 px
        jac_term(4) = dh_dx(1); // vy 对应 py
        SPDLOG_INFO("Step 4: Jf^T ∇h\n{}", jac_term);

        // === Step 5: 拼接最终结果 ∇(L_f h) = ∇²h·f + Jfᵀ∇h ===
        Eigen::VectorXd total = jac_term;
        total.head<2>() += hess_term;
        return total;
    }

    // TODO: ideally this should return a symbolic expression
    // return std::pair<GiNaC::matrix, GiNaC::ex>
    // Ac should be a vector of shape 1 x CONTROL_VARS, Bc should be a scalar
    std::pair<Eigen::VectorXd, double> ConnectivityCBF::initConnCBF(
        const Eigen::MatrixXd &robot_states, // N x 6 matrix
        const Eigen::VectorXd &x_self,       // 当前机器人的状态 6x1
        int self_idx)                        // 当前机器人在 robot_states 中的索引
    {
        const int N = robot_states.rows(); // Number of robots
        SPDLOG_INFO("Initializing Connectivity CBF with {} robots", N);

        // Step 1: h = λ₂ - λ₂_min (numerically)
        const double sigma_val = dmax * dmax * dmax * dmax / std::log(2.0); // σ = R_s^4 / ln(2)
        auto [lambda2_val, eigenvec] = getLambda2FromL(robot_states.leftCols(2), dmax, sigma_val);
        eigenvec = eigenvec / eigenvec.norm();
        double h = lambda2_val - epsilon;
        SPDLOG_INFO("Step 1: λ₂ = {}, λ₂_min = {}, h = {}", lambda2_val, epsilon, h);

        // Step 2: 符号构造 ∇h (gradient of h), shape = N×2
        initSymbolLists(N);
        GiNaC::matrix dh_dx_sym = compute_dh_dx(N, dmax, sigma_val);                                 // shape Nx2
        GiNaC::matrix dh_dx_ginac = matrixSubsMatrix(dh_dx_sym, robot_states.leftCols(2), eigenvec); // N×2 数值表达式, robot_states.leftCols(2) 只取 px, py
        SPDLOG_INFO("Step 2: ∇h evaluated\n{}", dh_dx_ginac);
        Eigen::VectorXd dh_dx = Eigen::VectorXd::Zero(STATE_VARS);
        dh_dx(0) = GiNaC::ex_to<GiNaC::numeric>(dh_dx_ginac(self_idx, 0)).to_double();
        dh_dx(1) = GiNaC::ex_to<GiNaC::numeric>(dh_dx_ginac(self_idx, 1)).to_double();

        // Step 3: L_f h = ∇h · f(x_self)
        // 由于 f = A*x = [vx, vy, w, 0, 0, 0]，直接从 x_self 构造 f 的数值向量 // TODO: this could potentially be replaced by f in fields (if using symbolic)
        Eigen::VectorXd f_x = Eigen::VectorXd::Zero(STATE_VARS);
        f_x(0) = x_self(3); // vx
        f_x(1) = x_self(4); // vy
        f_x(2) = x_self(5); // w
        double lfh = dh_dx.dot(f_x);
        SPDLOG_INFO("Step 3: L_f h = ∇h · f = {}", lfh);

        // Step 4: 用符号 Hessian 构造 ∇(L_f h)
        Eigen::VectorXd dlfh_dx = compute_dLf_h_dx(
            dh_dx_sym,
            self_idx,
            robot_states.leftCols(2),
            eigenvec,
            x_self);
        SPDLOG_INFO("Step 4: ∇(L_f h)\n{}", dlfh_dx);

        //  Step 5: L_f² h = ∇(L_f h) · f
        double lf2h = dlfh_dx.dot(f_x);
        SPDLOG_INFO("Step 5: L_f² h = ∇(L_f h) · f = {}", lf2h);

        //  Step 6: L_g L_f h = ∇(L_f h) · g
        Eigen::MatrixXd g = Eigen::MatrixXd::Zero(STATE_VARS, CONTROL_VARS); // shape 6x3 // TODO: this could potentially be replaced by g in fields (if using symbolic)
        g(3, 0) = 1.0;
        g(4, 1) = 1.0;
        Eigen::VectorXd Ac = (dlfh_dx.transpose() * g).transpose(); // size = 3
        Ac(2) = 0.0;                                                // 防止角速度影响控制约束 // TODO: why is this necessary
        SPDLOG_INFO("Step 6: Ac = L_g L_f h = ∇(L_f h) · g\n{}", Ac);

        // Step 7: Bc = L_f² h + L_f(α(h)) + α(L_f h + α(h)), where α is a class of K functions
        double psi1 = lfh + gamma * h; // psi1 = L_f h + alpha1(h), assuming linear alpha1
        double Bc = lf2h               // L_f^2 h
                    + gamma * lfh      // Assuming linear alpha11: L_f(α(h)) = γ L_f h
                    + gamma * psi1;    // Assuming linear alpha2:
        SPDLOG_INFO("Step 7: Bc = L_f² h + L_f(α(h)) + α(L_f h + α(h)) = {} + {} + {} = {}", lf2h, gamma * lfh, gamma * psi1, Bc);
        return std::make_pair(Ac, Bc);

        //         // ✅ Step 7: 使用一致方式计算 Bc
        // // alpha(h)
        // GiNaC::ex h_sym = h;  // 直接用 double 值代入
        // GiNaC::ex alpha_h_expr = alpha(h_sym, gamma);
        // double alpha_h = GiNaC::ex_to<GiNaC::numeric>(alpha_h_expr).to_double();

        // // L_f α(h) = α'(h) * L_f h, 以 fifthAlpha = γ h^5 为例，其导数为 5γ h^4
        // double d_alpha_h = 5.0 * gamma * std::pow(h, 4);
        // double lf_alpha_h = d_alpha_h * lfh;

        // // α(L_f h + α(h))
        // GiNaC::ex nested_expr = alpha(lfh + alpha_h, gamma);
        // double nested_alpha = GiNaC::ex_to<GiNaC::numeric>(nested_expr).to_double();

        // // Final Bc
        // double Bc = lf2h + lf_alpha_h + nested_alpha;
        // std::cout << "[CBF] Bc = L_f² h + L_f(α(h)) + α(L_f h + α(h)) = " << Bc << std::endl;

        // return std::make_pair(Ac, Bc);
    }

    // TODO: Deprecated: check initConnCBF
    // Eigen::VectorXd ConnectivityCBF::getConnConstraints(
    //     const Eigen::VectorXd &x_self,
    //     const std::vector<Eigen::VectorXd> &other_positions)
    // {
    // }
    // double ConnectivityCBF::getConnBound(
    //     const Eigen::VectorXd &x_self,
    //     const std::vector<Eigen::VectorXd> &other_positions)
    // {
    // }

    // Get the minimum distance constraint bound for the current state and agent
    // Parameters:
    //   state: Current state vector
    //   agent_state: Other agent state vector
    // Returns:
    //   Bound value for the minimum distance constraint
    double ConnectivityCBF::getSafetyBound(Eigen::VectorXd state, Eigen::VectorXd agent_state)
    {
        // Substitute numerical values and evaluate
        GiNaC::ex expr = valueSubs(Bc_safe, state, agent_state);
        double Bc = GiNaC::ex_to<GiNaC::numeric>(expr).to_double();

        return Bc;
    }

    // Get the maximum distance constraint bound for the current state and agent
    // Parameters:
    //   state: Current state vector
    //   agent_state: Other agent state vector
    // Returns:
    //   Bound value for the maximum distance constraint
    double ConnectivityCBF::getMaxDistBound(Eigen::VectorXd state, Eigen::VectorXd agent_state)
    {
        // Substitute numerical values and evaluate
        GiNaC::ex expr = valueSubs(Bc_connectivity, state, agent_state);
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
        Eigen::Vector2d dummy_agent;
        dummy_agent.setZero();

        // Get bounds for each velocity component
        std::vector<GiNaC::ex> expressions(CONTROL_VARS);
        expressions[0] = valueSubs(Bc_v1_max, state, dummy_agent);
        expressions[1] = valueSubs(Bc_v2_max, state, dummy_agent);
        expressions[2] = valueSubs(Bc_v3_max, state, dummy_agent);

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
        Eigen::Vector2d dummy_agent;
        dummy_agent.setZero();

        // Get bounds for each velocity component
        std::vector<GiNaC::ex> expressions(CONTROL_VARS);
        expressions[0] = valueSubs(Bc_v1_min, state, dummy_agent);
        expressions[1] = valueSubs(Bc_v2_min, state, dummy_agent);
        expressions[2] = valueSubs(Bc_v3_min, state, dummy_agent);

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