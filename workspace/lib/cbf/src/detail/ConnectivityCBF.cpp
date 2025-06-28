#include <vector>
#include <cbf/detail/ConnectivityCBF.h>
#include <fstream> 
#include <tuple>

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
        gamma = 0.1; // Convergence rate parameter

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

        // Set default alpha function to fifth-order for better performance
        alpha = fifthAlpha;

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
        GiNaC::ex B3 = alpha(lfb1 + alpha_b, gamma);
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
        //逐步替换，比如第一行先把a里的px替换为state(0)，然后后面是基于这一步替换，再做把tmp中的py替换为state(1)...
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
// 将已经构造好的符号表达式（即带有 px, py 等符号变量的 GiNaC::matrix）中，
// 替换变量为具体数值，从而得到可以求值的表达式或数值矩阵。
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

    Eigen::RowVectorXd cbf::ConnectivityCBF::getConnConstraints(const Eigen::VectorXd& x_self,
                                      const std::vector<Eigen::Vector3d>& other_positions,
                                      double& h_out)  // 新增引用参数
    {
        const int N = 1 + other_positions.size();
        Eigen::MatrixXd robot_states(N, 6);
        robot_states.setZero();  // 初始化为 0（可保留）

        // 当前机器人位置和速度完整设置
        robot_states(0, 0) = x_self(0);  // px
        robot_states(0, 1) = x_self(1);  // py
        robot_states(0, 2) = x_self(2);  // yaw
        robot_states(0, 3) = x_self(3);  // vx
        robot_states(0, 4) = x_self(4);  // vy
        robot_states(0, 5) = x_self(5);  // yaw rate

      for (int i = 0; i < other_positions.size(); ++i) {
            robot_states(i + 1, 0) = other_positions[i][0];  // or .x()
            robot_states(i + 1, 1) = other_positions[i][1];  // or .y()
        }

        const double Rs = 3.0;
        const double sigma = 120;
        const double lambda2_min = 0.1;
        const double gamma = 1.0;
        auto [Ac, Bc, h] = this->initConnCBF(robot_states, x_self, 0, Rs, sigma, lambda2_min, gamma);
        h_out = h;
        return Ac;
    }


    std::pair<double, Eigen::VectorXd> getLambda2FromL(
    const Eigen::MatrixXd& robot_positions, 
    double Rs_value, 
    double sigma_value)
    {
        const int N = robot_positions.rows();
        Eigen::MatrixXd A_num = Eigen::MatrixXd::Zero(N, N);
        for (int i = 0; i < N; ++i)
        {
            for (int j = 0; j < N; ++j)
            {
                if (i != j)
                {
                    double dx = robot_positions(i, 0) - robot_positions(j, 0);
                    double dy = robot_positions(i, 1) - robot_positions(j, 1);
                    double dij2 = dx * dx + dy * dy;

                    if (dij2 <= Rs_value * Rs_value)
                    {
                        double weight = std::exp(std::pow(Rs_value * Rs_value - dij2, 2) / sigma_value) - 1;
                        A_num(i, j) = weight;
                    }
                }
            }
        }
        Eigen::MatrixXd D_num = Eigen::MatrixXd::Zero(N, N);
        for (int i = 0; i < N; ++i)
            D_num(i, i) = A_num.row(i).sum();
        Eigen::MatrixXd L_num = D_num - A_num;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(L_num);
        Eigen::VectorXd eigenvals = solver.eigenvalues();
        Eigen::MatrixXd eigenvecs = solver.eigenvectors();
        return {eigenvals(1), eigenvecs.col(1)};
    }

    void ConnectivityCBF::initSymbolLists(int N) {
        px_list.clear(); py_list.clear(); v_list.clear();
        for (int i = 0; i < N; ++i) {
            px_list.emplace_back(GiNaC::symbol("px" + std::to_string(i)));
            py_list.emplace_back(GiNaC::symbol("py" + std::to_string(i)));
            v_list.emplace_back(GiNaC::symbol("v"  + std::to_string(i)));
        }
    }

    GiNaC::matrix ConnectivityCBF::matrixSubsMatrix(
        const GiNaC::matrix& expr_matrix,
        const Eigen::MatrixXd& robot_positions,
        const Eigen::VectorXd& eigenvec,
        const Eigen::Vector2d& self_position)
    {
        GiNaC::exmap substitutions;
        // 机器人位置和特征值向量替换
        for (int i = 0; i < robot_positions.rows(); ++i) {
            substitutions[px_list[i]] = robot_positions(i, 0);
            substitutions[py_list[i]] = robot_positions(i, 1);
            substitutions[v_list[i]]  = eigenvec(i);
        }
        // 当前机器人的自身位置
        substitutions[px] = self_position(0);
        substitutions[py] = self_position(1);
        // 对整个矩阵统一替换
        GiNaC::matrix result = GiNaC::ex_to<GiNaC::matrix>(expr_matrix.subs(substitutions));
        return result;
    }


    GiNaC::matrix ConnectivityCBF::compute_dh_dx(int N, const GiNaC::ex& Rs, const GiNaC::ex& sigma)
    {
        initSymbolLists(N);
        GiNaC::matrix grad(N, 2);
        for (int i = 0; i < N; ++i) {
            GiNaC::ex dLdx = 0, dLdy = 0;
            for (int j = 0; j < N; ++j) {
                if (i == j) continue;
                GiNaC::ex dx = px_list[i] - px_list[j];
                GiNaC::ex dy = py_list[i] - py_list[j];
                GiNaC::ex dij2 = dx*dx + dy*dy;
                GiNaC::ex diff = GiNaC::pow(Rs, 2) - dij2;
                GiNaC::ex Aij = GiNaC::exp(GiNaC::pow(diff, 2) / sigma);
                GiNaC::ex dAij_dx = -4 * Aij * diff / sigma * dx;
                GiNaC::ex dAij_dy = -4 * Aij * diff / sigma * dy;
                GiNaC::ex vdiff2 = GiNaC::pow(v_list[i] - v_list[j], 2);
                dLdx += dAij_dx * vdiff2;
                dLdy += dAij_dy * vdiff2;
            }
            grad(i, 0) = dLdx;
            grad(i, 1) = dLdy;
        }
        return grad;
    }



    GiNaC::matrix ConnectivityCBF::compute_d2h_dx2(const GiNaC::matrix& dh_dx_sym, int self_idx)
    {
        GiNaC::matrix hess(2, 2);
        hess(0, 0) = GiNaC::diff(dh_dx_sym(self_idx, 0), px_list[self_idx]);
        hess(0, 1) = GiNaC::diff(dh_dx_sym(self_idx, 0), py_list[self_idx]);
        hess(1, 0) = GiNaC::diff(dh_dx_sym(self_idx, 1), px_list[self_idx]);
        hess(1, 1) = GiNaC::diff(dh_dx_sym(self_idx, 1), py_list[self_idx]);
        return hess;
    }

    Eigen::VectorXd ConnectivityCBF::compute_dLf_h_dx(
        const GiNaC::matrix& dh_dx_sym,
        int self_idx,
        const Eigen::MatrixXd& robot_positions,
        const Eigen::VectorXd& eigenvec,
        const Eigen::VectorXd& x_self)
    {
        // === Step 1: 构造符号 Hessian 矩阵 ===
        GiNaC::matrix hess_sym = compute_d2h_dx2(dh_dx_sym, self_idx);
        // 输出符号 Hessian 表达式（调试用）
        GiNaC::matrix hess_eval = matrixSubsMatrix(hess_sym, robot_positions, eigenvec, x_self.head<2>());
        // 输出数值化后的 Hessian 矩阵
        //std::cout << "[CBF] Evaluated Hessian matrix (∇²h evaluated):" << std::endl;
        // for (int i = 0; i < 2; ++i) {
        //     for (int j = 0; j < 2; ++j) {
        //         std::cout << "  H_eval(" << i << "," << j << ") = "
        //                 << GiNaC::ex_to<GiNaC::numeric>(hess_eval(i, j)).to_double() << std::endl;
        //     }
        // }
        // === Step 2: 计算 Hessian 项：∇²h · f(x) ===
        double fx = x_self(3);
        double fy = x_self(4);
        Eigen::Vector2d hess_term;
        hess_term(0) = GiNaC::ex_to<GiNaC::numeric>(hess_eval(0, 0) * fx + hess_eval(0, 1) * fy).to_double();
        hess_term(1) = GiNaC::ex_to<GiNaC::numeric>(hess_eval(1, 0) * fx + hess_eval(1, 1) * fy).to_double();
        //std::cout << "[CBF] Step 2: Hessian contribution ∇²h·f = (" << hess_term(0) << ", " << hess_term(1) << ")" << std::endl;

        // === Step 3: ∇h 数值化，替换变量获得数值梯度 ===
        GiNaC::matrix dh_eval = matrixSubsMatrix(dh_dx_sym, robot_positions, eigenvec, x_self.head<2>());
        Eigen::VectorXd dh_dx = Eigen::VectorXd::Zero(6);
        // 这里只对 px 和 py 非零，其他导数为零
        dh_dx(0) = GiNaC::ex_to<GiNaC::numeric>(dh_eval(self_idx, 0)).to_double();
        dh_dx(1) = GiNaC::ex_to<GiNaC::numeric>(dh_eval(self_idx, 1)).to_double();
        //std::cout << "[CBF] Step 3: ∇h(px, py) = (" << dh_dx(0) << ", " << dh_dx(1) << ")" << std::endl;

         // === Step 4: 构造 Jf^T ∇h 项 ===
        // 由于 f(x) = [vx, vy, w]，其雅可比矩阵 Jf 关于状态向量 x 的非零偏导为：
        // ∂vx/∂x4=1, ∂vy/∂x5=1, ∂w/∂x6=1，对应 transpose 后影响的就是 dh_dx 中 vx, vy 项
        Eigen::VectorXd jac_term = Eigen::VectorXd::Zero(6);
        jac_term(0) = 0.0;
        jac_term(1) = 0.0;
        jac_term(2) = 0.0;
        jac_term(3) = dh_dx(0);  // vx 对应 px
        jac_term(4) = dh_dx(1);  // vy 对应 py
        jac_term(5) = 0.0;
        // std::cout << "[CBF] Step 4: Jf^T ∇h = (";
        // for (int i = 0; i < 6; ++i) std::cout << jac_term(i) << (i < 5 ? ", " : "");
        // std::cout << ")" << std::endl;

        // === Step 5: 拼接最终结果 ∇(L_f h) = ∇²h·f + Jfᵀ∇h ===
        Eigen::VectorXd total = Eigen::VectorXd::Zero(6);
        total(0) = hess_term(0);
        total(1) = hess_term(1);
        total.segment<6>(0) += jac_term;
        return total;
    }


    std::tuple<Eigen::RowVectorXd, double, double> ConnectivityCBF::initConnCBF(
        const Eigen::MatrixXd& robot_states,  // N x 6 matrix
        const Eigen::VectorXd& x_self,        // 当前机器人的状态 6x1
        int self_idx,                         // 当前机器人在 robot_states 中的索引
        double Rs_val,
        double sigma_val,
        double lambda2_min,
        double gamma)
    {
        const int N = robot_states.rows();
        const int STATE_VARS = 6;
        const int CONTROL_VARS = 3;
        // Step 1: h = λ₂ - λ₂_min
        auto [lambda2_val, eigenvec] = getLambda2FromL(robot_states.leftCols(2), Rs_val, sigma_val);
        //std::cout << "[CBF] λ₂ = " << lambda2_val << std::endl;
        eigenvec = eigenvec / eigenvec.norm();
        double h = lambda2_val - lambda2_min;
        std::cout << "1[CBF] h = λ₂ - λ₂_min = " << h << std::endl;
        // Step 2: 符号构造 ∇h
        initSymbolLists(N);
        GiNaC::matrix dh_dx_sym = compute_dh_dx(N, Rs_val, sigma_val);
        GiNaC::matrix dh_dx_ginac = matrixSubsMatrix(dh_dx_sym, robot_states.leftCols(2), eigenvec); // N×2 数值表达式
        //std::cout << "[CBF] ∇h evaluated:" << std::endl;
        // for (int i = 0; i < N; ++i)
        //     std::cout << "  grad[" << i << "] = (" 
        //             << GiNaC::ex_to<GiNaC::numeric>(dh_dx_ginac(i, 0)).to_double() << ", "
        //             << GiNaC::ex_to<GiNaC::numeric>(dh_dx_ginac(i, 1)).to_double() << ")" << std::endl;
        //转换为 Eigen 数值矩阵
        Eigen::MatrixXd dh_dx_all(N, 2);
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < 2; ++j)
                dh_dx_all(i, j) = GiNaC::ex_to<GiNaC::numeric>(dh_dx_ginac(i, j)).to_double();
        // 拼接 6 维梯度向量（只对 px, py 有非零导数）
        Eigen::VectorXd dh_dx = Eigen::VectorXd::Zero(6);
        dh_dx(0) = dh_dx_all(self_idx, 0);  // ∂h/∂px
        dh_dx(1) = dh_dx_all(self_idx, 1);  // ∂h/∂py
        // Step 3: L_f h = ∇h · f(x_self)
        // 由于 f = A*x = [vx, vy, w, 0, 0, 0]，直接从 x_self 构造 f 的数值向量
        Eigen::VectorXd f_x = Eigen::VectorXd::Zero(6);
        f_x(0) = x_self(3); // vx
        f_x(1) = x_self(4); // vy
        f_x(2) = x_self(5); // w
        double lfh = dh_dx.dot(f_x);
        //std::cout << "[CBF] L_f h = ∇h · f = " << lfh << std::endl;
        // ✅ Step 4: 用符号 Hessian 构造 ∇(L_f h)
        Eigen::VectorXd dlfh_dx = compute_dLf_h_dx(
            dh_dx_sym,
            self_idx,
            robot_states.leftCols(2),
            eigenvec,
            x_self
        );
        //std::cout << "[CBF] ∇(L_f h): dlfh_dx = ";
       //std::cout << dlfh_dx.transpose() << std::endl;
        // Step 5: L_f² h = ∇(L_f h) · f
        double lf2h = dlfh_dx.dot(f_x);
        //std::cout << "[CBF] L_f² h = ∇(L_f h) · f = " << lf2h << std::endl;
        // Step 6: L_g L_f h = ∇(L_f h) · g
        Eigen::MatrixXd g = Eigen::MatrixXd::Zero(STATE_VARS, CONTROL_VARS);
        g(3, 0) = 1.0;
        g(4, 1) = 1.0;
        Eigen::RowVectorXd Ac = dlfh_dx.transpose() * g; // size = 2
        Ac(2) = 0.0; // 防止角速度影响控制约束
        // Step 7: Bc = L_f² h + γ L_f h + γ² h
        double gamma1 = 1.0;
        double gamma2 = 1.0;
        double psi1 = lfh + gamma1 * h; // psi1 = L_f h + alpha1(h), assuming linear alpha1
        double Bc = lf2h                // L_f^2 h
                    + gamma1 * lfh      // Assuming linear alpha11: L_f(α(h)) = γ L_f h
                    + gamma2 * psi1;    // Assuming linear alpha2:
        //std::cout << "[CBF] Ac (∇(L_f h)·g): " << Ac << std::endl;
        //std::cout << "[CBF] Bc = L_f² h + γ L_f h + γ² h = " << Bc << std::endl;
        return std::make_tuple(Ac, Bc, h);

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


    double ConnectivityCBF::getConnBound(
        const Eigen::VectorXd& x_self,
        const std::vector<Vector3d>& other_positions)
    {
        // === 组装 robot_states: [self; others] ===
        const int N = 1 + other_positions.size();
        Eigen::MatrixXd robot_states(N, 6);
        robot_states.setZero();
        // 当前机器人：位置和速度（只填 vx, vy）
        robot_states.row(0).head<2>() = x_self.head<2>();         // px, py
        robot_states.row(0).segment<2>(3) = x_self.segment<2>(3); // vx, vy
        // 其他机器人：只填 px, py，速度设为 0
        for (int i = 0; i < other_positions.size(); ++i) {
            robot_states.row(i + 1).head<2>() = other_positions[i].head<2>();
        }
        // === 参数设置 ===
        const double Rs = 3.0;
        const double sigma = 120;
        const double lambda2_min = 0.1;
        const double gamma = 1.0;
        // === 计算约束 ===
        auto [Ac, Bc, h] = initConnCBF(robot_states, x_self, 0, Rs, sigma, lambda2_min, gamma);
        return Bc;
    }

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