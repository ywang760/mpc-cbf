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

        // Maximum distance CBF: Ensures agents stay within communication/connectivity range
        auto res2 = initConnectivityCBF();
        Ac_connectivity = res2.first;
        Bc_connectivity = res2.second;

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

    // Initialize the maximum distance Control Barrier Function
    // This CBF ensures that the robot maintains a maximum distance (dmax) to preserve connectivity
    // Returns:
    //   A pair containing the CBF constraint matrix (Ac) and bound (Bc)
    std::pair<GiNaC::matrix, GiNaC::ex> ConnectivityCBF::initConnectivityCBF()
    {
        // Calculate relative position to other agent
        GiNaC::matrix d = x_agent.sub(GiNaC::matrix{{px}, {py}});

        // Calculate squared distance to other agent
        GiNaC::ex norm2 = GiNaC::pow(d(0, 0), 2) + GiNaC::pow(d(1, 0), 2);

        // Define barrier function: h(x) = -||x - x_agent||^2 + dmax^2
        // h > 0 when distance is less than maximum distance
        GiNaC::ex b2 = -norm2 + GiNaC::pow(dmax, 2);

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
        for (int i = 0; i < STATE_VARS; i++)
        {
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
        for (int i = 0; i < STATE_VARS; i++)
        {
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
        for (int i = 0; i < STATE_VARS; i++)
        {
            lfb_c = lfb_c + grad_bc(i, 0) * f(i, 0);
        }

        // Calculate L_g L_f h = ∇(L_f h) · g for each control input
        GiNaC::matrix Ac2 = GiNaC::matrix(1, CONTROL_VARS);
        for (int j = 0; j < CONTROL_VARS; j++)
        {
            GiNaC::ex Ac2j = 0.0;
            for (int i = 0; i < STATE_VARS; i++)
            {
                Ac2j += grad2_b2(i, 0) * g(i, j);
            }
            Ac2(0, j) = Ac2j;
        }

        // Calculate the constraint bound for the maximum distance boundary
        GiNaC::ex B1 = lf2b2;
        GiNaC::ex B2 = lfb_c;
        GiNaC::ex B3 = alpha(lfb2 + alpha_b2, gamma);
        GiNaC::ex Bc2 = B1 + B2 + B3;

        return std::make_pair(Ac2, Bc2);
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

    Eigen::VectorXd ConnectivityCBF::getConnectivityConstraints(Eigen::VectorXd state, Eigen::VectorXd agent_state)
    {
        // TODO: Implement this function
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