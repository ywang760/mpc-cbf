class LaplacianTest : public ::testing::Test {
    protected:
        void SetUp() override {
        }
    
        void TearDown() override {
        }
    };

    TEST_F(LaplacianTest, SymbolicLaplacianEvaluateTriangle) {
        // 参数初始化
        N = 3;
        double Rs_val = 1.5;
        double sigma_val = 0.5;
        Rs = Rs_val;
        sigma = sigma_val;
    
        // 构造符号变量 px, py
        px.clear(); py.clear();
        for (int i = 0; i < N; ++i) {
            px.emplace_back(GiNaC::symbol("px" + std::to_string(i)));
            py.emplace_back(GiNaC::symbol("py" + std::to_string(i)));
        }
    
        // 构造符号邻接矩阵 A、度矩阵 D、拉普拉斯矩阵 L
        auto [A, D, L] = computeLaplacian();
    
        // 正三角形节点坐标
        std::vector<std::pair<double, double>> positions = {
            {0.0, 0.0},
            {1.0, 0.0},
            {0.5, 0.866}  // 高度为 sqrt(3)/2 ≈ 0.866
        };
    
        // 构造 GiNaC::lst 替代符号变量
        GiNaC::lst subs;
        for (int i = 0; i < N; ++i) {
            subs.append(px[i] == positions[i].first);
            subs.append(py[i] == positions[i].second);
        }
    
        // 手动计算理论 A(i,j) 值（所有边长 = 1）
        double dij2 = 1.0;
        double expected_Aij = std::exp(std::pow(Rs_val * Rs_val - dij2, 2) / sigma_val) - 1;
        double expected_Dii = 2 * expected_Aij;
    
        // 检查 A 的非对角线数值
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                GiNaC::ex a_val = A(i, j).subs(subs).evalf();
                double a_numeric = GiNaC::ex_to<GiNaC::numeric>(a_val).to_double();
                if (i == j) {
                    EXPECT_NEAR(a_numeric, 0.0, 1e-8);
                } else {
                    EXPECT_NEAR(a_numeric, expected_Aij, 1e-6);
                }
            }
        }
    
        // 检查 D 的对角线数值
        for (int i = 0; i < N; ++i) {
            GiNaC::ex d_val = D(i, i).subs(subs).evalf();
            double d_numeric = GiNaC::ex_to<GiNaC::numeric>(d_val).to_double();
            EXPECT_NEAR(d_numeric, expected_Dii, 1e-6);
        }
    
        // 检查 L = D - A
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                GiNaC::ex l_val = L(i, j).subs(subs).evalf();
                double l_numeric = GiNaC::ex_to<GiNaC::numeric>(l_val).to_double();
                if (i == j) {
                    EXPECT_NEAR(l_numeric, expected_Dii, 1e-6);
                } else {
                    EXPECT_NEAR(l_numeric, -expected_Aij, 1e-6);
                }
            }
        }
    }
    
    TEST_F(LaplacianTest, TriangleEigenvaluesCheckAgainstMatlab) {
        N = 3;
        double Rs_val = 1.5;
        double sigma_val = 0.5;
        epsilon = 0.01;
    
        // Define regular triangle positions
        Eigen::MatrixXd tri(3, 2);
        tri << 0.0, 0.0,
               1.0, 0.0,
               0.5, 0.866;
    
        // Compute numerical Laplacian matrix and its eigenvalues
        Eigen::MatrixXd L_num;
        auto [lambda2, eigenvals] = getLambda2FromL(tri, L_num, Rs_val, sigma_val);
    
        // 1. Verify the Laplacian matrix is symmetric
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                EXPECT_NEAR(L_num(i, j), L_num(j, i), 1e-8);
    
        // 2. Compare all 3 eigenvalues against MATLAB reference values
        std::vector<double> reference_eigenvals = {
            0.00000000,    // First eigenvalue (should be zero)
            65.28469310,   // Second smallest eigenvalue (λ2)
            65.29473988    // Largest eigenvalue
        };
    
        ASSERT_EQ(eigenvals.size(), reference_eigenvals.size());
    
        for (int i = 0; i < N; ++i) {
            EXPECT_NEAR(eigenvals(i), reference_eigenvals[i], 1e-6)
                << "Mismatch in eigenvalue[" << i << "]";
        }
    
        // 3. Verify that lambda2 matches the second smallest eigenvalue
        EXPECT_NEAR(lambda2, reference_eigenvals[1], 1e-6);
        EXPECT_NEAR(lambda2, eigenvals(1), 1e-8);
    }
    
    