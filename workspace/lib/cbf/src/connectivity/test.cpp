#include <ginac/ginac.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <ginac/ginac.h>
#include <fstream>
#include <iostream>

// 假设 N、Rs、sigma、epsilon、px、py、current_robot_positions 是外部可访问变量
int N;
GiNaC::ex Rs, sigma, epsilon;
std::vector<GiNaC::symbol> px, py;
Eigen::MatrixXd current_robot_positions;

std::tuple<GiNaC::matrix, GiNaC::matrix, GiNaC::matrix> computeLaplacian()
{
    GiNaC::matrix A(N, N);
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            if (i != j)
            {
                GiNaC::ex dij = GiNaC::sqrt(GiNaC::pow(px[i] - px[j], 2) + GiNaC::pow(py[i] - py[j], 2));
                GiNaC::ex expr = GiNaC::exp(GiNaC::pow((Rs * Rs - dij * dij) , 2) / sigma) - 1;
                A(i, j) = expr;
            }
            else
            {
                A(i, j) = 0;
            }
        }
    }
    GiNaC::matrix D(N, N);
    for (int i = 0; i < N; ++i)
    {
        GiNaC::ex sum_aij = 0;
        for (int j = 0; j < N; ++j)
        {
            sum_aij += A(i, j);
        }
        D(i, i) = sum_aij;
    }

    GiNaC::matrix L(N, N);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            L(i, j) = D(i, j) - A(i, j);

    return {A, D, L};
}

std::pair<double, Eigen::VectorXd> getLambda2FromL(
    const Eigen::MatrixXd& robot_positions,
    Eigen::MatrixXd& L_num_out,
    double Rs_value,
    double sigma_value)
{
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
    L_num_out = L_num;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(L_num);
    Eigen::VectorXd eigenvals = solver.eigenvalues();
    return {eigenvals(1), eigenvals};
}



void run_experiment(const std::string& name, const Eigen::MatrixXd& positions, double Rs_value, double sigma_value)

{
    N = positions.rows();
    px.clear(); py.clear();
    for (int i = 0; i < N; ++i)
    {
        px.emplace_back(GiNaC::symbol("px" + std::to_string(i)));
        py.emplace_back(GiNaC::symbol("py" + std::to_string(i)));
    }

    auto [A, D, L_symbolic] = computeLaplacian();

    Eigen::MatrixXd L_num;
    auto [lambda2, eigenvalues] = getLambda2FromL(positions, L_num, Rs_value, sigma_value);

    GiNaC::ex b5 = lambda2 - epsilon;

    std::ofstream out(name + ".txt");
    out << "=== Experiment: " << name << " ===\n";
    out << "Robot positions:\n" << positions << "\n\n";
    out << "=== Matrix A (Adjacency) ===\n" << A << "\n";
    out << "=== Matrix D (Degree) ===\n" << D << "\n";
    out << "=== Matrix L (Symbolic Laplacian) ===\n" << L_symbolic << "\n";
    out << "=== Matrix L_num (Numeric Laplacian) ===\n" << L_num << "\n";
    out << "=== Eigenvalues ===\n" << eigenvalues.transpose() << "\n";
    out << "=== lambda2 ===\n" << lambda2 << "\n";
    out << "=== CBF b5 = lambda2 - epsilon ===\n" << b5 << "\n";
    out.close();
}

int main()
{
    double Rs_value = 1.5; 
    Rs = Rs_value;
    double sigma_value = 0.5;    // 数值化的 sigma
    sigma = sigma_value;
    epsilon = 0.01;

    // 正三角形
    Eigen::MatrixXd tri(3, 2);
    tri << 0.0, 0.0,
           1.0, 0.0,
           0.5, 0.866;

    // 正方形
    Eigen::MatrixXd square(4, 2);
    square << 0.0, 0.0,
              2.0, 0.0,
              0.0, 2.0,
              2.0, 2.0;

    // 正五边形
    Eigen::MatrixXd pentagon(5, 2);
    double angle = 2 * M_PI / 5;
    for (int i = 0; i < 5; ++i)
    {
        pentagon(i, 0) = std::cos(i * angle);
        pentagon(i, 1) = std::sin(i * angle);
    }

    run_experiment("triangle", tri, Rs_value, sigma_value);
    run_experiment("square", square, Rs_value, sigma_value);
    run_experiment("pentagon", pentagon, Rs_value, sigma_value);
    

    return 0;
}
