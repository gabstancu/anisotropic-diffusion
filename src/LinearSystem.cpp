#include "LinearSystem.hpp"

LinearSystem::LinearSystem (SparseMatrix A_, Vector b_) : A(std::move(A_)), b(std::move(b_))
{
    if (A.rows() != b.size())
    {
        throw std::invalid_argument("Dimension mismatch in LinearSystem.\n");
    }
    this-> _u_.setZero();
    this->  u .setZero();
    this->  N = static_cast<int>(A.rows());
}

void LinearSystem::solve_directly ()
{
    Eigen::SparseLU<SparseMatrix> solver;
    solver.compute(this->A);

    if (solver.info() != Eigen::Success)
    {   
        throw std::runtime_error("Sparse decomposition failed");
    }
    this->_u_ = solver.solve(this->b);
}

void LinearSystem::estimate_extremes_eigs (int max_steps, unsigned seed)
{
    int n = A.rows();
    if (n == 0)
    {
        throw std::runtime_error("Empty matrix in Lanczos.");
    }

    std::mt19937 gen(seed);
    std::uniform_real_distribution<> dist(-1.0, 1.0);

    Vector q = Vector::NullaryExpr(n, [&]() { return dist(gen); } );
    q.normalize();

    Vector q_prev = Vector::Zero(n);
    Vector Aq     = Vector::Zero(n);

    std::vector<double> alpha, beta;
    alpha.reserve(max_steps);
    beta.reserve(max_steps);

    for (int j = 0; j < max_steps; ++j)
    {
        Aq.noalias() = A * q;

        double a_j = q.dot(Aq);
        alpha.push_back(a_j);

        if (j > 0) 
            Aq.noalias() -= beta.back() * q_prev;
        Aq.noalias() -= a_j * q;

        double b_j = Aq.norm();
        if (b_j < 1e-14) break;

        beta.push_back(b_j);
        q_prev = q;
        q      = Aq / b_j;
    }

    // build the tridiagonal T
    const int m = static_cast<int>(alpha.size());
    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(m, m);

    for (int i = 0; i < m; ++i) 
    {
        T(i,i) = alpha[i];
        // only set off-diagonals if both indices are inside the matrix
        if (i+1 < m && i < static_cast<int>(beta.size())) {
            T(i, i+1) = beta[i];
            T(i+1, i) = beta[i];
        }
    }

    // After the Lanczos loop:
    if (m == 0) {
        throw std::runtime_error("Lanczos produced no alpha; check input matrix or max_steps.");
    }
    if (m == 1) {
        // With a single step, T is [alpha1]; both extremes are the same
        this->lambda_min = this->lambda_max = alpha[0];
        return;
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigs(T);
    this->lambda_min = eigs.eigenvalues()(0);
    this->lambda_max = eigs.eigenvalues()(m-1);
    this->kappa_hat  = this->lambda_max / this->lambda_min; 
    
}

void LinearSystem::reset_solution ()
{
    this->u.setZero();
}

void LinearSystem::print ()
{
    std::cout << "System dimension: " << this->N      << '\n';
    std::cout << "omega_: "           << this->omega_ << '\n';
}

