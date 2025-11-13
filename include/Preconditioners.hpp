// #ifndef PRECONDITIONERS_HPP
// #define PRECONDITIONERS_HPP

// #include <Eigen/Sparse>
// #include <unsupported/Eigen/IterativeSolvers>

// using SparseMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;
// using Vector       = Eigen::VectorXd;

// /* slower ????? */
// struct IC0
// {
//     Eigen::IncompleteCholesky<double> ic;

//     void compute (const SparseMatrix& A)
//     {
//         ic.setInitialShift(1e-10);
//         ic.compute(A);
//         if (ic.info() != Eigen::Success) 
//         {
//             throw std::runtime_error("IC(0) factorisation failed");
//         }
//     }

//     void apply(const Vector& r, Vector& z) const 
//     {
//         z = ic.solve(r);                
//     }
// };

// #endif // PRECONDITIONERS_HPP



#ifndef PRECONDITIONERS_HPP
#define PRECONDITIONERS_HPP

#include <Eigen/Sparse>
#include <cmath>
#include <vector>
#include <stdexcept>

using RowMat = Eigen::SparseMatrix<double, Eigen::RowMajor>;
using ColMat = Eigen::SparseMatrix<double>;
using Vector = Eigen::VectorXd;


struct IC0
{
    RowMat L;                               // lower-triangular factor
    std::vector<std::vector<int>> colRows;  // pattern info (rows i>j)
    double damping = 1e-12;                 // stabilization constant

    // Dot product between row i and row j restricted to columns < j
    static double dot_rows_lt(const RowMat& M, int i, int j)
    {
        double s = 0.0;
        RowMat::InnerIterator it_i(M, i), it_j(M, j);

        while (it_i && it_j && it_i.col() < j && it_j.col() < j) 
        {
            const int ci = it_i.col();
            const int cj = it_j.col();

            if (ci == cj) 
            { 
                s += it_i.value() * it_j.value(); 
                ++it_i; 
                ++it_j; 
            }
            else if (ci  < cj) 
            { 
                ++it_i; 
            }
            else 
            { 
                ++it_j; 
            }
        }
        return s;
    }

    // Factorization: build L for M \approx A^{-1}
    void compute(const RowMat& A_row, double damping_eps = 1e-12)
    {
        damping = damping_eps;

        RowMat A_rm = A_row;
        A_rm.makeCompressed();
        ColMat A = A_rm; // convert to column-major

        const int n = static_cast<int>(A.rows());
        if (A.rows() != A.cols())
            throw std::invalid_argument("IC0::compute: A must be square.");

        // Copy diagonal + lower part from A
        std::vector<Eigen::Triplet<double>> trips;
        trips.reserve(A.nonZeros());
        for (int k = 0; k < A.outerSize(); ++k) 
        {
            for (ColMat::InnerIterator it(A, k); it; ++it) 
            {
                const int i = it.row();
                const int j = it.col();

                if (i >= j) 
                    trips.emplace_back(i, j, it.value());
            }
        }

        L.resize(n, n);
        L.setFromTriplets(trips.begin(), trips.end());
        L.makeCompressed();

        // Build pattern: for each col j, which rows i>j exist
        colRows.assign(n, {});
        for (int i = 0; i < n; ++i) 
        {
            for (RowMat::InnerIterator it(L, i); it; ++it) 
            {
                const int j = it.col();
                if (i > j) 
                    colRows[j].push_back(i);
            }
        }

        // ---------------- IC(0) factorization ----------------
        for (int j = 0; j < n; ++j)
        {
            const double Ajj0 = L.coeff(j, j);
                  double s    = Ajj0;

            for (RowMat::InnerIterator it(L, j); it; ++it) 
            {
                const int p = it.col();
                if (p >= j) 
                    break;
                const double Ljp = it.value();
                s -= Ljp * Ljp;
            }

            if (s <= 0.0) 
            {
                const double base = std::abs(Ajj0);
                s = std::max(base * damping, damping);
            }

            const double Ljj = std::sqrt(s);
            L.coeffRef(j, j) = Ljj;

            for (int i : colRows[j]) 
            {
                const double Aij0 = L.coeff(i, j);
                const double dot  = dot_rows_lt(L, i, j);
                L.coeffRef(i, j)  = (Aij0 - dot) / Ljj;
            }
        }
        L.makeCompressed();
    }

    // Apply: z = M^{-1} r  with M = L * L^{T}
    void apply(const Vector& r, Vector& z) const
    {
        Vector y = L.template triangularView<Eigen::Lower>().solve(r);
               z = L.transpose().template triangularView<Eigen::Upper>().solve(y);
    }

    // Convenience overload returning the vector
    Vector apply(const Vector& r) const
    {
        Vector z(r.size());
        apply(r, z);
        return z;
    }
};

#endif // PRECONDITIONERS_HPP