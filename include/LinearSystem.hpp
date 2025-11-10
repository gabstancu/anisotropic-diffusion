#ifndef LINEAR_SYSTEM_HPP
#define LINEAR_SYSTEM_HPP

#include "Eigen/Eigen"
#include "utils/SolverLog.hpp"
#include <iostream>
#include <random>
#include <chrono>
#include <cmath>

using SparseMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;
using Vector       = Eigen::VectorXd;

class LinearSystem
{
    private:

        double    lambda_min = 0.0;
        double    lambda_max = 0.0;
        double    kappa_hat  = 0.0;

    public:

        LinearSystem() = default;

        LinearSystem (SparseMatrix A, Vector b);

        SparseMatrix A;
        Vector       b;
        Vector       u;
        Vector      _u_;
        int          N{0};
        double       omega_{1.0};

        SparseMatrix get_A              () const { return this->A; }
        Vector       get_b              () const { return this->b; }
        Vector       get_u              () const { return this->u; }
        Vector       get_exact_solution () const { return this->_u_; }
        double       get_omega_         () const { return this->omega_; }


        template <class Solver>
        void solve (Solver& solver, SolverLog& logger)
        {   
            auto start = std::chrono::high_resolution_clock::now();
            solver.solve(*this, logger);
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;

            logger.set_time_elapsed(elapsed);
            logger.set_dim(int(std::pow(this->N, 2)));
            logger.set_grid(this->N);
            logger.set_lanczos(this->lambda_min, this->lambda_max);
        }

        template <class Solver, class Preconditioner>
        void solve (Solver& solver, Preconditioner& preconditioner, SolverLog& logger)
        {   
            auto start = std::chrono::high_resolution_clock::now();
            solver.solve(*this, preconditioner, logger);
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;

            logger.set_time_elapsed(elapsed);
            logger.set_dim(int(std::pow(this->N, 2)));
            logger.set_grid(this->N);
            logger.set_lanczos(this->lambda_min, this->lambda_max);
        }

        void estimate_extremes_eigs (int      max_steps = 20,
                                     unsigned seed      = 0);

        void solve_directly ();
        void reset_solution ();
        void print          ();
        
        /* {l_max, l_min} */
        std::pair<double, double> get_lambdas    () { return std::make_pair(this->lambda_max, this->lambda_min); }
        double                    get_kappa_hat  () { return this->kappa_hat; }
    };

#endif // LINEAR_SYSTEM_HPP