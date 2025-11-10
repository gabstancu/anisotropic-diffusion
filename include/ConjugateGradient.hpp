#ifndef CONJUGATE_GRADIENT_HPP
#define CONJUGATE_GRADIENT_HPP

#include "Eigen/Eigen"
#include "utils/SolverLog.hpp"

using SparseMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;
using DenseMatrix  = Eigen::MatrixXd;
using Vector       = Eigen::VectorXd;

class ConjugateGradient
{
    private:
        double      tol             = 1e-6;
        int         max_iters       = 1e6;
        std::string name            = "CG";
        Vector      final_solution;
        SolverLog   logger;

    public:

        ConjugateGradient() = default;
        ConjugateGradient(double tolerance, int max_iters);

        template<class System>
        void solve (System& system, SolverLog& logger)
        {
            const auto& A = system.A;
            const auto& b = system.b;
                  auto& u = system.u;

            logger.set_tol(this->tol);
            logger.set_max_iters(this->max_iters); 
            logger.set_labels(this->name, "");
            logger.set_dim(A.rows());

            Vector r = b - A * u;
            double b_norm = b.norm();
            double r_norm = r.norm();

            if (r_norm / b_norm <= tol) 
            {   
                this->final_solution = u;
                logger.set_final_solution(u);
                logger.mark_converged(true);
                return;
            }

            Vector d = r; // initial search direction
            Vector Ad(A.rows());

            for (int k = 0; k < max_iters; k++)
            {   
                // std::cout << "--------------------- iter. " << k+1 << " ---------------------\n";
                Ad.noalias() = A * d;

                double alpha      = ((r.transpose() * r) / (d.transpose() * Ad)).coeff(0); // step size
                double r_prev_dot = (r.transpose() * r).coeff(0); // to calculate beta

                u.noalias() += alpha * d;
                r.noalias() -= alpha * Ad;

                r_norm = r.norm();

                logger.inc_iter();
                logger.push_residual(r_norm / b_norm );

                if (r_norm / b_norm <= tol) 
                {   
                    this->final_solution = u;
                    logger.set_final_solution(u);
                    logger.mark_converged(true);
                    return;
                }

                double beta = r.dot(r) / r_prev_dot;
                d           = r + beta * d; // update direction
            }
            this->final_solution = u;
            logger.set_final_solution(u);
            return;
        }

        void print         ();
        void set_tol       (double tolerance) { this->tol       = tolerance; }
        void set_max_iters (int iters)        { this->max_iters = iters;}
        
};

#endif // CONJUGATE_GRADIENT_HPP