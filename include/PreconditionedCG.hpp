#ifndef PRECONDITIONED_CG_HPP
#define PRECONDITIONED_CG_HPP

#include "Eigen/Eigen"
#include "utils/SolverLog.hpp"
#include <iostream>

using SparseMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;
using DenseMatrix  = Eigen::MatrixXd;
using Vector       = Eigen::VectorXd;

class PreconditionedCG
{
    private:
        double      tol             = 1e-6;
        int         max_iters       = 1e6;
        std::string name            = "PCG";
        Vector      final_solution;

    public:
        PreconditionedCG() = default;
        PreconditionedCG(double tolerance, int max_iters);

        template <class System, class Preconditioner>
        void solve (System& system, Preconditioner& M, SolverLog& logger)
        {
            const auto& A = system.A;
            const auto& b = system.b;
                  auto& u = system.u;

            logger.set_tol(this->tol);
            logger.set_max_iters(this->max_iters); 
            logger.set_labels(this->name, "IC0");
            logger.set_dim(A.rows());
            

            Vector r(b.size()), zeta(b.size()), d(b.size());
            Vector Ad(A.rows());

            M.compute(A);                                 

            r.noalias()    = b - A * u;
            M.apply(r, zeta);                            
            d.noalias()    = zeta;

            double b_norm = b.norm();

            if (r.norm() / b_norm < tol) 
            {   
                this->final_solution = u;
                logger.set_final_solution(u);
                logger.mark_converged(true);
                return;
            }

            for (int k = 0; k < max_iters; k++)
            {   
                Ad.noalias()     = A * d;
                Vector r_prev    = r;
                Vector zeta_prev = zeta;

                double alpha = r.dot(zeta) / d.dot(Ad);
                u.noalias()     += alpha * d;
                r.noalias()     -= alpha * Ad;

                logger.inc_iter();
                logger.push_residual(r.norm() / b_norm );

                if (r.norm() / b_norm <= tol) 
                {   
                    this->final_solution = u;
                    logger.set_final_solution(u);
                    logger.mark_converged(true);
                    return;
                }

                M.apply(r, zeta);                         
                double beta = r.dot(zeta) / r_prev.dot(zeta_prev);
                d.noalias()  = zeta + beta * d;
            }
            this->final_solution = u;
            logger.set_final_solution(u);
            return;            
        }

        void print         ();
        void set_tol       (double tolerance) { this->tol       = tolerance; }
        void set_max_iters (int iters)        { this->max_iters = iters;}
};

#endif // PRECONDITIONED_CG_HPP