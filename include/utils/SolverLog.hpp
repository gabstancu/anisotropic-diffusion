#ifndef SOLVER_LOG_HPP
#define SOLVER_LOG_HPP

#include "Eigen/Eigen"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <filesystem>
#include <chrono>
#include <limits>
#include <vector>
#include <string>
#include "utils/helper.hpp"

using SparseMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;
using DenseMatrix  = Eigen::MatrixXd;
using Vector       = Eigen::VectorXd;

class SolverLog
{
    private:

        int    num_of_iterations = 0;
        int    max_iterations    = 0;
        double tolerance         = 0.0;
        int    converged         = 0;        // 1/0
        int    timed_out         = 0;        // 1/0
        int    system_dim        = 0;

        double lambda_min_hat    = 0.0;
        double lambda_max_hat    = 0.0;
        double kappa_hat         = 0.0;


        double ic0_shift         = 0.0;


        int      n               = 0;        // interior grid size
        double   epsilon         = 0.0;      // anisotropy
        unsigned seed            = 0;

        std::vector<double>           res_per_iteration{};
        std::vector<double>           diff_per_iteration{};
        std::chrono::duration<double> time_elapsed{0};

        Vector                        final_solution;
        Vector                        direct_solution;

        std::string solver_name;
        std::string precon_name;

    public:
    
        void set_labels       (std::string solver, std::string precond) { solver_name    = std::move(solver); precon_name = std::move(precond); }
        void set_grid         (int Ninner)                              { n              = Ninner; }
        void set_epsilon      (double eps)                              { epsilon        = eps; }
        void set_tol          (double tol)                              { tolerance      = tol; }
        void set_max_iters    (int m)                                   { max_iterations = m; }
        void set_dim          (int dim)                                 { system_dim     = dim; }
        void set_seed         (unsigned s)                              { seed           = s; }
        void set_time_elapsed (std::chrono::duration<double> time)      { time_elapsed   = time; }
        void set_ic0_shift    (double s)                                { ic0_shift      = s; }
        void set_lanczos      (double lmin, double lmax)                { lambda_min_hat = lmin; lambda_max_hat = lmax; kappa_hat = (lmin>0)? (lmax/lmin) : std::numeric_limits<double>::infinity(); }

        void inc_iter      ()        { ++num_of_iterations; }
        void mark_converged(bool ok) { converged = ok ? 1 : 0; }
        void mark_timed_out(bool to) { timed_out = to ? 1 : 0; }

        void push_residual(double rrel) { res_per_iteration.push_back(rrel); }


        void set_final_solution (const Vector& x) { final_solution  = x; }
        void set_direct_solution(const Vector& x) { direct_solution = x; }

        void print() {
            if (!solver_name.empty())   std::cout << "\nSolver: "     << solver_name  << '\n';
            if (!precon_name.empty())   std::cout << "Precon.: "    << precon_name  << '\n';
            std::cout << "System dim: " << system_dim                                                    << '\n'
                    << "n: "          << n                        << "  eps: "       << epsilon        << '\n'
                    << "iters: "      << num_of_iterations        << " / "           << max_iterations << '\n'
                    << "converged: "  << converged                << "  timed_out: " << timed_out      << '\n'
                    << "time: "       << time_elapsed.count()*1e3 << " ms\n";
        }

    void log ()
    {
        namespace fs = std::filesystem;

        // paths
        const std::string log_path  = get_current_working_directory() + "/logs/";
        const std::string file_path = log_path + "results.csv";

        // ensure directory
        fs::create_directories(log_path);

        const bool new_file = !fs::exists(file_path);

        std::ofstream file(file_path, std::ios::app);
        if (!file.is_open())
            throw std::runtime_error("[LOG] Failed to open: " + file_path);

        if (new_file) {
            file << "n,system_dim,epsilon,solver,tol,max_iters,iterations,time_ms,res_rel,"
                    "lambda_min_hat,lambda_max_hat,kappa_hat,precond,ic0_shift,seed,converged,timed_out\n";
        }

        const double time_ms       = time_elapsed.count() * 1e3;
        const double res_rel_final = res_per_iteration.empty()
                                     ? std::numeric_limits<double>::quiet_NaN()
                                     : res_per_iteration.back();

        file << std::setprecision(8) << std::scientific
             << n                 << ","
             << system_dim        << ","
             << epsilon           << ","
             << (solver_name.empty() ? std::string("unknown") : solver_name) << ","
             << tolerance         << ","
             << max_iterations    << ","
             << num_of_iterations << ","
             << time_ms           << ","
             << res_rel_final     << ","
             << lambda_min_hat    << ","
             << lambda_max_hat    << ","
             << kappa_hat         << ","
             << (precon_name.empty() ? std::string("none") : precon_name) << ","
             << ic0_shift         << ","
             << seed              << ","
             << converged         << ","
             << timed_out         
             << "\n";
    }
};

#endif // SOLVER_LOG_HPP