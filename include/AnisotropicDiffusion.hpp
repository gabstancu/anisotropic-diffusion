#ifndef ANISOTROPIC_DIFFUSION_HPP
#define ANISOTROPIC_DIFFUSION_HPP

#include "Eigen/Eigen"
#include <ginac/ginac.h>
#include <optional>
#include <variant>
#include <iostream>
#include "LinearSystem.hpp"
#include "BoundaryCondition.hpp"
#include "utils/SolverLog.hpp"

using SparseMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;
using DenseMatrix  = Eigen::MatrixXd;
using Vector       = Eigen::VectorXd;

class AnisotropicDiffusion
{
    private:

        /* Equation data */
        std::vector<GiNaC::symbol>                                      vars;
        GiNaC::ex                                                       analytical_expression{0.0};
        GiNaC::ex                                                       f{0.0};                     // source/sink
        std::pair<std::pair<double, double>, std::pair<double, double>> domain;
        DenseMatrix                                                     analytical_solution;
        double                                                          epsilon;
        double                                                          omega_;
        
        BoundaryConditions bc;
        
        /* Discretisation data */
        int         GRID_SIZE; // N + 2
        DenseMatrix grid;
        double      h;         // grid spacing

        /* Linear system */
        LinearSystem system;


        void   evaluate_analytical_solution  ();
        double evaluate_f                    (int i, int j);
        void   initialise_grid               ();
        void   initialise_linear_system      ();

        bool   has_source                    () { return !f.is_zero(); }   
        void   calc_omega_                   (); /* DO NOT USE */     

    public:
        AnisotropicDiffusion () = default;
        AnisotropicDiffusion (std::vector<GiNaC::symbol> vars, 
                              std::pair<std::pair<double, double>, std::pair<double, double>> domain, 
                              double eps, 
                              int GRID_SIZE);

        void print ();

        LinearSystem get_LinearSystem  ()               { return this->system; }
        void         set_source        (GiNaC::ex f)    { this->f = f; }
        void         set_analytical_u_ (GiNaC::ex _u_)  { this->analytical_expression = _u_;} 
        double       get_epsilon       ()               { return this->epsilon; }
        DenseMatrix  get_grid          ()               { return this->grid; }
        void         set_boundary      (Side side, Type type, std::variant<GiNaC::ex, double> expr);
        void         setup             ();
        void         fill_grid         (Vector u);
        void         save_grid         (std::string solver_name="");
};

#endif // ANISOTROPIC_DIFFUSION_HPP