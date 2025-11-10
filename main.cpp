#include "Eigen/Eigen"
#include <string>
#include <ginac/ginac.h>
#include <iostream>

#include "AnisotropicDiffusion.hpp"
#include "LinearSystem.hpp"
#include "ConjugateGradient.hpp"
#include "PreconditionedCG.hpp"
#include "Preconditioners.hpp"


std::vector<double> epsilon  = {
    1.0, 7e-1, 5e-1, 3e-1, 1e-1, 7e-2, 5e-2, 3e-2, 1e-2,
    7e-3, 5e-3, 3e-3, 1e-3, 7e-4, 5e-4, 3e-4, 1e-4,
    3e-5, 1e-5, 1e-6
};

std::vector<int>    GRID_DIM = {66, 130, 258, 514};


void eval_PCG ()
{
    std::vector<GiNaC::symbol> vars;
    GiNaC::symbol              x_("x"), y_("y"); 
    vars.push_back(x_); 
    vars.push_back(y_);
    
    std::pair<std::pair<double, double>, std::pair<double, double>> domain;
    domain = std::make_pair(std::make_pair(0.0, 1.0), std::make_pair(0.0, 1.0));

    double alpha_val = 0.8;   // top scale
    double b0        = 1.2;   // source amplitude

    GiNaC::ex f = b0 * GiNaC::sin(GiNaC::Pi*x_) * GiNaC::sin(GiNaC::Pi*y_);

    for (double e : epsilon)
    {   
        std::cout << "------------------------- epsilon: " << e << "-------------------------\n";
        for (int N_out : GRID_DIM)
        {   
            std::cout << "\nN: " << N_out << '\n';
            SolverLog logger;
            AnisotropicDiffusion AD_2D (vars, domain, e, N_out); 
            AD_2D.set_source(f);
            AD_2D.set_boundary(Side::Left,   Type::Dirichlet, 4.0);
            AD_2D.set_boundary(Side::Right,  Type::Dirichlet, 4.0);
            AD_2D.set_boundary(Side::Bottom, Type::Dirichlet, GiNaC::sin(GiNaC::Pi*x_));
            AD_2D.set_boundary(Side::Top,    Type::Dirichlet, alpha_val * GiNaC::sin(GiNaC::Pi*x_));
            AD_2D.setup();
            logger.set_epsilon(AD_2D.get_epsilon());

            LinearSystem system = AD_2D.get_LinearSystem();

            PreconditionedCG PCG;
            IC0 ic_0;
            system.solve(PCG, ic_0, logger);
            logger.print();
            logger.log();
        }
    }     
}


void eval_CG ()
{   
    std::vector<GiNaC::symbol> vars;
    GiNaC::symbol              x_("x"), y_("y"); 
    vars.push_back(x_); 
    vars.push_back(y_);

    std::pair<std::pair<double, double>, std::pair<double, double>> domain;
    domain = std::make_pair(std::make_pair(0.0, 1.0), std::make_pair(0.0, 1.0));

    double alpha_val = 0.8;   // top scale
    double b0        = 1.2;   // source amplitude

    GiNaC::ex f = b0 * GiNaC::sin(GiNaC::Pi*x_) * GiNaC::sin(GiNaC::Pi*y_);
    
    for (double e : epsilon)
    {   
        std::cout << "------------------------- epsilon: " << e << "-------------------------\n";       
        for (int N_out : GRID_DIM)
        {   
            std::cout << "\nN: " << N_out << '\n';
            SolverLog logger;
            AnisotropicDiffusion AD_2D (vars, domain, e, N_out); 
            AD_2D.set_source(f);
            AD_2D.set_boundary(Side::Left,   Type::Dirichlet, 4.0);
            AD_2D.set_boundary(Side::Right,  Type::Dirichlet, 4.0);
            AD_2D.set_boundary(Side::Bottom, Type::Dirichlet, GiNaC::sin(GiNaC::Pi*x_));
            AD_2D.set_boundary(Side::Top,    Type::Dirichlet, alpha_val * GiNaC::sin(GiNaC::Pi*x_));
            AD_2D.setup();
            logger.set_epsilon(AD_2D.get_epsilon());

            LinearSystem system = AD_2D.get_LinearSystem();
            system.estimate_extremes_eigs();
            std::pair<double, double> p = system.get_lambdas();
            logger.set_lanczos(p.second, p.first);
            
            ConjugateGradient CG;
            system.solve(CG, logger);
            logger.print();
            logger.log();
        }
    } 
}


int main ()
{     
    // std::vector<double> eps = {1.0, 4e-1, 2e-1, 4e-2, 2e-2, 1e-3};

    // std::vector<GiNaC::symbol> vars;
    // GiNaC::symbol              x_("x"), y_("y");    
    // vars.push_back(x_); 
    // vars.push_back(y_);

    // double alpha_val = 0.8;   // top scale
    // double b0        = 1.2;   // source amplitude

    // std::pair<std::pair<double, double>, std::pair<double, double>> domain;
    // domain = std::make_pair(std::make_pair(0.0, 1.0), std::make_pair(0.0, 1.0));

    // for (double e : eps)
    // {
    //     // GiNaC::ex analytical_solution = GiNaC::sin(GiNaC::Pi * x_) * GiNaC::cos(GiNaC::Pi * y_) + e * GiNaC::sin(2 * GiNaC::Pi * x_) * GiNaC::cos(GiNaC::Pi * y_);
    //     // GiNaC::ex f = GiNaC::pow(GiNaC::Pi, 2) * GiNaC::cos(GiNaC::Pi * y_) * ( (1 + e) * GiNaC::sin(GiNaC::Pi * x_) + e * (4 + e) * GiNaC::sin(2 * GiNaC::Pi * x_) );
    //     // ex k = Pi / sqrt(e);
    //     // ex u_BC = sin(Pi*x_) * (sinh(k*(1 - y_)) + alpha_val*sinh(k*y_)) / sinh(k);
    //     // ex w = (b0)/(pow(Pi,2)*(1 + e)) * sin(Pi*x_) * sin(Pi*y_);
    //     // ex analytical_solution = u_BC + w;
    //     ex f                   = b0 * sin(Pi*x_) * sin(Pi*y_);

    //     AnisotropicDiffusion AD_2D (vars, domain, e, 514);     
    //     AD_2D.set_source(f);
    //     AD_2D.set_boundary(Side::Left,   Type::Dirichlet, 4.0);
    //     AD_2D.set_boundary(Side::Right,  Type::Dirichlet, 4.0);
    //     AD_2D.set_boundary(Side::Bottom, Type::Dirichlet, sin(Pi*x_));
    //     AD_2D.set_boundary(Side::Top,    Type::Dirichlet, alpha_val * sin(Pi*x_));
    //     // AD_2D.set_boundary(Side::Bottom, Type::Dirichlet, GiNaC::sin(GiNaC::Pi * x_) + e * GiNaC::sin(2 * GiNaC::Pi * x_));
    //     // AD_2D.set_boundary(Side::Top,    Type::Dirichlet, -GiNaC::sin(GiNaC::Pi * x_) - e * GiNaC::sin(2 * GiNaC::Pi * x_));
    //     AD_2D.setup();

    //     LinearSystem system = AD_2D.get_LinearSystem();

    //     PreconditionedCG PCG; SolverLog logger;
    //     IC0 ic_0;
    //     system.solve(PCG, ic_0, logger);
    //     AD_2D.fill_grid(system.u);
    //     AD_2D.save_grid("PCG");
    //     logger.print();
    //     // logger.log();
    // }

   
    eval_CG();
    eval_PCG();


    

    // system.estimate_extremes_eigs();
    // logger_CG.set_epsilon(AD_2D.get_epsilon());
    // logger_PCG.set_epsilon(AD_2D.get_epsilon());
    // // std::pair<double, double> p = system.get_lambdas();
    // // std::cout << "l. max: "    << p.first                << " l. min: " << p.second << '\n';
    // // std::cout << "kappa hat: " << system.get_kappa_hat()                            << '\n';

    // // std::cout << "A:\n"    << system.A << '\n';
    // // std::cout << "b:\n"    << system.b << '\n';
    // // std::cout << "u:\n"    << system.u << '\n';
    // // std::cout << "grid:\n" << AD_2D.get_grid() << '\n';


    // system.reset_solution();

    // PreconditionedCG PCG;
    // IC0 ic_0;
    // system.solve(PCG, ic_0, logger_PCG);
    // logger_PCG.print();
    // logger_PCG.log();

    return 0;
}