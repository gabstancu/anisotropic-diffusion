#include "PreconditionedCG.hpp"

PreconditionedCG::PreconditionedCG (double tolerance, int max_iters)
{
    this->tol       = tolerance;
    this->max_iters = max_iters;
}


void PreconditionedCG::print ()
{
    std::cout << "Solver: "     << this->name      << '\n';
    std::cout << "Tolerance: "  << this->tol       << '\n';
    std::cout << "Max. iters: " << this->max_iters << '\n';
}