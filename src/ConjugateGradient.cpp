#include "ConjugateGradient.hpp"

ConjugateGradient::ConjugateGradient (double tolerance, int max_iters)
{
    this->tol       = tolerance;
    this->max_iters = max_iters;
    
    this->logger.set_tol(this->tol);
    this->logger.set_max_iters(this->max_iters);
    this->logger.set_labels(this->name, "");
}

void ConjugateGradient::print ()
{

}