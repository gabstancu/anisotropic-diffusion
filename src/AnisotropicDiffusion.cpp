#include "AnisotropicDiffusion.hpp"

AnisotropicDiffusion::AnisotropicDiffusion (std::vector<GiNaC::symbol> vars, 
                                            std::pair<std::pair<double, double>, std::pair<double, double>> domain, 
                                            double epsilon, 
                                            int GRID_SIZE)
{
    this->bc.bottom.expr = 0.0;
    this->bc.top.expr    = 0.0;
    this->bc.left.expr   = 0.0;
    this->bc.right.expr  = 0.0;

    this->vars                  = vars;
    this->epsilon               = epsilon;
    this->domain                = domain;
    this->GRID_SIZE             = GRID_SIZE;
    this->h                     = 1.0 / (GRID_SIZE - 2 + 1.0);

    // calc_omega_();

    // this->f = source;
}

void AnisotropicDiffusion::setup ()
{
    initialise_grid();
    evaluate_analytical_solution();
    initialise_linear_system();
}


void AnisotropicDiffusion::calc_omega_ () /* DO NOT USE  */
{
    this->omega_ = 2.0 / (1 + std::sin((M_PI) / this->GRID_SIZE-2 ));
    this->omega_ = std::min(this->omega_, 1.9 + (0.05) / this->GRID_SIZE-1);       
}


void AnisotropicDiffusion::evaluate_analytical_solution ()
{   
    if (analytical_expression.is_zero()) /* the expression of the analytical solution has not been specified */
    {
        return;
    }

    analytical_solution.resize(GRID_SIZE, GRID_SIZE);
            
    double x_ = 0.0, y_ = 0.0;

    for (int i = 0; i < GRID_SIZE; i++)
    {  
        y_ = h * (GRID_SIZE - 1 - i);
        for (int j = 0; j < GRID_SIZE; j++)
        {
            x_ = h * j;
            analytical_solution(i, j) = grid(i, j);

            if ((i>=1 && i <= GRID_SIZE - 2) && (j>=1 && j <= GRID_SIZE - 2))
            {
                GiNaC::exmap m;
                m[vars[0]] = x_; 
                m[vars[1]] = y_;
                GiNaC::ex evaluated_expression  = analytical_expression.subs(m).evalf();
                double    value                 = GiNaC::ex_to<GiNaC::numeric>(evaluated_expression).to_double();
                analytical_solution(i, j)       = value;
            }
        }
    }
}


double AnisotropicDiffusion::evaluate_f (int i, int j)
{
    if (!has_source())
        return 0.0;
    
    double x_ = h * j, y_ = h * (GRID_SIZE - 1 - i);

    GiNaC::exmap m;
    m[vars[0]] = x_;
    m[vars[1]] = y_;
    GiNaC::ex evaluated_expression = f.subs(m).evalf();
    return GiNaC::ex_to<GiNaC::numeric>(evaluated_expression).to_double();
}


void AnisotropicDiffusion::initialise_grid ()
{
    // int N = GRID_SIZE - 2;

    grid.resize(GRID_SIZE, GRID_SIZE);

    double boundary_val;

    for (int i = 0; i < GRID_SIZE; i++)
    {
        for (int j = 0; j < GRID_SIZE; j++)
        {
            grid(i, j) = 0.0;

            // top boundary
            if (i == 0)
            {   
                grid(i, j) = bc.top.evaluate(std::make_pair(j, GRID_SIZE - 1 - i), vars, h);
            }

            // bottom boundary
            if (i == GRID_SIZE - 1)
            {
                grid(i, j) = bc.bottom.evaluate(std::make_pair(j, GRID_SIZE - 1 - i), vars, h);
            }

            // left boundary
            if (j == 0 && (i > 0 && i < GRID_SIZE - 1)) // keep top value
            {
                grid(i, j) = bc.left.evaluate(std::make_pair(j, GRID_SIZE - 1 - i), vars, h);
            }

            // right boundary
            if (j == GRID_SIZE - 1 && (i > 0 && i < GRID_SIZE - 1))
            {
                grid(i, j) = bc.right.evaluate(std::make_pair(j, GRID_SIZE - 1 - i), vars, h);
            }
        }
    }
}


void AnisotropicDiffusion::set_boundary (Side side, Type type, std::variant<GiNaC::ex, double> expr)
{
    if (side == Side::Top)
    {
        this->bc.top.expr = expr;
        this->bc.top.side = side;
        this->bc.top.type = type;
    }
    else if (side == Side::Bottom)
    {
        this->bc.bottom.expr = expr;
        this->bc.bottom.side = side;
        this->bc.bottom.type = type;
    }
    else if (side == Side::Left)
    {
        this->bc.left.expr = expr;
        this->bc.left.side = side;
        this->bc.left.type = type;
    }
    else
    {
       this->bc.right.expr = expr;
       this->bc.right.side = side;
       this->bc.right.type = type;
    }
}


void AnisotropicDiffusion::fill_grid (Vector u)
{
    int N = grid.rows();
    int K = u.size();
    int i, j;

    for (int k = 0; k < K; k++)
    {
        i = k / (N-2) + 1;
        j = k % (N-2) + 1;
        grid(i, j) = u(k);
    }
}


void AnisotropicDiffusion::initialise_linear_system ()
{
    /* - u_{xx} - epsilon * u_{yy]} = f */
    int N   = GRID_SIZE - 2;
    int dim = N * N;

    system.A.resize(dim, dim);
    system.A.reserve(Eigen::VectorXi::Constant(dim, 5));

    system.b = Eigen::VectorXd::Zero(dim);
    system.u = Eigen::VectorXd::Zero(dim);
    system.N = N;

    double cC =  2.0 * (1.0 + epsilon);      // center
    double cX = -1.0;                        // west/east
    double cY = -epsilon;                    // south/north

    int k;

    for (int i = 1; i <= N; i++)
    {
        for (int j = 1; j <= N; j++)
        {
            k = (i - 1) * N + (j - 1);

            system.A.insert(k, k) = cC;

            system.b(k) += h * h * evaluate_f(i, j);

            if (i == N) // top neighbor
                system.b(k) += (-cY) * grid(i + 1, j);
            else
                system.A.insert(k, k + N) = cY;

            if (i == 1) // bottom neighbor
                system.b(k) += (-cY) * grid(i - 1, j);
            else
                system.A.insert(k, k - N) = cY; 

            if (j == N) // right neighbor
                system.b(k) += (-cX) * grid(i, j + 1);
            else
                system.A.insert(k, k + 1) = cX; 

            if (j == 1) // left neighbor
                system.b(k) += (-cX) * grid(i, j - 1);
            else
                system.A.insert(k, k - 1) = cX;
        }
    }

    // Compress the sparse matrix for efficient computation
    system.A.makeCompressed();
}


void AnisotropicDiffusion::print ()
{

}


void AnisotropicDiffusion::save_grid (std::string solver_name)
{   
    std::string filename = std::to_string(this->GRID_SIZE) + "_" + std::to_string(this->epsilon) + ".dat";
    if (!solver_name.empty())
        filename = solver_name + "_" + filename;

    std::filesystem::path full_path = std::filesystem::current_path() / filename;
    std::filesystem::create_directories(full_path.parent_path());

    std::ofstream file(full_path);

    for (int i = 0; i < this->grid.rows(); i++)
    {
        for (int j = 0; j < this->grid.rows(); j++)
        {   
            if (!solver_name.empty())
                file << this->grid(i, j) << ' ';
            else
                file << this->analytical_solution(i, j) << ' ';
        }
        file << '\n';
    }
    std::cout << "[OK] Solution successfully saved to " << full_path << ".\n";
}