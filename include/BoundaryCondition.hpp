#ifndef BOUNDARY_CONDITION_HPP
#define BOUNDARY_CONDITION_HPP

#include <ginac/ginac.h>
#include <optional>
#include <variant>
#include <iostream>

enum class Side { Top, Bottom, Left, Right };;

inline constexpr std::string_view to_string(Side side) {
    switch (side) {
        case Side::Top:    return "Top";
        case Side::Bottom: return "Bottom";
        case Side::Left:   return "Left";
        case Side::Right:  return "Right";
    }
    return "Unknown side";
}

enum class Type { Dirichlet, Neumann, Robin };

inline constexpr std::string_view to_string(Type type) {
    switch (type) {
        case Type::Dirichlet: return "Dirichlet";
        case Type::Neumann:   return "Neumann";
        case Type::Robin:     return "Robin";
    }
    return "Unknown type";
}

using BoundaryType = std::variant<GiNaC::ex, double>;

struct BoundaryCondition
{
    std::optional<BoundaryType> expr = std::nullopt;
    Side side;
    Type type;

    bool is_set()
    {
        return expr.has_value();
    }

    bool is_expression()
    {
        return expr.has_value() && std::holds_alternative<GiNaC::ex>(expr.value());
    }

    bool is_constant()
    {
        return expr.has_value() && std::holds_alternative<double>(expr.value());   
    }

    double evaluate (std::pair<int, int> coordinates, 
                     std::vector<GiNaC::symbol> vars, 
                     double h)
    {
        double val = 0.0;
        int x = coordinates.first, y = coordinates.second;

        if (this->is_set())
        {
            if (this->is_constant())
            {   
                val = std::get<double>(*expr);
                // printf("i: %d j: %d | x: %f y: %f val: %f\n", i, j, j*h, i*h, val);
                // printf("at point (x', y') -> (%f, %f)\n", x * h, y * h);
                // std::cout << val << " " << to_string(this->side) << "\n\n";
                return val;
            }
            else if (this->is_expression())
            {   
                GiNaC::exmap m; 
                m[vars[0]]                   = x * h; m[vars[1]] = y * h;
                GiNaC::ex b_expr             = std::get<GiNaC::ex>(*expr);
                GiNaC::ex evaluated_boundary = b_expr.subs(m).evalf();
                val                          = GiNaC::ex_to<GiNaC::numeric>(evaluated_boundary).to_double();
                return val;
            }
        }
        return val;
    }

    void print ()
    {
        if (!is_set())
        {
            std::cout << "Side " << to_string(side) << " was not initialised\n.";
            return;
        }
        else
        {
            std::cout << "Side: " << to_string(side) << '\n';
            if (is_constant())
            {
                std::cout << "Value: " << std::get<double>(*expr) << '\n';
            }
            else if (is_expression())
            {
                std::cout << "Expression: " << std::get<GiNaC::ex>(*expr) << '\n';
            }
        }
    }
};


struct BoundaryConditions
{
    BoundaryCondition top    {std::nullopt, Side::Top,    Type::Dirichlet};
    BoundaryCondition bottom {std::nullopt, Side::Bottom, Type::Dirichlet};
    BoundaryCondition left   {std::nullopt, Side::Left,   Type::Dirichlet}; 
    BoundaryCondition right  {std::nullopt, Side::Right,  Type::Dirichlet};
};

#endif // BOUNDARY_CONDITION_HPP