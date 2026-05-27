# Nested Quadrature Strategies

Some quadrature strategies are only useble for a subest of all simplex-simplex interactions. If two simplicces do fullfill the right conditions this quadrature strategie returns a specified quadrature rule. If not, it calls the nested quadrature strategy on the simplex-simplex pair. This nested quadrature strategy is stored in the nested_strat field that a nested quadrature strategy must have.
An example of a nested quadstrat is given by
```julia
struct SauterQStrat{NestedStrat,S} <: NestedQuadStrat
    nested_strat::NestedStrat
    sauter_schwab_common_tetr::S
    sauter_schwab_common_face::S
    sauter_schwab_common_edge::S
    sauter_schwab_common_vert::S
end
```