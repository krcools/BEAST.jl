struct DoubleNumWiltonSauterQStrat{R,S}
    outer_rule_far::R
    inner_rule_far::R
    outer_rule_near::R
    inner_rule_near::R
    sauter_schwab_common_tetr::S
    sauter_schwab_common_face::S
    sauter_schwab_common_edge::S
    sauter_schwab_common_vert::S
end

struct DoubleNumQStrat{R}
    outer_rule::R
    inner_rule::R
end


# defaultquadstrat(op, tfs, bfs) = DoubleNumWiltonSauterQStrat(2,3,6,7,6,5,4,3)
defaultquadstrat(op, tfs, bfs) = defaultquadstrat(op, refspace(tfs), refspace(bfs))
# defaultquadstrat(op, ::RefSpace, ::RefSpace) = error("No default quadstrat set.")

struct SingleNumQStrat
    quad_rule::Int
end