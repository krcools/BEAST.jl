using InteractiveUtils

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

function quadinfo(op, tfs, bfs; quadstrat=defaultquadstrat(op, tfs, bfs))

    tels, tad = assemblydata(tfs)
    bels, bad = assemblydata(bfs)

    tref = refspace(tfs)
    bref = refspace(bfs)

    i, τ = 1, first(tels)
    j, σ = 1, first(bels)

    @show quadstrat
    println(@which BEAST.quaddata(op,tref,bref,tels,bels,quadstrat))

    qd = quaddata(op,tref,bref,tels,bels,quadstrat)
    println(@which quadrule(op,tref,bref,i,τ,j,σ,qd,quadstrat))

    nothing
end