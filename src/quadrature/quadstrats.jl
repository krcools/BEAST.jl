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

struct DoubleNumSauterQstrat{R,S}
    outer_rule::R
    inner_rule::R
    sauter_schwab_common_tetr::S
    sauter_schwab_common_face::S
    sauter_schwab_common_edge::S
    sauter_schwab_common_vert::S
end

struct DoubleNumQStrat{R}
    outer_rule::R
    inner_rule::R
end

struct SauterSchwab3DQStrat{R,S}
    outer_rule::R
    inner_rule::R
    sauter_schwab_1D::S
    sauter_schwab_2D::S
    sauter_schwab_3D::S
    sauter_schwab_4D::S
end

struct OuterNumInnerAnalyticQStrat{R}
    outer_rule::R
end


function defaultquadstrat(op, tfs, bfs) 
    @warn "made use of general defaultquadstrat routine which cals defaultquadstrat on refspaces"
    defaultquadstrat(op, refspace(tfs), refspace(bfs))
end
macro defaultquadstrat(dop, body)
    @assert dop.head == :tuple
    @assert length(dop.args) == 3
    op = dop.args[1]
    tfs = dop.args[2]
    bfs = dop.args[3]
    ex = quote
        function BEAST.defaultquadstrat(::typeof($op), ::typeof($tfs), ::typeof($bfs))
            $body
        end
    end
    return esc(ex)
end

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