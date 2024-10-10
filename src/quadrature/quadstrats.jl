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


defaultquadstrat(op, tfs, bfs) = defaultquadstrat(op, refspace(tfs), refspace(bfs))
macro defaultquadstrat(dop, body)
    @assert dop.head == :tuple
    if length(dop.args) == 3
        op = dop.args[1]
        tfs = dop.args[2]
        bfs = dop.args[3]
        ex = quote
            function BEAST.defaultquadstrat(::typeof($op), ::typeof($tfs), ::typeof($bfs))
                $body
            end
        end
        return esc(ex)
    elseif length(dop.args) == 2
        lin = dop.args[1]
        tfs = dop.args[2]
        ex = quote
            function BEAST.defaultquadstrat(::typeof($lin), ::typeof($tfs))
                $body
            end
        end
        return esc(ex)
    end
    error("@defaultquadstrat expects a first argument of the for (op,tfs,bfs) or (linform,tfs)")
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

"""
    quaddata(operator, test_refspace, trial_refspace, test_elements, trial_elements)

Returns an object cashing data required for the computation of boundary element
interactions. It is up to the client programmer to decide what (if any) data is
cached. For double numberical quadrature, storing the integration points for
example can significantly speed up matrix assembly.

- `operator` is an integration kernel.
- `test_refspace` and `trial_refspace` are reference space objects. `quadata`
is typically overloaded on the type of these local spaces of shape functions.
(See the implementation in `maxwell.jl` for an example).
- `test_elements` and `trial_elements` are iterable collections of the geometric
elements on which the finite element space are defined. These are provided to
allow computation of the actual integrations points - as opposed to only their
coordinates.
"""
function quaddata end


"""
    quadrule(operator, test_refspace, trial_refspace, test_index, test_chart, trial_index, trial_chart, quad_data)

Based on the operator kernel and the test and trial elements, this function builds
an object whose type and data fields specify the quadrature rule that needs to be
used to accurately compute the interaction integrals. The `quad_data` object created
by `quaddata` is passed to allow reuse of any precomputed data such as quadrature
points and weights, geometric quantities, etc.

The type of the returned quadrature rule will help in deciding which method of
`momintegrals` to dispatch to.
"""
function quadrule end