struct GWPDivSpace{T,M,P} <: Space{T}
    geo::M
    fns::Vector{Vector{Shape{T}}}
    pos::P
    degree::Int
end

function refspace(s::GWPDivSpace{T}) where {T} GWPDivRefSpace{T,s.degree}() end
function subset(s::S,I) where {S<:GWPDivSpace} S(s.geo, s.fns[I], s.pos[I], s.degree) end

function gwpdiv(mesh, edges=nothing; order)

    T = coordtype(mesh)
    S = Shape{T}

    space = BEAST.gwpcurl(mesh, edges; order)
    fns = Vector{Vector{S}}(undef, length(space.fns))
    for (i,fn) in enumerate(space.fns)
        fns[i] = [S(s.cellid, s.refid, -s.coeff) for s in fn]
    end

    return GWPDivSpace(space.geo, fns, space.pos, order)
end

@testitem "GWPcurl global: numfunctions" begin
    using CompScienceMeshes

    h = 0.5
    mesh = meshrectangle(1.0, 1.0, 0.5)
    edges = setminus(skeleton(mesh,1), boundary(mesh))

    order = 2
    gwp = BEAST.gwpdiv(mesh, edges; order=order)

    ne = order+1
    nf = order * (order+1)
    Nt = length(edges)*ne + length(mesh)*nf
    @test numfunctions(gwp) == Nt
end


function divergence(X::GWPDivSpace, geo, fns)
    ch = chart(geo, first(geo))
    dim = dimension(ch)
    dom = domain(ch)
    # dim = dimension(dom)

    P = X.degree
    Cont = -1
    NF = binomial(dim+P, dim)

    LagrangeBasis{P,Cont,NF}(geo, fns, deepcopy(positions(X)))
end


@testitem "divergence - global" begin
    using CompScienceMeshes
    const CSM = CompScienceMeshes

    T = Float64

    m = CSM.meshrectangle(1.0, 1.0, 0.5, 3)
    X = BEAST.gwpdiv(m; order=4)
    divX = BEAST.divergence(X)

    x = BEAST.refspace(X)
    divx = BEAST.refspace(divX)

    err = zero(T)
    for i in eachindex(X.fns)
        fn = X.fns[i]
        for j in eachindex(fn)
            cellid = X.fns[i][j].cellid
            ch = chart(m, cellid)
            
            u = (0.2341, 0.4312)
            p = neighborhood(ch, u)

            r1 = zero(T)
            ϕp = x(p)
            for sh in X.fns[i]
                sh.cellid == cellid || continue
                r1 += sh.coeff * ϕp[sh.refid].divergence
            end

            r2 = zero(T)
            ϕp = divx(p)
            for sh in divX.fns[i]
                sh.cellid == cellid || continue
                r2 += sh.coeff * ϕp[sh.refid].value
            end
            global err = max(err, abs(r1-r2))
        end
    end

    @test err < sqrt(eps(T))
end