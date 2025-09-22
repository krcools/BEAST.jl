
abstract type AbstractCompInt <: AbstractOperator end
abstract type Kernel{T} end

#Make these structures, those are not yet pulled trough the smartmath function
struct CompDoubleInt{T} <: AbstractCompInt
    tfunc
    op1
    kernel::Kernel{T}
    op2
    bfunc
end

struct CompSingleInt{T} <: AbstractCompInt
    tfunc
    op1
    kernel::Kernel{T}
    op2
    bfunc
end
scalartype(op::CompDoubleInt{T}) where {T} = T 
scalartype(op::CompSingleInt{T}) where {T} = T 


struct CompDoubleKern{T,U,V,W} <: IntegralOperator
    op1::U
    kernel::W
    op2::V
end
function CompDoubleKern(p1,k::Kernel{T},p2) where {T}
    CompDoubleKern{T,typeof(p1),typeof(p2),typeof(k)}(p1,k,p2)
end
struct CompSingleKern{T,U,V,W} <: LocalOperator
    op1::U
    kernel::W
    op2::V
end
function CompSingleKern(p1,k::Kernel{T},p2) where {T}
    CompSingleKern{T,typeof(p1),typeof(p2),typeof(k)}(p1,k,p2)
end
scalartype(op::CompDoubleKern{T}) where {T} = T 
scalartype(op::CompSingleKern{T}) where {T} = T 

integralop(a::CompDoubleInt) = CompDoubleKern(a.op1,a.kernel,a.op2)
integralop(a::CompSingleInt) = CompSingleKern(a.op1,a.kernel,a.op2)


function assemble!(op::AbstractCompInt, tfs::Space, bfs::Space,
    store, threading = Threading{:multi};
    quadstrat=defaultquadstrat)

    assemble!(integralop(op),op.tfunc(tfs),op.bfunc(bfs),store,threading;quadstrat)

end



####### kernel definition (also posibilities to do symbolic actions)

# struct NullKernel <: Kernel{Int} end
# struct NxDotNy{T} <: Kernel{T} 
#     mult::Kernel{T}
# end
# NxDotNy() = NxDotNy(IdentityKernel())
# NxDotNy(a::NxDotNy) = @error "dont do nx⋅ny*nx⋅ny (yet)"

struct HH3DGreen{T} <: Kernel{T}
    gamma::T
end
struct HH3DGradGreen{T} <: Kernel{T}
    gamma::T
end

function (op::HH3DGreen)(x,y)
    gamma = op.gamma

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1/R
    green = exp(-gamma*R)*(i4pi*iR)
    green
end
function (op::HH3DGreen)(x::Union{SVector,Vector},y)
    gamma = op.gamma

    r = x - cartesian(y)
    R = norm(r)
    iR = 1/R
    green = exp(-gamma*R)*(i4pi*iR)
    green
end

function (op::HH3DGradGreen)(x,y)
    gamma = op.gamma

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1/R
    green = exp(-gamma*R)*(iR*i4pi)
    gradgreen = -(gamma + iR) * green * (iR * r)

    return gradgreen
end
function (op::HH3DGradGreen)(x::Union{SVector,Vector},y)
    gamma = op.gamma

    r = x - cartesian(y)
    R = norm(r)
    iR = 1/R
    green = exp(-gamma*R)*(iR*i4pi)
    gradgreen = -(gamma + iR) * green * (iR * r)
    tt = gradgreen

    return tt
end

##### integrand evaluation
function integrand(op::CompSingleKern,kerneldata,x,f,g,disp)
    kernel = op.kernel(x,disp)
    op1 = op.op1
    op2 = op.op2

    op1(f.value,op2(kernel,g.value))

end

function (op::Integrand{<:CompDoubleKern})(x,y,f,g)
    kernel = op.operator.kernel(x,y)
    op1 = op.operator.op1
    op2 = op.operator.op2
    _integrands(f,g) do fi, gi
        op1(fi.value,op2(kernel,gi.value))
    end
end

defaultquadstrat(op::CompDoubleInt,tfs::Space,bfs::Space) = DoubleNumSauterQstrat(5,5,5,5,5,5)
defaultquadstrat(op::CompSingleInt,tfs::Space,bfs::Space) = SingleNumQStrat(6)
##### assembling the single integral 
function assemble!(biop::CompSingleKern, tfs::Space, bfs::Space, store,
    threading::Type{Threading{:multi}};
    quadstrat=defaultquadstrat(biop, tfs, bfs))

    return assemble_local_mixed!(biop, tfs, bfs, store; quadstrat)
end

function assemble_local_mixed!(biop::CompSingleKern, tfs::Space{T}, bfs::Space{T}, store;
    quadstrat=defaultquadstrat) where {T}
    tol = sqrt(eps(T))

    quadstrat = quadstrat(biop, tfs, bfs)

    trefs = refspace(tfs)
    brefs = refspace(bfs)

    tels, tad, tact_to_global = assemblydata(tfs)
    bels, bad, bact_to_global = assemblydata(bfs)

    tgeo = geometry(tfs)
    bgeo = geometry(bfs)

    tdom = domain(chart(tgeo, first(tgeo)))
    bdom = domain(chart(bgeo, first(bgeo)))

    num_trefs = numfunctions(trefs, tdom)
    num_brefs = numfunctions(brefs, bdom)

    qd = quaddata(biop, trefs, brefs, tels, bels, quadstrat)

    # store the bcells in an octree
    tree = elementstree(bels)

    print("dots out of 10: ")
    todo, done, pctg = length(tels), 0, 0
    for (p,tcell) in enumerate(tels)

        tc, ts = boundingbox(tcell.vertices)
        pred = (c,s) -> boxesoverlap(c,s,tc,ts)

        for box in boxes(tree, pred)
            for q in box
                bcell = bels[q]

                if overlap(tcell, bcell)

                    isct = intersection2(tcell, bcell)
                    disp = displacement(displacementchart(geometry(tfs),tact_to_global[p]),displacementchart(geometry(bfs),bact_to_global[q])) #number to multiply with test normal, if zero use the orientation of displacementvector
                    for (cellt,cellb) in isct
                        volume(cellb) < tol && continue
                        # @assert sign(dot(normal(chart(geometry(tfs),p)),normal(cellt))) == 1
                        # @assert sign(dot(normal(chart(geometry(bfs),q)),normal(cellb))) == 1

                        P = restrict(brefs, bcell, cellb)
                        Q = restrict(trefs, tcell, cellt)

                        qr = quadrule(biop, trefs, brefs, cellt,cellb, qd, quadstrat)
                        zlocal = cellinteractions(biop, trefs, brefs, cellt, qr, disp)
                        zlocal = Q * zlocal * P'

                        for i in 1 : num_trefs
                            for j in 1 : num_brefs
                                for (m,a) in tad[p,i]
                                    for (n,b) in bad[q,j]
                                        store(a * zlocal[i,j] * b, m, n)
                                    end # next basis function this cell supports
                                end # next test function this cell supports
                            end # next refshape on basis side
                        end # next refshape on test side

                    end # next cell in intersection
                end # if overlap
            end # next cell in the basis geometry
        end # next box in the octree

        done += 1
        new_pctg = round(Int, done / todo * 100)
        if new_pctg > pctg + 9
            print(".")
            pctg = new_pctg
        end
    end # next cell in the test geometry

    println("")
end
kernelvals(::CompSingleKern,x) = nothing
function cellinteractions(biop::CompSingleKern, trefs::U, brefs::V, cell, qr, disp) where {U<:RefSpace{T},V<:RefSpace{T}} where {T}

    num_tshs = length(qr[1][3])
    num_bshs = length(qr[1][4])

    zlocal = zeros(T, num_tshs, num_bshs)
    for q in qr

        w, mp, tvals, bvals = q[1], q[2], q[3], q[4]
        j = w * jacobian(mp)
        kernel = kernelvals(biop, mp)

        for m in 1 : num_tshs
            tval = tvals[m]

            for n in 1 : num_bshs
                bval = bvals[n]

                igd = integrand(biop, kernel, mp, tval, bval,disp)
                zlocal[m,n] += j * igd

            end
        end
    end

    return zlocal
end

function quadrule(op::CompSingleKern, ψ::RefSpace, ϕ::RefSpace, cellt,cellb, (qd,A), qs::SingleNumQStrat)
    for i in eachindex(qd)
        q = qd[i]
        w, pt = q[1], neighborhood(cellt,q[2])
        pb = neighborhood(cellb,carttobary(cellb,cartesian(pt)))
        A[i] = (w, pb, ψ(pt), ϕ(pb))
    end
    return A
end

