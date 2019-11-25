using CompScienceMeshes
using BEAST

dir = dirname(pathof(BEAST))
Γ = readmesh(joinpath(dir,"../examples/sphere2.in"))

X = raviartthomas(Γ)
Y = buffachristiansen(Γ)

bels, bad = BEAST.assemblydata(X)
tels, tad = BEAST.assemblydata(Y)

Y = BEAST.RTBasis(Y.geo.mesh, Y.fns, Y.pos)

Δt, Nt = 0.6, 200
T2 = timebasisshiftedlagrange(Δt, Nt, 2)
T3 = timebasisshiftedlagrange(Δt, Nt, 3)
δ = timebasisdelta(Δt, Nt)

DL = TDMaxwell3D.doublelayer(speedoflight=1.0)
SL = TDMaxwell3D.singlelayer(speedoflight=1.0, numdiffs=1)

# BEAST.allocatestorage(K,W,V,Val{:bandedstorage},BEAST.LongDelays{:ignore})

function as(op, W, V)

    X = spatialbasis(V)
    Y = spatialbasis(W)

    aux = BEAST.EmptyRP(1.0)

    M = numfunctions(Y)
    N = numfunctions(X)

    K0 = fill(typemax(Int), M, N)
    K1 = zeros(Int, M, N)

    function store(v,m,n,k)
        K0[m,n] = min(K0[m,n],k)
        K1[m,n] = max(K1[m,n],k)
    end

    BEAST.assemble_chunk!(aux, W, V, store)

end
