import Plotly
function showfn(space,i)
    geo = geometry(space)
    T = coordtype(geo)
    X = Dict{Int,T}()
    Y = Dict{Int,T}()
    Z = Dict{Int,T}()
    U = Dict{Int,T}()
    V = Dict{Int,T}()
    W = Dict{Int,T}()
    for (s,sh) in enumerate(space.fns[i])
        chrt = chart(geo, sh.cellid)
        #nbd = CompScienceMeshes.center(chrt)
        len = length(Iterators.product((0.0:0.1:1.0), (0.0:0.1:1.0)))
        for (m,n) in enumerate(Iterators.product((0.0:0.1:1.0), (0.0:0.1:1.0)))
            n[1]+n[2] > 1 && continue
            nbd = neighborhood(chrt,[n[1],n[2]])
            vals = refspace(space)(nbd)
            x,y,z = cartesian(nbd)
            #@show vals[sh.refid].value
            u,v,w = vals[sh.refid].value
            # @show x, y, z
            # @show u, v, w
            X[sh.cellid*len+m] = x
            Y[sh.cellid*len+m] = y
            Z[sh.cellid*len+m] = z
            U[sh.cellid*len+m] = get(U,sh.cellid*len+m,zero(T)) + sh.coeff * u
            V[sh.cellid*len+m] = get(V,sh.cellid*len+m,zero(T)) + sh.coeff * v
            W[sh.cellid*len+m] = get(W,sh.cellid*len+m,zero(T)) + sh.coeff * w
        end
    end
    X = collect(values(X))
    Y = collect(values(Y))
    Z = collect(values(Z))
    U = collect(values(U))
    V = collect(values(V))
    W = collect(values(W))
    Plotly.cone(x=X,y=Y,z=Z,u=U,v=V,w=W)
end