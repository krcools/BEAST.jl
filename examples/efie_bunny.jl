using CompScienceMeshes
using BEAST
using LinearAlgebra
using PlotlyJS

# Flat triangle mesh through the corners of a curvilinear mesh, with the same cell order
function flatmesh(m::CompScienceMeshes.CurvilinearMesh)
    corners = [c.indices[1:3] for c in cells(m)]
    used = sort(unique(reduce(vcat, corners)))
    newid = Dict(v => i for (i, v) in enumerate(used))
    faces = [CompScienceMeshes.SimplexGraph(index((newid[v] for v in c)...)) for c in corners]
    return Mesh(vertices(m)[used], faces)
end

κ, η = 15.0, 1.0
λ = 2π / κ
h = λ / 10

# Stanford bunny scaled to a largest dimension of 1, meshed with cubic triangles.
# Γf is the flat mesh through the same corner points, with the same cells.
Γ = CompScienceMeshes.meshbunny(h; order=3)
Γ.vertices .= [point(v[1], -v[3], v[2]) for v in Γ.vertices]
zmin = minimum(v[3] for v in Γ.vertices)
Γ.vertices .= [v - point(0, 0, zmin) for v in Γ.vertices]

Γf = flatmesh(Γ)

# Lowest order div-conforming (GWP order 0, i.e. RT) functions on the flat mesh, and
# their counterparts on the cubic mesh, reordered to share the numbering and orientation
# of the DoFs. (Xf keeps its natural order: it matches the order of the BC functions,
# which keeps the duality matrix below close to diagonal and easy for GMRES.)
X = BEAST.gwpdiv(Γ; order=2)
@show numfunctions(X)
# Xf = BEAST.gwpdiv(Γf; order=0)
# edgeof(fn) = sort([(s.cellid, s.refid) for s in fn])
# perm = getindex.(Ref(Dict(edgeof(fn) => i for (i, fn) in enumerate(X.fns))), edgeof.(Xf.fns))
# X = BEAST.GWPDivSpace(X.geo, X.fns[perm], X.pos[perm], X.degree)
# shapes(fn) = sort([(s.cellid, s.refid, s.coeff) for s in fn])
# @assert all(all(a[1:2] == b[1:2] && a[3] ≈ b[3] for (a, b) in zip(shapes(x), shapes(y)))
#     for (x, y) in zip(X.fns, Xf.fns))
@info "Mesh: $(numcells(Γ)) cubic triangles, $(numfunctions(X)) unknowns, h = λ/10 = $(round(h, digits=4))"

t = Maxwell3D.singlelayer(wavenumber=κ)
E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
e = (n × E) × n

# Quadrature suitable for curved elements (cf. BEAST's examples/efie_geop2.jl)
BEAST.@defaultquadstrat (t,X,X) BEAST.DoubleNumSauterQstrat(4,5,6,6,6,6)
BEAST.@defaultquadstrat (e,X) BEAST.SingleNumQStrat(12)

A = assemble(t, X, X; threading=:cellcoloring)
b = assemble(e, X)

# # Calderón preconditioner P = N⁻ᵀ T_YY N⁻¹ with Buffa-Christiansen functions Y.
# # BEAST cannot yet build BC functions on curvilinear meshes (that requires a curved
# # barycentric refinement), so P is built on the flat mesh through the same vertices.
# # P only has to be spectrally equivalent to the inverse of A; the accuracy of the
# # solution is determined by A and b on the cubic geometry.
# Y = buffachristiansen(Γf)
# Tyy = assemble(t, Y, Y)
# Nxy = assemble(NCross(), Xf, Y)
# NYX = BEAST.GMRESSolver(Nxy; restart=1_500, abstol=1e-10, reltol=1e-10, maxiter=1_500, verbose=false)
# NXY = BEAST.GMRESSolver(transpose(Nxy); restart=1_500, abstol=1e-10, reltol=1e-10, maxiter=1_500, verbose=false)
# P = NXY * Tyy * NYX
# tol = 1e-6
# solver = BEAST.GMRESSolver(A; left_preconditioner=P, restart=500, reltol=tol, maxiter=2_000, verbose=false)

solver = BEAST.lu(A)

u, ch = solve(solver, b)
# @info "Calderón preconditioned GMRES: $(ch.iters) iterations, converged = $(ch.isconverged)"

fcr, geo = facecurrents(u, X; order=5)
@show length(geo)
plt = PlotlyJS.plot(
    patch(geo, norm.(fcr), caxis=(0,2)),
    Layout(title="EFIE on cubic Stanford bunny, κ = $κ: |j|",
        template="plotly_white",
        scene=attr(aspectmode="data",
        # camera=attr(up=attr(x=0, y=1, z=0),
        # eye=attr(x=0.9, y=0.6, z=1.6),)
    )),
)

# outfile = joinpath(@__DIR__, "efie_bunny_cubic.png")
# savefig(plt, outfile; width=900, height=800)
# println("Plot written to $outfile")
