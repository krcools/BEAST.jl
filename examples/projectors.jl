using BEAST, CompScienceMeshes, LinearAlgebra
using Plots
# using JLD2

setminus(A,B) = submesh(!in(B), A)

radius = 1.0

nearstrat = BEAST.DoubleNumWiltonSauterQStrat(6, 7, 6, 7, 7, 7, 7, 7)
farstrat  = BEAST.DoubleNumQStrat(1,2)

dmat(op,tfs,bfs) = BEAST.assemble(op,tfs,bfs; quadstrat=nearstrat)
# hmat(op,tfs,bfs) = AdaptiveCrossApproximation.h1compress(op,tfs,bfs; nearstrat=nearstrat,farstrat=farstrat)
mat = dmat

h = [0.1, 0.05, 0.025, 0.0125]
κ = [1.0, 10.0]

h = 0.3
κ = 0.00001
γ = im*κ

# function runsim(;h, κ)

SL = Maxwell3D.singlelayer(wavenumber=κ)
WS = Maxwell3D.weaklysingular(wavenumber=κ)
HS = Maxwell3D.hypersingular(wavenumber=κ)
N = NCross()

E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
e = (n × E) × n;

Γ = meshsphere(;radius, h)
∂Γ = boundary(Γ)

edges = setminus(skeleton(Γ,1), ∂Γ)
verts = setminus(skeleton(Γ,0), skeleton(∂Γ,0))

Σ = Matrix(connectivity(Γ, edges, sign))
Λ = Matrix(connectivity(Γ, edges, sign))

PΣ = Σ * pinv(Σ'*Σ) * Σ'
PΛH = I - PΣ

ℙΛ = Λ * pinv(Λ'*Λ) * Λ'
ℙHΣ = I - ℙΛ

MR = γ * PΣ + PΛH
ML = PΣ + 1/γ * PΛH

MRΣ = γ * PΣ
MRΛH = PΛH
MLΣ = PΣ
MLΛH = 1/γ * PΛH

𝕄R = γ * ℙΛ + ℙHΣ
𝕄L = ℙΛ + 1/γ * ℙHΣ

X = raviartthomas(Γ)
Y = buffachristiansen(Γ)

@hilbertspace p
@hilbertspace q

SLxx = assemble(@discretise(SL[p,q], p∈X, q∈X), materialize=mat)
WSxx = assemble(@discretise(WS[p,q], p∈X, q∈X), materialize=mat)

ex = BEAST.assemble(@discretise e[p] p∈X)

sys0 = SLxx
sys1 = MLΣ * SLxx * MRΣ + MLΛH * WSxx * MRΣ + MLΣ * WSxx * MRΛH + MLΛH * WSxx * MRΛH

rhs0 = ex
rhs1 = ML * ex

u0, ch0 = solve(BEAST.GMRESSolver(sys0, tol=2e-5, restart=250), rhs0)
v1, ch1 = solve(BEAST.GMRESSolver(sys1, tol=2e-5, restart=250), rhs1)

u1 = MR * v1
# u2 = MR * v2

# error()
# @show ch1.iters
# @show ch2.iters

Φ, Θ = [0.0], range(0,stop=π,length=40)
pts = [point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for ϕ in Φ for θ in Θ]

near0 = potential(MWFarField3D(wavenumber=κ), pts, u0, X)
near1 = potential(MWFarField3D(wavenumber=κ), pts, u1, X)
    # near2 = potential(MWFarField3D(wavenumber=κ), pts, u2, X)

    # u1, ch1.iters, u2, ch2.iters, near1, near2, X
# end

plot();
plot!(Θ, norm.(near0));
scatter!(Θ, norm.(near1))
# scatter!(Θ, norm.(near2))

error()

using LinearAlgebra
using Plots
plotly()
w0 = eigvals(Matrix(Sxx))
w1 = eigvals(Matrix(sys))
w2 = eigvals(Matrix(P * sys))
plot(exp.(2pi*im*range(0,1,length=200)))
scatter!(w0)
scatter!(w1)
scatter!(w2)


function makesim(d::Dict)
    @unpack h, κ = d
    u1, ch1, u2, ch2, near1, near2, X = runsim(;h, κ)
    fulld = merge(d, Dict(
        "u1" => u1,
        "u2" => u2,
        "ch1" => ch1,
        "ch2" => ch2,
        "near1" => near1,
        "near2" => near2
    ))
end

method = splitext(basename(@__FILE__))[1]


params = @strdict h κ
dicts = dict_list(params)
for (i,d) in enumerate(dicts)
    @show d
    f = makesim(d)
    @tagsave(datadir("simulations", method, savename(d,"jld2")), f)
end

#' Visualise the spectrum

# mSxx = BEAST.convert_to_dense(Sxx)
# mSyy = BEAST.convert_to_dense(Syy)

# Z = mSxx
# W = iN' * mSyy * iN * mSxx;

# wZ = eigvals(Matrix(Z))
# wW = eigvals(Matrix(W))

# plot(exp.(im*range(0,2pi,length=200)))
# scatter!(wZ)
# scatter!(wW)


# Study the various kernels
# HS = Maxwell3D.singlelayer(gamma=0.0, alpha=0.0, beta=1.0)
# Id = BEAST.Identity()

# Z12 = BEAST.lagrangecxd0(G12)
# Z23 = BEAST.lagrangecxd0(Ĝ23)
# Z = Z12 × Z23

# W12 = BEAST.duallagrangecxd0(G12)
# W23 = BEAST.duallagrangecxd0(Ĝ23)
# W = W12 × W23

# DX = assemble(Id, Z, divergence(X))
# HX = assemble(HS, X, X)

# DY = assemble(Id, W, divergence(Y))
# HY = assemble(HS, Y, Y)

# Nx = BEAST.NCross()
# NYX = assemble(Nx, Y, X)

# Q = HY * iN * HX

using AlgebraicMultigrid
A = poisson(100)
b = rand(100);
solve(A, b, RugeStubenAMG(), maxiter=1, abstol=1e-6)
