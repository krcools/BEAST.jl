struct QuasiLocalSingleLayerMW{R<:Real,F} <: QuasiLocalOperator
    gamma::R
    alpha::R
    beta::R
    range::R
    nominator::F
end

scalartype(op::QuasiLocalSingleLayerMW{T}) where {T} = T
defaultquadstrat(op::QuasiLocalSingleLayerMW,
    tfs::RTRefSpace, bfs::RTRefSpace) = DoubleNumSauterQstrat(4,5,5,5,4,3)

oprange(op::QuasiLocalSingleLayerMW) = op.range


struct ExponentialNominator{R}
    gamma::R
end
function (f::ExponentialNominator)(x,y)
    exp(-f.gamma * norm(x-y))
end

function QuasiLocalSingleLayerMW(gamma, range)
    f = ExponentialNominator(gamma)
    QuasiLocalSingleLayerMW(gamma, -gamma, -1/(gamma), range, f)
end

function QuasiLocalSingleLayerMW(gamma, range, f)
    QuasiLocalSingleLayerMW(gamma, -gamma, -1/(gamma), range, f)
end

function (igd::Integrand{<:QuasiLocalSingleLayerMW})(x,y,f,g)
    α = igd.operator.alpha
    β = igd.operator.beta
    γ = igd.operator.gamma
    F = igd.operator.nominator

    cx = cartesian(x)
    cy = cartesian(y)
    r = cx - cy
    R = norm(r)
    # iR = 1 / R
    # green = exp(-γ*R) / (4π*R)
    green = F(cx,cy) / (4π*R)

    αG = α * green
    βG = β * green

    _integrands(f,g) do fi,gj
        αG * dot(fi.value, gj.value) + βG * dot(fi.divergence, gj.divergence)
    end
end

@testitem "quasi-local singlelayer MW" begin
    using CompScienceMeshes, LinearAlgebra
    using SparseArrays

    fn = joinpath(dirname(pathof(BEAST)),"../test","assets","sphere45.in")
    Γ1 = readmesh(fn)
    Γ2 = CompScienceMeshes.translate(Γ1, point(0,0,3.0))
    # Γ = meshsphere(radius=1.0, h=0.075)

    δ = 0.025
    γ = 1/δ
    t1 = BEAST.QuasiLocalSingleLayerMW(γ, 12*δ)
    t2 = BEAST.MWSingleLayer3D(γ)
    X1 = raviartthomas(Γ1)
    X2 = raviartthomas(Γ2)

    @time Z11 = assemble(t1, X1, X1)
    @time W11 = assemble(t2, X1, X1)
    @time Z12 = assemble(t1, X1, X2)
    @time W12 = assemble(t2, X1, X2)
    # @show norm(Z11-W11)
    # @show norm(W11)
    # @show norm(Z12-W12)
    # @show norm(W12)
    # @show norm(Z12)
    @show length(nonzeros(Z11)) / length(Z11)

    @test norm(Z11-W11) < 1e-5
    @test norm(Z12-W12) < 1e-19
    @test norm(Z12) == 0
    @test length(nonzeros(Z11)) / length(Z11) < 0.86
end