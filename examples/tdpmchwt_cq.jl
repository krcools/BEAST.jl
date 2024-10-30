using BEAST, LinearAlgebra, CompScienceMeshes, ConvolutionOperators, Printf

# Physical parameters
ϵ, μ = 1.0, 1.0
ϵ′, μ′ = 1.0, 1.0
η, η′ = √(μ/ϵ), √(μ′/ϵ′)

T = BEAST.LaplaceDomainOperator(s::ComplexF64 -> MWSingleLayer3D(s*√(ϵ*μ), -s*√(ϵ*μ), ComplexF64(1/√(ϵ*μ))/s))
K = BEAST.LaplaceDomainOperator(s::ComplexF64 -> MWDoubleLayer3D(s*√(ϵ*μ)))
T′ = BEAST.LaplaceDomainOperator(s::ComplexF64 -> MWSingleLayer3D(s*√(ϵ′*μ′), -s*√(ϵ′*μ′), ComplexF64(1/√(ϵ′*μ′))/s))
K′ = BEAST.LaplaceDomainOperator(s::ComplexF64 -> MWDoubleLayer3D(s*√(ϵ′*μ′)))

Γ = meshsphere(1.0, 0.3)
X = raviartthomas(Γ)
Y = buffachristiansen(Γ)

@hilbertspace k l
@hilbertspace j m

Nt, Δt = 200, 10.0
(A, b, c) = butcher_tableau_radau_2stages()
CQ = StagedTimeStep(Δt, Nt, c, A, b, 30, 1.0001)

duration = 40 * Δt                                       
delay = 60 * Δt                                        
amplitude = 1.0
gaussian = creategaussian(duration, delay, amplitude)
fgaussian = fouriertransform(gaussian)
polarisation, direction = x̂, ẑ
E = planewave(polarisation, direction, gaussian, 1.0)
H = direction × E

BEAST.@defaultquadstrat (T, X⊗CQ, X⊗CQ) BEAST.DoubleNumWiltonSauterQStrat(6, 7, 6, 7, 7, 7, 7, 7)
BEAST.@defaultquadstrat (T′, X⊗CQ, X⊗CQ) BEAST.DoubleNumWiltonSauterQStrat(6, 7, 6, 7, 7, 7, 7, 7)
BEAST.@defaultquadstrat (K, X⊗CQ, X⊗CQ) BEAST.DoubleNumWiltonSauterQStrat(6, 7, 6, 7, 7, 7, 7, 7)
BEAST.@defaultquadstrat (K′, X⊗CQ, X⊗CQ) BEAST.DoubleNumWiltonSauterQStrat(6, 7, 6, 7, 7, 7, 7, 7)

pmchwt = @discretise(η*T[k,j] + η′*T′[k,j] + K[l,j] + K′[l,j] 
                        - K[k,m] - K′[k,m] + 1/η*T[l,m] + 1/η′*T′[l,m] == E[k] + H[l],
                        k∈X⊗CQ, l∈X⊗CQ, j∈X⊗CQ, m∈X⊗CQ)
je = solve(pmchwt)


# Exact solution
N = NCross()
Nyx = assemble(N, Y, X)
iNyx = inv(Matrix(Nyx))
Zr = zeros(size(iNyx))

δ = timebasisdelta(Δt, Nt)	                			                       
linform = @discretise(H[k] + E[l], k∈Y⊗δ, l∈Y⊗δ)
rhs = BEAST.assemble(linform.linform, linform.test_space_dict)
j_ext = [iNyx Zr; Zr iNyx] * rhs


# Plot the numerical solutions
using Plots
x = 1e-2 * Δt/3 * [1:1:Nt;]

plt = Plots.plot(
    size = (600, 400),
    grid = false,
    xscale = :identity, 
    # xlims = (0, 1620),
    # xticks = [400, 800, 1200, 1600],
    xtickfont = Plots.font(10, "Computer Modern"), 
    yscale = :log10, 
    ylims = (1e-15, 1e0), 
    yticks = [1e-20, 1e-15, 1e-10, 1e-5, 1e0, 1e5],
    ytickfont = Plots.font(10, "Computer Modern"),
    xlabel = "Time ({\\mu}s)",
    ylabel = "Current density intensity (A/m)",
    titlefont = Plots.font(10, "Computer Modern"),
    guidefont = Plots.font(11, "Computer Modern"),
    colorbar_titlefont = Plots.font(10, "Computer Modern"),
    legendfont = Plots.font(11, "Computer Modern"),
    legend = :topright,
    dpi = 300)
    
Plots.plot!(x, abs.(j_ext[1, :]), label="Exact solution", linecolor=1, lw=1.2)
Plots.plot!(x, abs.(je[1, :]), label="TD-PMCHWT", linecolor=2, lw=1.2)