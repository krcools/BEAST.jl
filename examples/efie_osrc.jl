# # OSRC-preconditioned EFIE
# In this example, the preconditioning ability of the On-Surface Radiation Condition (OSRC) Magnetic-to-Electric (MtE) map is shown on the Electric Field Integral Equation (EFIE), on the cuboid.

# We begin with loading needed packages
using BEAST, CompScienceMeshes

using Makeitso
using PlotlyDocumenter

using LinearAlgebra

using Plots
BLAS.set_num_threads(32)

# First, define the cuboid geometry. The unpreconditioned EFIE can then be assembled and solved.

@target geo (;h) -> begin
      (; Γ = CompScienceMeshes.meshcuboid(1.0, 1.0, 1.0,  h))
end

@target formulation_EFIE (geo,; κ) -> begin
    X = raviartthomas(geo.Γ)
    T = Maxwell3D.singlelayer(wavenumber=κ)

    E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ);
    e = (n × E) × n;

    bx = assemble(e, X)
    Axx = assemble(T,X,X); 
    return (;Axx=Axx, bx=bx)
end

@target solution_EFIE (formulation_EFIE,) -> begin
    (;Axx, bx) = formulation_EFIE
    iT = BEAST.GMRESSolver(Axx; restart=1_500, maxiter=1_500)
    u, ch = BEAST.solve(iT, bx)
    return (;iters=ch.iters, u=u)
end;

# Next, assemble the OSRC MtE operator on the primal grid, using the original basis functions of the EFIE.
# We then use it as a preconditioner for the EFIE, during the solution proces.

@target OSRC_MtE (geo,;κ, Np, θ_p, curvature) -> begin
    MtE_OSRC_operator = BEAST.MteOsrc(κ, Np, θ_p, curvature)
    X = BEAST.raviartthomas(geo.Γ)
    MtE_map = assemble(MtE_OSRC_operator, X, X)
    return (;MtE=MtE_map)
end

@target solution_OSRC_EFIE (formulation_EFIE, OSRC_MtE;) -> begin
    (;Axx, bx) = formulation_EFIE
    MtE_map = OSRC_MtE.MtE
    PAxx = MtE_map*Axx
    Pbx = MtE_map*bx
    iT = BEAST.GMRESSolver(PAxx; restart=1_500, reltol=1e-6, maxiter=1_500)
    u, ch = BEAST.solve(iT, Pbx)
    return (;iters=ch.iters, u=u)
end;

# Set now the simulation parameters and run for different mesh sizes `h`.

@target benchmark_OSRC (solution_EFIE, solution_OSRC_EFIE, ;) -> begin
      return (;iters_EFIE = solution_EFIE.iters, iters_OSRC_EFIE = solution_OSRC_EFIE.iters)
end
@sweep benchmark_OSRC_sweep (!benchmark_OSRC,; h=[], κ=[]) -> benchmark_OSRC

κ = 1.0*pi
h_values = [0.4, 0.3, 0.2, 0.1]
Np = 10                                 # number of Padé terms
θ_p = pi/2                              # rotation angle of the poles
curvature = 1.0                         # curvature parameter

df = make(benchmark_OSRC_sweep; h=h_values, κ=κ, Np=Np, θ_p=θ_p,curvature=curvature);

# Finally, the GMRES iterations for the unpreconditioned and OSRC-preconditioned EFIE are plotted in function of `h`.
# The reduced GMRES iterations for the OSRC-preconditioned EFIE indicates the speed-up in solving time.
sort!(df, "h")
Plots.plot(df.h, df.iters_EFIE, label="EFIE", l=2)
Plots.plot!(df.h, df.iters_OSRC_EFIE, label="OSRC-EFIE", l=2)
Plots.xlabel!("h")
Plots.ylabel!("GMRES iterations")
