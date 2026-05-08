# In this example, the preconditioning ability of the On-Surface Radiation Condition (OSRC) Magnetic-to-Electric (MtE) map is shown on the Electric Field Integral Equation (EFIE), on the sphere.
# # OSRC operators 
# The MtE map ``\text{V}^{+}`` is a map of the trace of the magnetic field ``n \times h`` to the trace of the electric field ``n \times e`` on a scattering surface.
# The approximation of the local MtE surface operator as an On-Surface Radiation Condition (OSRC) operator on an arbitrary surface ``\Gamma`` is given by:
# ```math
# \begin{equation}
#     \text{V}^{+} \left( n \times h \right) + n \times e = 0, \quad \text{on} \ \Gamma,
# \end{equation}
# ```
# In the implementation a highly accurate pseudo-differential operator ``\text{V}_{\epsilon}^{+}`` is used to represent the Magnetic-to-Electric map.
# ```math
# \begin{aligned}
#     \text{V}_{\epsilon}^{+}  &:= - \pmb{\Theta}^{-1} \Lambda_{2, \epsilon}^{-1} \Lambda_{1, \epsilon}, \\
#     \pmb{\Theta} \left( \cdot \right) &:= \left( n \times \cdot \right), \\
#     \Lambda_{2, \epsilon} &:= I - \textbf{curl}_{\Gamma} \frac{1}{\kappa_{\epsilon}^2} \textup{curl}_{\Gamma}, \\
#     \Lambda_{1, \epsilon} &:= \left(I + \mathcal{J} \right)^{1/2}, \\
#     \mathcal{J} &:=\frac{\Delta_{\Gamma}}{\kappa_{\epsilon}^2} = \textbf{Grad}_{\Gamma} \frac{1}{\kappa_{\epsilon}^2} \textup{Div}_{\Gamma} - \textbf{curl}_{\Gamma} \frac{1}{\kappa_{\epsilon}^2} \textup{curl}_{\Gamma}. 
# \end{aligned}
# ```
# More details on this OSRC operator and its implementation are given in
# [An OSRC Preconditioner for the EFIE (Betcke et al. (2021))](https://arxiv.org/abs/2111.10761)
# The Electric-to-Magnetic (EtM) OSRC operator is also implemented, as given in 
# [Approximate local magnetic-to-electric surface operators for time-harmonic Maxwell’s equations (C. Geuzaine et al. (2014))](https://www.researchgate.net/publication/261636307_Approximate_local_magnetic-to-electric_surface_operators_for_time-harmonic_Maxwell's_equations)
#
# # OSRC EFIE preconditioning example
# We begin with loading needed packages
using BEAST, CompScienceMeshes

using Makeitso
using LinearAlgebra

using Plots
using PlotlyDocumenter

# Now, assemble and solve the unpreconditioned EFIE

BLAS.set_num_threads(8)

@target geo (;h) -> begin
      (; Γ = CompScienceMeshes.meshsphere(radius=1.0, h=h))
end

@target formulation_EFIE (geo,; κ) -> begin
    X = raviartthomas(geo.Γ)
    T = Maxwell3D.singlelayer(wavenumber=κ)

    E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ);
    e = (n × E) × n;

    bx = assemble(e, X)
    A = assemble(T,X,X); 
    return (;bilforms=(;A), linforms=(;bx))
end

@target solution_EFIE (formulation_EFIE,; residual) -> begin
      (;bilforms, linforms) = formulation_EFIE
      (;A) = bilforms; (;bx) = linforms;
      iT = BEAST.GMRESSolver(A; restart=1_500, reltol=residual, maxiter=1_500)
      u, ch = BEAST.solve(iT, bx)
      return (;iters=ch.iters, u=u)
end

# Next, assemble the OSRC MtE operator on the primal grid and use it to precondition the EFIE

@target OSRC_MtE_op (geo,;κ, Np, curvature) -> begin
    MtE_OSRC_operator = BEAST.MtE_OSRC_op(κ, Np, pi/2, curvature)
    Nd = BEAST.nedelec(geo.Γ)
    MtE_map = assemble(MtE_OSRC_operator, Nd, Nd)
    return (;MtE=MtE_map)
end

@target solution_OSRC_EFIE (formulation_EFIE, OSRC_MtE_op; residual) -> begin
      (;bilforms, linforms) = formulation_EFIE
      (;A) = bilforms; (;bx) = linforms;
      P_OSRC = OSRC_MtE_op.MtE
      iT = BEAST.GMRESSolver(A; restart=1_500, reltol=residual, maxiter=1_500, left_preconditioner=P_OSRC)
      u, ch = BEAST.solve(iT, bx)
      return (;iters=ch.iters, u=u)
end

# Finally, assemble and use the Calderón preconditioner, for a comparison with OSRC preconditioning

@target calderon_preconditioner (geo,;κ) -> begin
    Γ = geo.Γ
    X = raviartthomas(Γ)
    Y = BEAST.buffachristiansen(Γ)

    T = Maxwell3D.singlelayer(wavenumber=κ)
    N = NCross()

    Tyy = assemble(T,Y,Y);
    Nxy = Matrix(assemble(N,X,Y));
    iNxy = inv(Nxy);
    P = iNxy' * Tyy * iNxy
    return (;P=P)
end

@target solution_cald_EFIE (formulation_EFIE, calderon_preconditioner; residual) -> begin
      (;bilforms, linforms) = formulation_EFIE
      (;A) = bilforms; (;bx) = linforms;
      P_calderon = calderon_preconditioner.P
      iT = BEAST.GMRESSolver(A; restart=1_500, reltol=residual, maxiter=1_500, left_preconditioner=P_calderon)
      u, ch = BEAST.solve(iT, bx)
      return (;iters=ch.iters, u=u)
end

# Set the simulation parameters and run for different mesh sizes ``h``

@target benchmark_OSRC (solution_EFIE, solution_OSRC_EFIE, solution_cald_EFIE, ;) -> begin
      return (;iters_EFIE = solution_EFIE.iters, iters_OSRC_EFIE = solution_OSRC_EFIE.iters, iters_cald_EFIE = solution_cald_EFIE.iters)
end
@sweep benchmark_OSRC_sweep (!benchmark_OSRC,; h=[], κ=[], residual=[], Np=[], curvature=[]) -> benchmark_OSRC

h_values = [0.3, 0.2, 0.15]
κ = 1.0*pi
residual = 1e-6
Np = 4
curvature = 1.0

df = make(benchmark_OSRC_sweep; h=h_values, κ=κ, residual=residual, Np=Np, curvature=curvature)

# Finally, the GMRES iterations for the different methods are plotted
using Plots 
plot(df.h, df.iters_EFIE, label="EFIE", l=2)
plot!(df.h, df.iters_OSRC_EFIE, label="OSRC-EFIE", l=2)
plot!(df.h, df.iters_cald_EFIE, label="Calderón-EFIE", l=2)
xlabel!("h")
ylabel!("GMRES iterations")