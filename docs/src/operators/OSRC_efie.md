```@meta
EditURL = "../../../examples/OSRC_efie.jl"
```

In this example, the preconditioning ability of the On-Surface Radiation Condition (OSRC) Magnetic-to-Electric (MtE) map is shown on the Electric Field Integral Equation (EFIE), on the sphere.
# OSRC operators
The MtE map ``\text{V}^{+}`` is a map of the trace of the magnetic field ``n \times h`` to the trace of the electric field ``n \times e`` on a scattering surface.
The approximation of the local MtE surface operator as an On-Surface Radiation Condition (OSRC) operator on an arbitrary surface ``\Gamma`` is given by:
```math
\begin{equation}
    \text{V}^{+} \left( n \times h \right) + n \times e = 0, \quad \text{on} \ \Gamma,
\end{equation}
```
In the implementation a highly accurate pseudo-differential operator ``\text{V}_{\epsilon}^{+}`` is used to represent the Magnetic-to-Electric map.
```math
\begin{aligned}
    \text{V}_{\epsilon}^{+}  &:= - \pmb{\Theta}^{-1} \Lambda_{2, \epsilon}^{-1} \Lambda_{1, \epsilon}, \\
    \pmb{\Theta} \left( \cdot \right) &:= \left( n \times \cdot \right), \\
    \Lambda_{2, \epsilon} &:= I - \textbf{curl}_{\Gamma} \frac{1}{\kappa_{\epsilon}^2} \textup{curl}_{\Gamma}, \\
    \Lambda_{1, \epsilon} &:= \left(I + \mathcal{J} \right)^{1/2}, \\
    \mathcal{J} &:=\frac{\Delta_{\Gamma}}{\kappa_{\epsilon}^2} = \textbf{Grad}_{\Gamma} \frac{1}{\kappa_{\epsilon}^2} \textup{Div}_{\Gamma} - \textbf{curl}_{\Gamma} \frac{1}{\kappa_{\epsilon}^2} \textup{curl}_{\Gamma}.
\end{aligned}
```
More details on this OSRC operator and its implementation are given in
[An OSRC Preconditioner for the EFIE (Betcke et al. (2021))](https://arxiv.org/abs/2111.10761)
The Electric-to-Magnetic (EtM) OSRC operator is also implemented, as given in
[Approximate local magnetic-to-electric surface operators for time-harmonic Maxwell’s equations (C. Geuzaine et al. (2014))](https://www.researchgate.net/publication/261636307_Approximate_local_magnetic-to-electric_surface_operators_for_time-harmonic_Maxwell's_equations)

# OSRC EFIE preconditioning example
We begin with loading needed packages

````julia
using BEAST, CompScienceMeshes

using Makeitso
using LinearAlgebra

using Plots
using PlotlyDocumenter
````

Now, assemble and solve the unpreconditioned EFIE

````julia
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
````

````
0x4d338de5db5a428f
````

Next, assemble the OSRC MtE operator on the primal grid and use it to precondition the EFIE

````julia
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
````

````
0x48b1ab5f44e88685
````

Finally, assemble and use the Calderón preconditioner, for a comparison with OSRC preconditioning

````julia
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
````

````
0x3b488b9b8fce1697
````

Set the simulation parameters and run for different mesh sizes ``h``

````julia
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
````

```@raw html
<div><div style = "float: left;"><span>3×8 DataFrame</span></div><div style = "clear: both;"></div></div><div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">iters_EFIE</th><th style = "text-align: left;">iters_cald_EFIE</th><th style = "text-align: left;">h</th><th style = "text-align: left;">residual</th><th style = "text-align: left;">curvature</th><th style = "text-align: left;">κ</th><th style = "text-align: left;">Np</th><th style = "text-align: left;">iters_OSRC_EFIE</th></tr><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;"></th><th title = "Union{Missing, Int64}" style = "text-align: left;">Int64?</th><th title = "Union{Missing, Int64}" style = "text-align: left;">Int64?</th><th title = "Union{Missing, Float64}" style = "text-align: left;">Float64?</th><th title = "Union{Missing, Float64}" style = "text-align: left;">Float64?</th><th title = "Union{Missing, Float64}" style = "text-align: left;">Float64?</th><th title = "Union{Missing, Float64}" style = "text-align: left;">Float64?</th><th title = "Union{Missing, Int64}" style = "text-align: left;">Int64?</th><th title = "Union{Missing, Int64}" style = "text-align: left;">Int64?</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: right;">156</td><td style = "text-align: right;">9</td><td style = "text-align: right;">0.15</td><td style = "text-align: right;">1.0e-6</td><td style = "text-align: right;">1.0</td><td style = "text-align: right;">3.14159</td><td style = "text-align: right;">4</td><td style = "text-align: right;">48</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: right;">129</td><td style = "text-align: right;">9</td><td style = "text-align: right;">0.2</td><td style = "text-align: right;">1.0e-6</td><td style = "text-align: right;">1.0</td><td style = "text-align: right;">3.14159</td><td style = "text-align: right;">4</td><td style = "text-align: right;">47</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">3</td><td style = "text-align: right;">106</td><td style = "text-align: right;">10</td><td style = "text-align: right;">0.3</td><td style = "text-align: right;">1.0e-6</td><td style = "text-align: right;">1.0</td><td style = "text-align: right;">3.14159</td><td style = "text-align: right;">4</td><td style = "text-align: right;">45</td></tr></tbody></table></div>
```

Finally, the GMRES iterations for the different methods are plotted

````julia
using Plots
plot(df.h, df.iters_EFIE, label="EFIE", l=2)
plot!(df.h, df.iters_OSRC_EFIE, label="OSRC-EFIE", l=2)
plot!(df.h, df.iters_cald_EFIE, label="Calderón-EFIE", l=2)
xlabel!("h")
ylabel!("GMRES iterations")
````

```@raw html
<img src="data:image/png;base64,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" />
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

