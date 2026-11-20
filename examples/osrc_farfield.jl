# # OSRC operators
# The On-Surface Radiation Condition (OSRC) operators are approximations of the mappings between the tangential trace of the magnetic field ``n \times h`` and the tangential trace of the electric field ``n \times e`` on a scattering surface ``\Gamma``.
# In particular, the OSRC Magnetic-to-Electric (MtE) map ``V_{\epsilon}`` approximates the relation
# ```math
# \begin{aligned}
#     V_{\epsilon} \left( n \times h \right) + n \times e = 0, \quad \text{on} \ \Gamma.
# \end{aligned}
# ```
# The OSRC MtE map is based on the pseudodifferential expression for the exact MtE map on the plane.
# This pseudodifferential operator provides a highly accurate approximation on general surfaces, in the high-frequency regime:
# ```math
# \begin{aligned}
#     V_{\epsilon}  &:= - \pmb{\Theta}^{-1} \Lambda_{2, \epsilon}^{-1} \Lambda_{1, \epsilon}, \\
#  \end{aligned}
# ```
# where
# ```math
# \begin{aligned}
#     \pmb{\Theta} \left( \cdot \right) &:= \left( n \times \cdot \right), \\
#     \Lambda_{2, \epsilon} &:= I - \textbf{curl}_{\Gamma} \frac{1}{\kappa_{\epsilon}^2} \textup{curl}_{\Gamma}, \\
#     \Lambda_{1, \epsilon} &:= \left(I + \mathcal{J} \right)^{1/2}, \\
#     \mathcal{J} &:=\frac{\Delta_{\Gamma}}{\kappa_{\epsilon}^2} = \textbf{grad}_{\Gamma} \frac{1}{\kappa_{\epsilon}^2} \textup{div}_{\Gamma} - \textbf{curl}_{\Gamma} \frac{1}{\kappa_{\epsilon}^2} \textup{curl}_{\Gamma}. 
# \end{aligned}
# ```
# Here, the complex wavenumber ``\kappa_{\epsilon}`` is defined as
# ```math
# \begin{aligned}
#     \kappa_{\epsilon} := \kappa - i \epsilon, \quad \epsilon = 0.39 \kappa^{1/3} (\mathcal{H})^{2/3},
# \end{aligned}
# ```
# where ``\mathcal{H}`` denotes the mean curvature of the surface.
# The square root in ``\Lambda_{1, \epsilon}`` is approximated using a rational rotated branch-cut Padé series:
# ```math
# \begin{aligned}
#   \Lambda_{1, \epsilon} &\approx \left( R_{0} I - \sum_{j=1}^{N_{p}} \frac{A_j}{B_j} \Pi_{j}^{-1} \right), \quad  \Pi_{j} := I + B_j \mathcal{J},\\
# \end{aligned}
# ```
# where ``A_j, B_j`` and ``R_0`` are the parameters of the Padé series.
# The Padé series is determined by two free parameters: the number of Padé terms ``N_p`` and the rotation angle of the poles ``\theta_p``.
# Further details on this OSRC operator and its discretization are provided in [An OSRC Preconditioner for the EFIE (I. Fierro-Piccardo et al., 2021)](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=10022032). 
# The Electric-to-Magnetic (EtM) OSRC operator ``\Lambda_{\epsilon}`` is defined similarly, see [Approximate local magnetic-to-electric surface operators for time-harmonic Maxwell’s equations (C. Geuzaine et al., 2014)](https://www.researchgate.net/publication/261636307_Approximate_local_magnetic-to-electric_surface_operators_for_time-harmonic_Maxwell's_equations)
# ```math
# \begin{aligned}
#     \Lambda_{\epsilon}  &:= - \Lambda_{1, \epsilon}^{-1} \Lambda_{2, \epsilon}  \pmb{\Theta}. 
#  \end{aligned}
# ```

# # OSRC Far Field
# In this example, we demonstrate the use of an OSRC operator as approximate mapping in the high-frequency regime, by validating the corresponding far field.
# We begin with loading needed packages
using BEAST, CompScienceMeshes
using SphericalScattering

using LinearAlgebra

using Plots 
# Start with the definiton of the needed physical properties:
c = 2.99792458e8      #speed of light
f = 8/2*c             #frequency
κ = 2π * f / c        #wavenumber
λ = c / f             #wavelength
η = 1.0;

# Next, we define the geometry and create the mesh.
h = λ/8
Γ = meshsphere(1.0, h);

# Define also basis functions on the mesh. 
# For the OSRC operators, we can use basis functions on the primal mesh, namely the Raviart-Thomas and Nédélec basis functions.

X = raviartthomas(Γ);
Nd = n×X;

# Define now the operators and the excitation

T = Maxwell3D.singlelayer(wavenumber=κ)
Einc = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
Hinc = -1/(im*κ*η) * curl(Einc);

# Next, assemble the OSRC EtM operator. The OSRC operator can be defined by its number of Padé terms `N_p`, 
# the curvature used in the complex damping parameter and the rotation angle of the poles θ_p (commonly chosen as π/2).
Np = 10                                         # number of Padé terms
θ_p = pi/2                                      # rotation angle of the poles
curvature = 1.0                                 # curvature parameter
etm = BEAST.EtmOsrc(κ, Np, θ_p, curvature)
etm_map = assemble(etm, X, X);

# We now approximate the scattered magnetic field by mapping the trace of the scattered electric field (in primal coefficients) onto the trace of the scattered magnetic field (in primal coefficients) using the OSRC EtM map.
e = (n × Einc)
exsc = -assemble(e, X)                          # trace of the scattered electric field on PEC surface
Gxx = assemble(BEAST.Identity(), X, X)
iG = BEAST.lu(Gxx)
exsc_primal = iG * exsc
mapped_h_Nd_primal = etm_map * exsc_primal;

# Also add the incident magnetic field (in primal coefficients), to obtain the total magnetic field
h = (n × Hinc)
h_Nd = assemble(h, Nd)
h_Nd_primal = iG * h_Nd
h_tot_OSRC = mapped_h_Nd_primal .+ h_Nd_primal;

# Compute now the far field using the OSRC solution.
Θ, Φ = range(0.0,stop=π,length=1002), 0.0
points_cartFF = [point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for θ in Θ for ϕ in Φ]
points_cartFF_OSRC = points_cartFF[1:10:end]
FF_OSRC = -im * f / (2 * c) * potential(MWFarField3D(T), points_cartFF_OSRC, h_tot_OSRC, Nd);

# For reference, compute the analytical far field using the Mie series (from SphericalScattering.jl)
sp = PECSphere(radius=1.0)
ex = planeWave(frequency=f)
FF = scatteredfield(sp, ex, FarField(points_cartFF));

# Finally, compare the OSRC far field with the Mie series.
# We observe that the OSRC solution is a very good approximation of the reference solution, considering its cheap computational cost.
# This corresponds to Fig. 4 (right) from [Approximate local magnetic-to-electric surface operators for time-harmonic Maxwell’s equations (C. Geuzaine et al. (2014))](https://www.researchgate.net/publication/261636307_Approximate_local_magnetic-to-electric_surface_operators_for_time-harmonic_Maxwell's_equations)
Θ_plot = collect(180/pi * Θ)
Plots.plot(Θ_plot,  norm.(FF), label="MIE series", lw=2, yscale=:log10)
Plots.plot!(Θ_plot[1:10:end], norm.(FF_OSRC), label="OSRC", marker=:utriangle)
Plots.xlabel!("θ (degrees)")
Plots.ylabel!("Far field")
