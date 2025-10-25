# # Higher-order Multiplicative Calderon Preconditioning
# 
# We begin with loading the needed packages  
using CompScienceMeshes, BEAST
using LinearAlgebra
using PlotlyJS
using PlotlyDocumenter #hide
# Then we start with defining the geometry and creating the mesh.
h=0.5
M = meshsphere(1.0, h; generator=:gmsh);
# Once the geometry has been created, the basis function space ``X`` can be defined on the mesh. 
# In this example, we use basis functions of order ``p=1``. Further, higher-order basis space ``Y`` will be needed to build a preconditioner.
p=1
X = BEAST.gwpdiv(M;order=p)
Y = BEAST.gwpdiv(M;order=p+2);
# Defintion of material properties
κ,  η  = 1.0, 1.0
κ′, η′ = √2.0κ, η/√2.0
α, α′ = 1/η, 1/η′;
# Defintion of integral operators and excitation.
T = Maxwell3D.singlelayer(wavenumber=κ)
T′ = Maxwell3D.singlelayer(wavenumber=κ′)
K  = Maxwell3D.doublelayer(wavenumber=κ)
K′ = Maxwell3D.doublelayer(wavenumber=κ′)
E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ);
H = -1/(im*κ*η)*curl(E)

e = (n × E) × n
h = (n × H) × n;
# Next, we increase the default number of quadrature points, since higher order basis are used.  
BEAST.@defaultquadstrat (T,X,X) BEAST.DoubleNumSauterQstrat(7,7,6,6,6,6)
BEAST.@defaultquadstrat (K,X,X) BEAST.DoubleNumSauterQstrat(7,7,6,6,6,6)
BEAST.@defaultquadstrat (e,X) BEAST.SingleNumQStrat(10)
BEAST.@defaultquadstrat (h,X) BEAST.SingleNumQStrat(10);

# Here is the core that allows Calderon preconditioning for a higher order basis, without using dual functions.
# The ``\Theta``-operators enable mapping a lower-order solenoidal subspace into a higher-order nonsolenoidal subspace and vice-versa.
# They are defined like this:
#
# ```math
# \begin{aligned}
#   \Theta_\Sigma  &= P_\Sigma^r N_{rp} P_\Lambda^p \\
#   \Theta_\Lambda &= P_\Lambda^r N_{rp} P_\Sigma^p
# \end{aligned}
# ```
#
# where PΣ and PΛ are the quasi-Helmholtz projectors, N the mixed Gram matrix of the lower and higher order basis, and the super-/subscripts indicate the order of the basis.

ΘΣ = BEAST.ThetaStars()
ΘΛ = BEAST.ThetaLoops();
# With everthing defined and setup, the linear system can be built.
@hilbertspace j m
@hilbertspace r s;
# In a first step, the standard PMCHWT system is assembled based on the basis ``X``. 
Ah = assemble( (η*T+η′*T′)[r,j] -      (K+K′)[r,m] +
                    (K+K′)[s,j] + (α*T+α′*T′)[s,m], s∈X, r∈X, j∈X, m∈X);
# Next, we assemble the block diagonal preconditioner using the higher order basis ``Y``.
Ch = assemble(T[r,j]+T[s,m],s∈Y, r∈Y, j∈Y, m∈Y);
# The ``\Theta``-Operator can be assembled like other operators. Mapping from the lowere order basis ``X`` to the higher order basis ``Y``. 
Θh = assemble( ΘΣ[r,j] + ΘΛ[r,j] + ΘΣ[s,m] + ΘΛ[s,m] ,s∈Y, r∈Y, j∈X, m∈X);

# And last but not least, the right hand side is assembled too.
bh = assemble(-e[j]-h[m],j ∈ X, m∈X);

# After everything has been assembled the linear system can be solved using GMRES.
# Once using the higher order preconditioner 
u, stats = solve(BEAST.GMRES(Ah; M=Θh'*Ch*Θh, restart=true, atol=1e-8, rtol=1e-8, verbose=0, memory=50),bh);
# And without using a preconditioner
uref, statsref = solve(BEAST.GMRES(Ah; restart=true, atol=1e-8, rtol=1e-8, verbose=0, memory=500),bh);
# The number of iterations for the preconditioned system: 
@show stats.niter #hide
# And the regular system:
@show statsref.niter #hide
# Visualizationf of the solution
fcrj, _ = facecurrents(u[j],X)
fcrm, _ = facecurrents(u[m],X);

# Plot magnetic surface current
plt = plot(patch(M, norm.(fcrj)),Layout(title="PMCHWT m"))
PlotlyDocumenter.to_documenter(plt) #hide
# Plot electric surface current
plt = plot(patch(M, norm.(fcrm)),Layout(title="PMCHWT j"))
PlotlyDocumenter.to_documenter(plt) #hide
