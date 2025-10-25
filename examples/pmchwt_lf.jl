# # Low frequency stabilized PMCHWT
#
# We begin with loading needed packages
using CompScienceMeshes, BEAST
using SphericalScattering

using SparseArrays
using LinearAlgebra
using Printf

using Plots
using PlotlyDocumenter #hide

# We start with the definiton of the needed physical and material properties
c = 2.99792458e8      #speed of light
f = 10.0^(-40)*c/2π   #frequency
μ = 4π * 1e-7         #permeability
ϵ = 1/(c^2*μ)         #permittivity
κ = 2π * f / c        #wavenumber
λ = c / f             #wavelength

η = 1.0
ϵr′ = 3.0
μr′ = 1.0
  
κ′, η′ =√(ϵr′*μr′)*κ, η*√(μr′/ϵr′)
α, α′ = 1/η, 1/η′;
 
# Next, we define the geometry and create the mesh.
h=0.25
M = meshsphere(1.0, h; generator=:gmsh);

# Once the geometry has been created, the basis function space ``X`` can be defined on the mesh. 
# In this example, we use the Raviart-Thomas basis ``X`` and we will also needed a dual basis ``Y`` in the form of Buffa-Christiansen functions to build a preconditioner.
X = raviartthomas(M)
Y = buffachristiansen(M);

# Then comes the definition of the integral operators. We define separate operators for the weakly- and the hypersingular part of the Maxwell singlelayer operator, since they behave differently when acting on loops and stars. 
# In partiuclar we have that ``P_\Lambda T_h = 0`` and ``T_h P_\Lambda = 0``. Simimlar we have that the first term in the expansion of the Maxwell doublelayer vanishes when acted upon with loops. ``K=K_0+K_L`` with ``P_\Lambda K_0 P_\Lambda=0``.
Ts  = Maxwell3D.weaklysingular(wavenumber=κ)
Ts′ = Maxwell3D.weaklysingular(wavenumber=κ′)
Th  = Maxwell3D.hypersingular(wavenumber=κ)
Th′ = Maxwell3D.hypersingular(wavenumber=κ′)
K  = Maxwell3D.doublelayer(wavenumber=κ)
K′ = Maxwell3D.doublelayer(wavenumber=κ′)
KL  = BEAST.MWDoubleLayer3DLoop(im*κ)
KL′ = BEAST.MWDoubleLayer3DLoop(im*κ′);

# Definition of the excitation.
E = Maxwell3D.planewave(direction=ŷ, polarization=x̂, wavenumber=κ)
H = -1/(im*κ*η)*curl(E)

e = (n × E) × n
h = (n × H) × n;

# Definition of low frequency excitation. Similar to the doublelayer the lowest term in the expansion of the planewave vanishes when interacting with a loop.
Eex = Maxwell3D.planewaveExtractedKernel(direction=ŷ, polarization=x̂, wavenumber=κ)
Hex = -1/(im*κ*η)*curl(Eex)
  
ex = (n × Eex) × n
hx = (n × Hex) × n;

# Assemble local matrices
G = assemble(Identity(),X,X)
invG = BEAST.cholesky(G)

Gmix = assemble(NCross(),X,Y)
invGmix = BEAST.lu(Gmix)

nX =  numfunctions(X)
Z = spzeros(nX,nX);

# Assembling of the quasi-Helmholtz projectors: Here, the quasi-Helmholtz projectors are directly computed from the star matrix ``\Sigma`` and loop matrix ``\Lambda``. 
# We have the following ``P\Sigma = \Sigma ( \Sigma^T \Sigma )^+ \Sigma^T`` and ``P\Lambda = I - P\Sigma`` for the primal projectors, with the superscript ``+`` indicating the pseudo-inverse.
# The dual projectors are defined like this: ``\mathbb{P}\Lambda =  \Lambda ( \Lambda^T \Lambda )^+ \Lambda^T`` and ``\mathbb{P}\Sigma = I - \mathbb{P}\Lambda``.
PΣ = assemble(BEAST.PΣ(;compStrat = BEAST.Direct),X)
PΛ = assemble(BEAST.PΛ(;compStrat = BEAST.Direct),X)
ℙΣ = assemble(BEAST.ℙΣ(;compStrat = BEAST.Direct),X)
ℙΛ = assemble(BEAST.ℙΛ(;compStrat = BEAST.Direct),X);  

# Assemble the integral operators for PMCHWT
Tsn = assemble(Ts,X,X)    
Tsn′ = assemble(Ts′,X,X)    
Thn = assemble(Th,X,X)    
Thn′ = assemble(Th′,X,X)  

Kn = assemble(K,X,X)    
Kn′ = assemble(K′,X,X)    
KLn = assemble(KL,X,X)    
KLn′ = assemble(KL′,X,X);   

# Assemble the integral operators for preconditioner based on the dual basis ``Y``.
𝕋sn = assemble(Ts,Y,Y)    
𝕋hn = assemble(Th,Y,Y)    
𝕋sn′ = assemble(Ts′,Y,Y)    
𝕋hn′ = assemble(Th′,Y,Y);  

# Assemble the excitation vectors
eh = assemble(e,X)
hh = assemble(h,X)

exh = assemble(ex,X)
hxh = assemble(hx,X);

# ## Regular preconditioned PMCHWT
# 
# Building the preconditioner
M = [ invGmix'  Z
      Z         invGmix' ] *
    [ 𝕋sn+𝕋hn    Z
      Z          𝕋sn+𝕋hn ] * 
    [ invGmix Z
      Z       invGmix ]; 
# PMCHWT 
A = [ (η*(Tsn+Thn)+η′*(Tsn′+Thn′))  (-Kn-Kn′)  
      (Kn+Kn′)                      (α*(Tsn+Thn)+α′*(Tsn′+Thn′))]; 
# Right hand side
b = [ eh
      hh ]; 
# Solving the linear system using GMRES
u, stats= BEAST.solve(BEAST.GMRES(A; M=M, rtol=1e-8, verbose=0), (b));

# ## Low frequency stabilized PMCHWT
#
# Low frequency scaling factors
k = sqrt(κ)
ik = 1/k;
# Low-frequency stabilized preconditioner. The additional factor of ``1/\sqrt(\kappa)`` is to rescale the residual close to 1.
M =   ik *
    [ k*PΛ*invGmix'  -im*ik*PΣ*invGmix'  Z              Z                  
      Z               Z                  k*PΛ*invGmix' -im*ik*PΣ*invGmix'  ] *
          
    [ 𝕋sn+𝕋hn   𝕋sn  Z         Z  
      𝕋sn       𝕋sn  Z         Z 
      Z         Z    𝕋sn+𝕋hn   𝕋sn
      Z         Z    𝕋sn       𝕋sn  ] * 

    [ im*k*ℙΛ  Z       
      ik*ℙΣ    Z      
      Z        im*k*ℙΛ  
      Z        ik*ℙΣ     ]; 
#  Low-frequency stabilized PMCHWT system
A = [-im*ik*ℙΛ*invGmix   k*ℙΣ*invGmix   Z                    Z             
      Z                  Z             -im*ik*ℙΛ*invGmix     k*ℙΣ*invGmix ] *
   
    [ (η*(Tsn)+η′*(Tsn′))  (η*(Tsn)+η′*(Tsn′))              (-KLn-KLn′)           (-Kn-Kn′)                         
      (η*(Tsn)+η′*(Tsn′))  (η*(Tsn+Thn)+η′*(Tsn′+Thn′))     (-Kn-Kn′)             (-Kn-Kn′)                         
      (KLn+KLn′)           (Kn+Kn′)                         (α*(Tsn)+α′*(Tsn′))   (α*(Tsn)+α′*(Tsn′))              
      (Kn+Kn′)             (Kn+Kn′)                         (α*(Tsn)+α′*(Tsn′))   (α*(Tsn+Thn)+α′*(Tsn′+Thn′)) ] *
          
    [ ik*PΛ    Z         
      im*k*PΣ  Z         
      Z        ik*PΛ    
      Z        im*k*PΣ ]; 

# And the right hand side
b = [-im*ik*ℙΛ*invGmix   k*ℙΣ*invGmix   Z                    Z             
      Z                  Z             -im*ik*ℙΛ*invGmix     k*ℙΣ*invGmix ] *
   
    [ exh
      eh
      hxh
      hh  ]; 

# Solving the linear system using GMRES
y, stats_lf = BEAST.solve(BEAST.GMRES(A;M=M,rtol=1e-8,verbose=0), (b));

# Recovering actual solution
u_lf =[ ik*PΛ+im*k*PΣ Z 
        Z             ik*PΛ+im*k*PΣ]  * y;

# ## Postprocessing
# Generate far field points
Φ, Θ = [0.0], range(0,stop=π,length=180)
ffpts = [point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for ϕ in Φ for θ in Θ];
    
# Compute far field for regular PMCHWT
j = u[1:nX]
m = u[nX+1:2nX]

ffm = potential(MWFarField3D(wavenumber=κ), ffpts, m, X)
ffj = potential(MWFarField3D(wavenumber=κ), ffpts, j, X)
ff =  -im * κ / 4π *( -η[1]*ffj + cross.(ffpts, ffm))

biRCS_mom = 4π*norm.(ff).^2/λ^2;

# Compute far field for low frequency PMCHWT. Keeping loop and star contributions seperate until the end.
jl = ik*PΛ* y[1:nX]
js = im*k*PΣ*  y[1:nX]
ms = im*k*PΣ* y[nX+1:2nX]
ml = ik*PΛ* y[nX+1:2nX]

ffms = potential(MWFarField3D(wavenumber=κ), ffpts, ms, X)
ffjs = potential(MWFarField3D(wavenumber=κ), ffpts, js, X)

ffml = potential(BEAST.MWFarField3DDropConstant(im*κ,1.0), ffpts, ml, X)
ffjl = potential(BEAST.MWFarField3DDropConstant(im*κ,1.0), ffpts, jl, X)
ffsl =  -im * κ / 4π *( -η[1]*(ffjs+ffjl) + cross.(ffpts, (ffms+ffml)))

biRCS_mom_sl = 4π*norm.(ffsl).^2/λ^2;

# Computing reference solution (Mie series)
sp = DielectricSphere(
        radius      = 1.0,
        filling     = Medium(ϵ*ϵr′,μ*μr′)
    )

ex = planeWave(frequency=f)

biRCS = rcs(sp,ex,ffpts) / λ^2;

# Plot radar cross section
Plots.plot(Θ*180/π,10*log10.(biRCS),label="Mie Series")
Plots.scatter!(Θ*180/π,10*log10.(biRCS_mom), label="PMCHWT",markershape=:star)
Plots.scatter!(Θ*180/π,10*log10.(biRCS_mom_sl), label="Low frequency PMCHWT",markershape=:star)
Plots.plot!(xlabel="θ in degree",ylabel="RCS/λ² in dB",xlim=(0,180),legend=:outerbottom, title="Bistatic Radar Cross-Section at $(@sprintf("%.2e",f)) Hz")
