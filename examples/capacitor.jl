using BEAST
using CompScienceMeshes
using LinearAlgebra

# define geometry & material parameters
f₀= 3.0e6 # 3 GHz
c = 3.0e8
ϵ₀ = 8.85418782e-12
μ₀ = 1.25663706e-6
ϵᵣ=1; μᵣ=1
l = w = 1.0 #Length and width of capacitor plates
d = 0.1     #seperation of plates
h = 1/12    #size of meshes
ω = 2π*f₀
vₑ = c/sqrt(ϵᵣ)
λₑ = vₑ/f₀
η = sqrt(μ₀*μᵣ/ϵ₀*ϵᵣ)
κ = ω/vₑ

println("wavenumber = ", κ)
println("freq = ", f₀, " Hz")
println("λ = ", λₑ, "m; l = ", l ,"m")

# Build and mesh the structure
Γ₁ = meshrectangle(l,w,h);
translate!(Γ₁, point(0.0,0.0,d))
Γ₀ = meshrectangle(l,w,h)

γ₁ = meshsegment(l, h, 3)
translate!(γ₁, point(0.0,0.0,d))
γ₀ = meshsegment(l, h, 3)

Γ = weld(Γ₁, Γ₀)

# define the excitation
V₀ = 1.0
f = ScalarTrace(p -> V₀)

#define basis function
RT = rt_ports(Γ,[γ₁,γ₀])

# define the equations
@hilbertspace j
@hilbertspace k
T = MWSingleLayer3D(im*κ)
EFIE = @varform η*T[k,j] == f[ntrace(k,γ₁)]

# discretise & solve the equation
efie = @discretise EFIE j∈RT k∈RT
u = solve(efie)

#Calculate & plot face currents
fcr,geo = facecurrents(u, RT); println("Face currents calculated")
# patch_mat(geo, fcr); println("Facecurrents plotted")

# Compute the Scalar Potential across the ports
z = range(-0.5, stop=0.5, length=100)
pts1 = point.(0.5,0.5,z)
# pts1 = [point(0.5,0.5,a) for a in range(-0.5,stop=0.5,length=100)]
volt = η*potential(MWSingleLayerPotential3D(κ), pts1, u, RT)

# Plot the Scalar potential
using Plots
Φ = real.(getindex.(volt,1))
display(plot(z,Φ,xlabel="height",ylabel="scalar potential",label=""))

# Compare numerical capacitance with parallel plate approximation
idx = getindex_rtg(RT) #get index of global RWG
Q = norm(-u[idx]/(im*ω))
C = (ϵ₀*ϵᵣ*(l*w))/d
println("Relative different w.r.t. large plate capacitance: ", 2*abs(Q-C)/abs(Q+C))
