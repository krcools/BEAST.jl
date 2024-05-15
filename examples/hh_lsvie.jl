using BEAST
using CompScienceMeshes
using StaticArrays
using LinearAlgebra
using Test
using Plots
using SphericalScattering

# Layered Dielectic Sphere
# MoM (Lippmann Schwinger VIE) vs. Analytical Solution (SphericalScattering)



# Environment
ε1 = 1.0*ε0
μ1 = μ0

# Layered Dielectic Sphere
ε2 = 5.0*ε0 # outer shell
μ2 = μ0

ε3 = 20.0*ε0
μ3 = μ0

ε4 = 7.0*ε0
μ4 = μ0

ε5 = 30.0*ε0 # inner core
μ5 = μ0

r = 1.0
radii = SVector(0.2*r, 0.6*r, 0.8*r, r) # inner core stops at first radius
filling = SVector(Medium(ε5, μ5), Medium(ε4, μ4), Medium(ε3, μ3), Medium(ε2, μ2))




# Mesh, Basis
h = 0.1
mesh = CompScienceMeshes.tetmeshsphere(r,h)
X = lagrangec0d1(mesh; dirichlet = false)
@show numfunctions(X)

bnd = boundary(mesh)
strc = X -> strace(X, bnd)


# VIE Operators
function generate_tau(ε_env, radii, filling)

    function tau(x::SVector{U,T}) where {U,T}
        for (i, radius) in enumerate(radii)
            norm(x) <= radius && (return T(1.0 - filling[i].ε/ε_env))
        end
        #return 0.0
        error("Evaluated contrast outside the sphere!")
    end

    return tau
end
τ = generate_tau(ε1, radii, filling) #contrast function

I = Identity()
V = VIE.hhvolume(tau = τ, wavenumber = 0.0)
B = VIE.hhboundary(tau = τ, wavenumber = 0.0)
Y = VIE.hhvolumegradG(tau = τ, wavenumber = 0.0)


# Exitation
dirE = SVector(1.0, 0.0, 0.0)
dirgradΦ = - dirE
amp = 1.0
Φ_inc = VIE.linearpotential(direction = dirgradΦ, amplitude = amp)


# Assembly
b = assemble(Φ_inc, X)

Z_I = assemble(I, X, X)

Z_V = assemble(V, X, X)
Z_B = assemble(B, strc(X), X)

Z_Y = assemble(Y, X, X)


# System matrix
Z_version1 = Z_I + Z_V + Z_B
Z_version2 = Z_I + Z_Y


# Solution of the linear system of equations
u_version1 = Z_version1 \ Vector(b)
u_version2 = Z_version2 \ Vector(b)



## Observe the potential on a x-line in dirgradΦ direction

# x-line at y0, z0 - Potential only inside the sphere mesh valid!
y0 = 0.0
z0 = 0.0
x_range = range(0.0, stop = 1.0*r, length = 200)
points_x = [point(x, y0, z0) for x in x_range]
x = collect(x_range)

# SphericalScattering solution on x-line
sp = LayeredSphere(radii = radii, filling = filling)
ex = UniformField(direction = dirE, amplitude = amp)
Φ_x = real.(field(sp, ex, ScalarPotential(points_x)))

# MoM solution on x-line
Φ_MoM_version1_x = BEAST.grideval(points_x, u_version1, X, type = Float64)
Φ_MoM_version2_x = BEAST.grideval(points_x, u_version2, X, type = Float64)

# Plot
plot(x, Φ_x, label = "SphericalScattering")
plot!(x, Φ_MoM_version1_x, label = "MoM: S=I+B+V, boundary+volume")
plot!(x, Φ_MoM_version2_x, label = "MoM: S=I+Y, gradgreen")
xlims!(0.0, 1.0) # sphere center to radius r
ylims!(-0.3, 0.0)
title!("Potential Φ(x, y0, z0)")
xlabel!("x")