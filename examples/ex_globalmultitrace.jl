using CompScienceMeshes
using BEAST
using LinearAlgebra
import Plotly

# Exterior wavenumber
κ₀ = 3.0

# This is where the number of domains enters the problem description
κ = [1.5κ₀, 2.5κ₀]

# Description of the domain boundaries
h = 0.125
Γ1 = meshcuboid(0.5, 1.0, 1.0, h)
Γ2 =  -Mesh([point(-x,y,z) for (x,y,z) in vertices(Γ1)], deepcopy(cells(Γ1)))
Γ = [Γ1, Γ2]

# Beyond this point the problem description is completely independent
# of the number of domains and their relative positioning

# Incident field
Einc = Maxwell3D.planewave(direction=(x̂+ẑ)/√2, polarization=ŷ, wavenumber=κ₀)
Hinc = -1/(im*κ₀)*curl(Einc)

# Definition of the boundary integral operators
T0 = Maxwell3D.singlelayer(wavenumber=κ₀)
K0 = Maxwell3D.doublelayer(wavenumber=κ₀)

T = [Maxwell3D.singlelayer(wavenumber=κᵢ) for κᵢ ∈ κ]
K = [Maxwell3D.doublelayer(wavenumber=κᵢ) for κᵢ ∈ κ]

# Definition of per-domain bilinear forms
@hilbertspace m j
@hilbertspace k l
A0 = K0[k,m] - T0[k,j] + T0[l,m] + K0[l,j]
A = [Kᵢ[k,m] - Tᵢ[k,j] + Tᵢ[l,m] + Kᵢ[l,j] for (Tᵢ,Kᵢ) ∈ zip(T,K)]
N = BEAST.NCross()
Nᵢ = -0.5*N[k,m] - 0.5*N[l,j]

# Building the global system
p = BEAST.hilbertspace(:p, length(κ))
q = BEAST.hilbertspace(:q, length(κ))
Adiag = BEAST.Variational.DirectProductKernel(A)
Ndiag = BEAST.Variational.BlockDiagKernel(Nᵢ)
B = A0[p,q] + Adiag[p,q] + Ndiag[p,q] - Nᵢ[p,q]

# Also for the right hand side, first per domain contributions
# are constructed, followed by the global synthesis
e = (n × Einc) × n
h = (n × Hinc) × n
bᵢ = e[k] - h[l]
b = bᵢ[p]

Xₕ = raviartthomas.(Γ)
Yₕ = raviartthomas.(Γ)

Pₕ = [Xᵢ×Yᵢ for (Xᵢ,Yᵢ) ∈ zip(Xₕ,Yₕ)]
Qₕ = [Xᵢ×Yᵢ for (Xᵢ,Yᵢ) ∈ zip(Xₕ,Yₕ)]
deq = BEAST.discretise(B==b, (p .∈ Pₕ)..., (q .∈ Qₕ)...)
u = solve(deq)
# u = gmres(deq, tol=1e-4, maxiter=2500)

fcrm = [facecurrents(u[qᵢ][m], Xᵢ)[1] for (qᵢ,Xᵢ) ∈ zip(q,Xₕ)]
fcrj = [facecurrents(u[qᵢ][j], Xᵢ)[1] for (qᵢ,Xᵢ) ∈ zip(q,Xₕ)]

Plotly.plot([patch(Γᵢ, norm.(fcrᵢ), caxis=(0,2)) for (fcrᵢ,Γᵢ) ∈ zip(fcrm,Γ)])
Plotly.plot([patch(Γᵢ, norm.(fcrᵢ), caxis=(0,2)) for (fcrᵢ,Γᵢ) ∈ zip(fcrj,Γ)])

function nearfield(um,Xm,uj,Xj,κ,η,points)

    K = BEAST.MWDoubleLayerField3D(wavenumber=κ)
    T = BEAST.MWSingleLayerField3D(wavenumber=κ)

    Em = potential(K, points, um, Xm)
    Ej = potential(T, points, uj, Xj)

    Hm = potential(T, points, um, Xm)
    Hj = potential(K, points, uj, Xj)
    
    return -Em + η * Ej, 1/η*Hm + Hj
end

Xs = range(-2.0,2.0,length=150)
Zs = range(-1.5,2.5,length=100)
pts = [point(x,0.5,z) for z in Zs, x in Xs]

EH0 = [nearfield(u[qᵢ][m], Xᵢ, u[qᵢ][j], Yᵢ, κ₀, 1.0, pts) for (qᵢ,Xᵢ,Yᵢ) ∈ zip(q,Xₕ,Yₕ)]
EHi = [nearfield(-u[qᵢ][m], Xᵢ, -u[qᵢ][j], Yᵢ, κᵢ, 1.0, pts) for (qᵢ,Xᵢ,Yᵢ,κᵢ) ∈ zip(q,Xₕ,Yₕ,κ)]

Eo = sum(getindex.(EH0,1)) + Einc.(pts)
Ho = sum(getindex.(EH0,2)) + Hinc.(pts)

Ei = getindex.(EHi,1)
Hi = getindex.(EHi,2)

Etot = Eo + sum(Ei)
Htot = Ho + sum(Hi)

import Plots
Plots.heatmap(Xs, Zs, clamp.(real.(getindex.(Etot,2)),-2.0,2.0); colormap=:viridis)