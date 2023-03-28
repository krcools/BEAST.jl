using CompScienceMeshes
using BEAST
using LinearAlgebra
import Plotly

# Exterior wavenumber
κ₀ = 3.0

# This is where the number of domains enters the problem description
κ = [1.5κ₀, 2.5κ₀]

# Description of the domain boundaries
h = 0.15
Γ1 = meshcuboid(0.5, 1.0, 1.0, h)
Γ2 =  -Mesh([point(-x,y,z) for (x,y,z) in vertices(Γ1)], deepcopy(cells(Γ1)))
Γ = [Γ1, Γ2]

# Beyond this point the problem description is completely independent
# of the number of domains and their relative positioning

# Incident field
E = Maxwell3D.planewave(direction=(x̂+ẑ)/√2, polarization=ŷ, wavenumber=κ₀)
H = -1/(im*κ₀)*curl(E)

# Definition of the boundary integral operators
T0 = Maxwell3D.singlelayer(wavenumber=κ₀)
K0 = Maxwell3D.doublelayer(wavenumber=κ₀)

T = [Maxwell3D.singlelayer(wavenumber=κᵢ) for κᵢ ∈ κ]
K = [Maxwell3D.doublelayer(wavenumber=κᵢ) for κᵢ ∈ κ]

@hilbertspace m j
@hilbertspace k l
A0 = K0[k,m] - T0[k,j] + T0[l,m] + K0[l,j]
A = [Kᵢ[k,m] - Tᵢ[k,j] + Tᵢ[l,m] + Kᵢ[l,j] for (Tᵢ,Kᵢ) ∈ zip(T,K)]
Adiag = BEAST.Variational.DirectProductKernel(A)

N = BEAST.NCross()
Nᵢ = -0.5*N[k,m] - 0.5*N[l,j]
Ndiag = BEAST.Variational.BlockDiagKernel(Nᵢ)

# This needs a clean interface but for now it is the
# only way to keep the number of domains a runtime variable
psyms = [Symbol(:p,i) for i in eachindex(κ)]
qsyms = [Symbol(:q,i) for i in eachindex(κ)]
p = [BEAST.Variational.HilbertVector(i,psyms,[]) for i ∈ eachindex(psyms)]
q = [BEAST.Variational.HilbertVector(i,qsyms,[]) for i ∈ eachindex(qsyms)]

B = A0[p,q] + Adiag[p,q] + Ndiag[p,q] - Nᵢ[p,q]

e = (n × E) × n
h = (n × H) × n
bᵢ = e[k] - h[l]
b = bᵢ[p]

Xₕ = raviartthomas.(Γ)
Yₕ = raviartthomas.(Γ)

Pₕ = [Xᵢ×Yᵢ for (Xᵢ,Yᵢ) ∈ zip(Xₕ,Yₕ)]
Qₕ = [Xᵢ×Yᵢ for (Xᵢ,Yᵢ) ∈ zip(Xₕ,Yₕ)]
deq = BEAST.discretise(B==b, (p .∈ Pₕ)..., (q .∈ Qₕ)...)
# u = solve(deq)
u = gmres(deq, tol=1e-4, maxiter=2500)

fcrm = [facecurrents(u[qᵢ][m], Xᵢ)[1] for (qᵢ,Xᵢ) ∈ zip(q,Xₕ)]
fcrj = [facecurrents(u[qᵢ][j], Xᵢ)[1] for (qᵢ,Xᵢ) ∈ zip(q,Xₕ)]

Plotly.plot([patch(Γᵢ, norm.(fcrᵢ), caxis=(0,2)) for (fcrᵢ,Γᵢ) in zip(fcrm,Γ)])
Plotly.plot([patch(Γᵢ, norm.(fcrᵢ), caxis=(0,2)) for (fcrᵢ,Γᵢ) in zip(fcrj,Γ)])

function nearfield(um,Xm,uj,Xj,κ,η,points)

    K = BEAST.MWDoubleLayerField3D(wavenumber=κ)
    T = BEAST.MWSingleLayerField3D(wavenumber=κ)

    Em = potential(K, points, um, Xm)
    Ej = potential(T, points, uj, Xj)
    E = -Em + η * Ej

    Hm = potential(T, points, um, Xm)
    Hj = potential(K, points, uj, Xj)
    H = 1/η*Hm + Hj

    return E, H
end

Xs = range(-1.5,1.5,length=150)
Zs = range(-0.5,1.5,length=100)
pts = [point(x,0.5,z) for z in Zs, x in Xs]

EHi = [nearfield(-u[qᵢ][m], Xᵢ, -u[qᵢ][j], Yᵢ, κᵢ, 1.0, pts) for (qᵢ,Xᵢ,Yᵢ,κᵢ) ∈ zip(q,Xₕ,Yₕ,κ)]
EH0 = [nearfield(u[qᵢ][m], Xᵢ, u[qᵢ][j], Yᵢ, κ₀, 1.0, pts) for (qᵢ,Xᵢ,Yᵢ) ∈ zip(q,Xₕ,Yₕ)]

Eo = sum(getindex.(EH0,1)) + E.(pts)
Ho = sum(getindex.(EH0,2)) + H.(pts)

Ei = getindex.(EHi,1)
Hi = getindex.(EHi,2)

Etot = Eo + sum(Ei)
Htot = Ho + sum(Hi)

import Plots
Plots.heatmap(Xs, Zs, clamp.(abs.(getindex.(Etot,2)),-2.0,2.0); colormap=:viridis)