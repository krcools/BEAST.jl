using CompScienceMeshes
using BEAST
using LinearAlgebra
using Test

# Γ = meshcuboid(1.0, 1.0, 1.0, 0.25)
# Γ = meshsphere(1.0, 0.3)
Γ = meshrectangle(1.0, 1.0, 1/3, 3)
Γ1 = CompScienceMeshes.rotate(Γ, pi/2*point(1,0,0))
Γ2 = CompScienceMeshes.weld(Γ,Γ1)

sol = 1.0
SL0 = TDMaxwell3D.singlelayer(speedoflight=sol, numdiffs=0)
SL1 = TDMaxwell3D.singlelayer(speedoflight=sol, numdiffs=1)
DL0 = TDMaxwell3D.doublelayer(speedoflight=sol, numdiffs=0)

Δt, Nt = 0.6, 80
duration = 20 * Δt
delay = 1.5 * duration
amplitude = 1.0
gaussian = creategaussian(duration, delay, amplitude)
direction, polarisation = ẑ, x̂
E = BEAST.planewave(polarisation, direction, gaussian, 1.0)
E1 = BEAST.planewave(polarisation, direction, derive(gaussian), 1.0)


δ = timebasisdelta(Δt, Nt)
T0 = timebasisshiftedlagrange(Δt, Nt, 0)
T1 = timebasisshiftedlagrange(Δt, Nt, 1)
T2 = timebasisshiftedlagrange(Δt, Nt, 2)
T3 = timebasisshiftedlagrange(Δt, Nt, 3)

iT1 = integrate(T1)
iT2 = integrate(T2)

X = raviartthomas(Γ)
Y = buffachristiansen(Γ, ibscaled=false)

Z1, store1 = BEAST.allocatestorage(SL0, X⊗δ, X⊗T1,
    BEAST.Val{:bandedstorage},
    BEAST.LongDelays{:compress})
Z2, store2 = BEAST.allocatestorage(SL0, X⊗δ, X⊗T1,
    BEAST.Val{:bandedstorage},
    BEAST.LongDelays{:ignore})

BEAST.assemble!(SL0, X⊗δ, X⊗T1, store1)
BEAST.assemble!(SL0, X⊗δ, X⊗T1, store2)

@test norm(Z1-Z2,Inf) < 1e-12

Z3, store7 = BEAST.allocatestorage(SL1, X⊗δ, X⊗T2,
    BEAST.Val{:bandedstorage},
    BEAST.LongDelays{:compress})
Z4, store8 = BEAST.allocatestorage(SL1, X⊗δ, X⊗T2,
    BEAST.Val{:bandedstorage},
    BEAST.LongDelays{:ignore})

BEAST.assemble!(SL1, X⊗δ, X⊗T2, store7)
BEAST.assemble!(SL1, X⊗δ, X⊗T2, store8)



# x = rand(size(Z1)[1:2]...)
# csx = cumsum(x, dims = 2)
#
# j, k_start = 10, 2
# y1 = zeros(eltype(Z1), size(Z1,1))
# y2 = zeros(eltype(Z2), size(Z2,1))
# BEAST.convolve!(y1, Z1, x, csx, j, k_start)
# y2 = convolve(Z2, x, j, k_start)
# @show norm(y2 - y1)


b = assemble(E, X⊗δ)
b1 = assemble(E1, X⊗δ)

W1 = inv(Z1[:,:,1])
W2 = inv(Z2[:,:,1])
@test norm(W1-W2,Inf) < 1e-12

W3 = inv(Z3[:,:,1])
W4 = inv(Z4[:,:,1])
@test norm(W3-W4) < 1e-12

x1 = marchonintime(W1, Z1, b, Nt)
x2 = marchonintime(W1, Z2, b, Nt)
@test norm(x1-x2,Inf) < 1e-12

x3 = marchonintime(W3, Z3, b1, Nt)
x4 = marchonintime(W4, Z4, b1, Nt)
@test norm(x3-x4,Inf) < 1e-12

@test norm(vec(x1)-vec(x3)) / norm(vec(x3)) < 0.1
@test norm(vec(x1)-vec(x4)) / norm(vec(x4)) < 0.1

X2 = raviartthomas(Γ2)
DLh = (DL0, X2⊗δ, X2⊗iT2)

Z9, store9 = BEAST.allocatestorage(DLh...,
    BEAST.Val{:bandedstorage},
    BEAST.LongDelays{:compress})
Z10, store10 = BEAST.allocatestorage(DLh...,
    BEAST.Val{:bandedstorage},
    BEAST.LongDelays{:ignore})

@time BEAST.assemble!(DLh..., store9)
@time BEAST.assemble!(DLh..., store10)

@test norm(Z9-Z10,Inf) < 1e-12
