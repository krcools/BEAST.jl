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

fr1, store1 = BEAST.allocatestorage(SL0, X⊗δ, X⊗T1,
    BEAST.Val{:bandedstorage},
    BEAST.LongDelays{:compress})
fr2, store2 = BEAST.allocatestorage(SL0, X⊗δ, X⊗T1,
    BEAST.Val{:densestorage},
    BEAST.LongDelays{:ignore})

BEAST.assemble!(SL0, X⊗δ, X⊗T1, store1); Z1 = fr1()
BEAST.assemble!(SL0, X⊗δ, X⊗T1, store2); Z2 = fr2()

for k in axes(Z1,3)
    local T = scalartype(SL0, X⊗δ, X⊗T1)
    Z1k = zeros(T, size(Z1)[1:2])
    Z2k = zeros(T, size(Z2)[1:2])
    BEAST.ConvolutionOperators.timeslice!(Z1k, Z1, k)
    BEAST.ConvolutionOperators.timeslice!(Z2k, Z2, k)
    @test norm(Z1k .- Z2k, Inf) < 1e-12
end
# for m in axes(Z1,1)
#     for n in axes(Z1,2)
#         for k in axes(Z1,3)
#             @test norm(Z1[m,n,k] - Z2[m,n,k]) < 1e-12
# end end end
# @test norm(Z1-Z2,Inf) < 1e-12

fr3, store7 = BEAST.allocatestorage(SL1, X⊗δ, X⊗T2,
    BEAST.Val{:bandedstorage},
    BEAST.LongDelays{:compress})
fr4, store8 = BEAST.allocatestorage(SL1, X⊗δ, X⊗T2,
    BEAST.Val{:densestorage},
    BEAST.LongDelays{:ignore})

BEAST.assemble!(SL1, X⊗δ, X⊗T2, store7); Z3 = fr3()
BEAST.assemble!(SL1, X⊗δ, X⊗T2, store8); Z4 = fr4()



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

function timeslice(T,Z,k)
    Zk = zeros(T, size(Z)[1:2])
    BEAST.ConvolutionOperators.timeslice!(Zk, Z, k)
end

T = scalartype(SL0, X⊗δ, X⊗T1)
W1 = inv(timeslice(T,Z1,1))
W2 = inv(timeslice(T,Z2,1))

# W1 = inv(AbstractArray(Z1)[:,:,1])
# W2 = inv(AbstractArray(Z2)[:,:,1])
# W1 = inv(BEAST.timeslice(Z1,1))
# W2 = inv(BEAST.timeslice(Z2,1))
@test norm(W1-W2,Inf) < 1e-12

T = scalartype(SL1, X⊗δ, X⊗T2)
W3 = inv(timeslice(T,Z3,1))
W4 = inv(timeslice(T,Z4,1))

# W3 = inv(AbstractArray(Z3)[:,:,1])
# W4 = inv(AbstractArray(Z4)[:,:,1])
# W3 = inv(BEAST.timeslice(Z3,1))
# W4 = inv(BEAST.timeslice(Z4,1))
@test norm(W3-W4) < 1e-12

x1 = marchonintime(W1, Z1, b, Nt)
x2 = marchonintime(W2, Z2, b, Nt)
@test norm(x1-x2,Inf) < 1e-12

x3 = marchonintime(W3, Z3, b1, Nt)
x4 = marchonintime(W4, Z4, b1, Nt)
@test norm(x3-x4,Inf) < 1e-12

@test norm(vec(x1)-vec(x3)) / norm(vec(x3)) < 0.1
@test norm(vec(x1)-vec(x4)) / norm(vec(x4)) < 0.1

X2 = raviartthomas(Γ2)
Y2 = raviartthomas(Γ2)
DLh = (DL0, X2⊗δ, X2⊗iT2)

fr9, store9 = BEAST.allocatestorage(DLh...,
    BEAST.Val{:bandedstorage},
    BEAST.LongDelays{:compress})
fr10, store10 = BEAST.allocatestorage(DLh...,
    BEAST.Val{:densestorage},
    BEAST.LongDelays{:ignore})

@time BEAST.assemble!(DLh..., store9)
@time BEAST.assemble!(DLh..., store10)

Z9 = fr9()
Z10 = fr10()

for k in axes(Z9,3)
    local T = scalartype(DLh...)
    Z9k = zeros(T, size(Z9)[1:2])
    Z10k = zeros(T, size(Z10)[1:2])
    BEAST.ConvolutionOperators.timeslice!(Z9k, Z9, k)
    BEAST.ConvolutionOperators.timeslice!(Z10k, Z10, k)
    @test norm(Z9k .- Z10k, Inf) < 1e-12
end

# # @test norm(Z9-Z10,Inf) < 1e-12
# for m in axes(Z9,1)
#     for n in axes(Z9,2)
#         for k in axes(Z9,3)
#             @test norm(Z9[m,n,k] - Z10[m,n,k]) < 1e-12
# end end end

iDLh = (integrate(DL0), Y2⊗δ, X2⊗T1)
fr11, store11 = BEAST.allocatestorage(iDLh..., BEAST.Val{:densestorage}, BEAST.LongDelays{:ignore})
fr12, store12 = BEAST.allocatestorage(iDLh..., BEAST.Val{:bandedstorage}, BEAST.LongDelays{:compress})
# fr13, store13 = BEAST.allocatestorage(iDLh..., BEAST.Val{:bandedstorage}, BEAST.LongDelays{:ignore})


BEAST.assemble!(iDLh..., store11); Z11 = fr11()
BEAST.assemble!(iDLh..., store12); Z12 = fr12()
# BEAST.assemble!(iDLh..., store13); Z13 = fr13()

# @test norm(Z11-Z12,Inf) < 1e-8
# for m in axes(Z11,1)
#     for n in axes(Z11,2)
#         for k in axes(Z11,3)
#             @test norm(Z11[m,n,k] - Z12[m,n,k]) < 1e-12
# end end end
# @test norm(Z11-Z13,Inf) < 1e-8

for k in axes(Z11,3)
    local T = scalartype(iDLh...)
    Z11k = zeros(T, size(Z11)[1:2])
    Z12k = zeros(T, size(Z12)[1:2])
    BEAST.ConvolutionOperators.timeslice!(Z11k, Z11, k)
    BEAST.ConvolutionOperators.timeslice!(Z12k, Z12, k)
    @test norm(Z11k .- Z12k, Inf) < 1e-12
end