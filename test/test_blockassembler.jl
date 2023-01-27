using BEAST
using CompScienceMeshes
using LinearAlgebra
using Test

r = 10.0
λ = 20 * r
k = 2 * π / λ

sphere = readmesh(joinpath(dirname(@__FILE__),"assets","sphere5.in"), T=Float64)

D = Maxwell3D.doublelayer(wavenumber=k)
X = raviartthomas(sphere)
Y = buffachristiansen(sphere)

A = assemble(D, X, X)

@views blkasm = BEAST.blockassembler(D, X, X)

@views function assembler(Z, tdata, sdata)
    @views store(v,m,n) = (Z[m,n] += v)
    blkasm(tdata,sdata,store)
end

A_blk = zeros(ComplexF64, length(X.fns), length(Y.fns))
assembler(A_blk, [1:length(X.fns);], [1:length(Y.fns);])

@test norm(A - A_blk) ≈ 0 atol=eps(Float64) 
