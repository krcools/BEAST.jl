# test resolutuion of #66

using BEAST
using Test
using CompScienceMeshes
using LinearAlgebra

function hassemble(operator::BEAST.AbstractOperator,
    test_functions,
    trial_functions)

    blkasm = BEAST.blockassembler(operator, test_functions, trial_functions)

    function assembler(Z, tdata, sdata)
        store(v,m,n) = (Z[m,n] += v)
        blkasm(tdata,sdata,store)
    end

    mat = zeros(scalartype(operator), 
                numfunctions(test_functions), 
                numfunctions(trial_functions))

    assembler(mat, 1:numfunctions(test_functions), 1:numfunctions(trial_functions))
    return mat
end

c = 3e8
μ = 4π * 1e-7
ε = 1/(μ*c^2)
f = 1e8
λ = c/f
k = 2π/λ
ω = k*c
η = sqrt(μ/ε)

a = 1.0
Γ = CompScienceMeshes.meshcuboid(a,a,a,0.2)

𝓣 = Maxwell3D.singlelayer(wavenumber=k)
𝓚 = Maxwell3D.doublelayer(wavenumber=k)

X = raviartthomas(Γ)
Y = buffachristiansen(Γ)

println("Number of RWG functions: ", numfunctions(X))

T_blockassembler = hassemble(𝓣, X, X)
T_standardassembler = assemble(𝓣, X, X)

@test norm(T_blockassembler - T_standardassembler)/norm(T_standardassembler) ≈ 0.0 atol=1e-14

T_bc_blockassembler = hassemble(𝓣, Y, Y)
T_bc_standardassembler = assemble(𝓣, Y, Y)

@test norm(T_bc_blockassembler - T_bc_standardassembler)/norm(T_bc_standardassembler) ≈ 0.0 atol=1e-14

K_mix_blockassembler = hassemble(𝓚,Y,X)
K_mix_standardassembler = assemble(𝓚,Y,X)

T_mix_blockassembler = hassemble(𝓣, Y, X)
T_mix_standardassembler = assemble(𝓣, Y, X)

@test norm(K_mix_blockassembler - K_mix_standardassembler)/norm(K_mix_standardassembler) ≈ 0.0 atol=1e-14
@test norm(T_mix_blockassembler - T_mix_standardassembler)/norm(T_mix_standardassembler) ≈ 0.0 atol=1e-14