
# test resolutuion of #66

using BEAST
using Test
using CompScienceMeshes
using LinearAlgebra

function hassemble2(operator::BEAST.AbstractOperator,
    test_functions,
    trial_functions)

    blkasm = BEAST.blockassembler(operator, test_functions, trial_functions)

    function assembler2(Z, tdata, sdata)
        store(v,m,n) = (Z[m,n] += v)
        blkasm(tdata,sdata,store)
    end

    mat = zeros(scalartype(operator), 
                numfunctions(test_functions), 
                numfunctions(trial_functions))

    assembler2(mat, 1:numfunctions(test_functions), 1:numfunctions(trial_functions))
    return mat
end

for T in [Float32, Float64]
    c = T(3e8)
    μ = T(4π * 1e-7)
    ε = T(1/(μ*c^2))
    f = T(1e8)
    λ = T(c/f)
    k = T(2π/λ)
    ω = T(k*c)
    η = T(sqrt(μ/ε))

    a = T(1)
    Γ = CompScienceMeshes.meshcuboid(a,a,a,T(0.2); generator=:compsciencemeshes)

    𝓣 = Maxwell3D.singlelayer(wavenumber=k)
    𝓚 = Maxwell3D.doublelayer(wavenumber=k)

    X = raviartthomas(Γ)
    Y = buffachristiansen(Γ)

    # println("Number of RWG functions: ", numfunctions(X))

    T_blockassembler = hassemble2(𝓣, X, X)
    T_standardassembler = assemble(𝓣, X, X)

    @test norm(T_blockassembler - T_standardassembler)/norm(T_standardassembler) ≈ 0.0 atol=100*eps(T)

    T_bc_blockassembler = hassemble2(𝓣, Y, Y)
    T_bc_standardassembler = assemble(𝓣, Y, Y)

    @test norm(T_bc_blockassembler - T_bc_standardassembler)/norm(T_bc_standardassembler) ≈ 0.0 atol=100*eps(T)

    K_mix_blockassembler = hassemble2(𝓚,Y,X)
    K_mix_standardassembler = assemble(𝓚,Y,X)

    T_mix_blockassembler = hassemble2(𝓣, Y, X)
    T_mix_standardassembler = assemble(𝓣, Y, X)

    if T==Float64
        @test norm(K_mix_blockassembler - K_mix_standardassembler)/norm(K_mix_standardassembler) ≈ 0.0 atol=100*eps(T)
        @test norm(T_mix_blockassembler - T_mix_standardassembler)/norm(T_mix_standardassembler) ≈ 0.0 atol=100*eps(T)
        @info "Tests executed!"
    end
end
