@testitem "testing TDFunctional" begin

    using CompScienceMeshes
    
    Γ = meshcuboid(1.0, 1.0, 1.0, 2.0; generator=:gmsh)
    X = raviartthomas(Γ)

    struct ConstFunctional{T} <: Functional{T}
        constant::T
    end

    function(f::ConstFunctional)(r)
        f.constant * [1.0, 0.0, 0.0]
    end

    struct FuncXGaussian{T} <: TDFunctional{T}
        functional::ConstFunctional{T}
        gaussian::BEAST.Gaussian{T}
    end

    function(f::FuncXGaussian)(r,t)
        r = cartesian(r)
        t = cartesian(t)[1]

        f.functional(r) * f.gaussian(t)     
    end

    Δt, Nt = 0.5, 20
    duration = 4 * Δt                                
    delay = 10 * Δt                                        
    
    func = ConstFunctional(1.0)
    gaussian = creategaussian(duration, delay)
    constxgaussian = FuncXGaussian(func, gaussian)


    ### Testing with Dirac delta 
    δ = timebasisdelta(Δt, Nt)
    A = assemble(constxgaussian, X⊗δ)
    # Analytically temporal testing
    Ats = assemble(func, X) * gaussian.(Δt*[1:1:Nt;])'

    @test maximum(abs.(A - Ats)) < 1e-12


    ### Testing with pulse functions
    p = timebasiscxd0(Δt, Nt) 	                			                      
    B = assemble(constxgaussian, X⊗p)
    igaussian = integrate(gaussian)
    Bts = assemble(func, X) * (igaussian.(Δt*[1:1:Nt;]) - igaussian.(Δt*[0:1:Nt-1;]))'
    
    @test maximum(abs.(B - Bts)) < 1e-12

    
    ### Testing with hat functions
    import SpecialFunctions: erf 
    function integratewithh(g::BEAST.Gaussian, ts, Δt)
        A = g.scaling
        t0 = g.delay
        w = g.width
        
        r = zeros(length(ts))

        for i in ts
            a = (i-1)*Δt
            x = i*Δt
            b = (i+1)*Δt

            ua = 4*(a-t0)/w
            ux = 4*(x-t0)/w
            ub = 4*(b-t0)/w

            r[i] = A/(2*Δt) * ((t0 - a)*(erf(ux) - erf(ua)) + (b - t0)*(erf(ub) - erf(ux))) + A*w/(8*√(π)*Δt)*(exp(-ua^2) + exp(-ub^2) - 2*exp(-ux^2))
        end
        return r
    end

    h = timebasisc0d1(Δt, Nt) 	 
    C = assemble(constxgaussian, X⊗h)
    Cts = assemble(func, X) * integratewithh(gaussian, [1:1:Nt;], Δt)'

    @test maximum(abs.(C - Cts)) < 1e-12
end
