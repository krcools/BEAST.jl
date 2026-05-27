using LoopVectorization

abstract type Operations end
struct Cross <: Operations end
struct Dot <: Operations end
struct Times <: Operations end
struct STimesS <: Operations end
struct VTimesS <: Operations end
struct STimesV <: Operations end
function (::Cross)(a,b)
    return cross(a,b)
end
function (::Dot)(a,b)
    return transpose(a)*b
end
function (::Union{Times,STimesS,VTimesS,STimesV})(a,b) 
    return a*b
end

function generate_cacheSIMDDoubleNum(operator::CompDoubleKern{T,Op1,Op2,Kern}, test_quad_data, trial_quad_data) where {T,Op1,Kern,Op2}
    a = length(test_quad_data[1,1].weights)
    b = length(trial_quad_data[1,1].weights)
    return [zeros(a,b) for _ in 1:cachesizeSIMDDoubleNum(Kern)]
end
function generate_cacheSIMDDoubleNum(operator::MWSingleLayer3D, test_quad_data, trial_quad_data) 
    a = length(test_quad_data[1,1].weights)
    b = length(trial_quad_data[1,1].weights)
    return [zeros(a,b) for _ in 1:2]
end
function generate_cacheSIMDDoubleNum(operator::MWDoubleLayer3D, test_quad_data, trial_quad_data) 
    a = length(test_quad_data[1,1].weights)
    b = length(trial_quad_data[1,1].weights)
    return [zeros(a,b) for _ in 1:6]
end
cachesizeSIMDDoubleNum(::Type{<:HH3DGreen}) = 2
cachesizeSIMDDoubleNum(::Type{<:HH3DGradGreen}) = 6
function quadpoints_simd(f, els, rules, derivative_needed)
    if derivative_needed 
        return map(Iterators.product(rules, els)) do (rule, el)
        pws = quadpoints(el, rule)
        (weights = [w for (p,w) in pws], 
        points = [cartesian(p)[i] for i in 1:length(cartesian(pws[1][1])), (p,w) in pws], 
        values = [f(p)[j][1][i] for  i in 1:length(f(pws[1][1])[1][1]), (p,w) in pws,j in 1:length(f(pws[1][1]))], 
        derivatives =[f(p)[j][2][i] for  i in 1:length(f(pws[1][1])[1][2]), (p,w) in pws,j in 1:length(f(pws[1][1]))])
        end
    else
        return map(Iterators.product(rules, els)) do (rule, el)
            pws = quadpoints(el, rule)
            (weights = [w for (p,w) in pws], 
            points = [cartesian(p)[i] for i in 1:length(cartesian(pws[1][1])), (p,w) in pws], 
            values = [f(p)[j][1][i] for  i in 1:length(f(pws[1][1])[1][1]), (p,w) in pws,j in 1:length(f(pws[1][1]))])
        end
    end
    
end

requires_derivative(::MWSingleLayer3D) = true
requires_derivative(a) = false


struct SIMDDoubleQuadRule{P,Q}
  outer_quad_points::P
  inner_quad_points::Q
  g_cache::Vector{Matrix{Float64}}
end


struct SIMDDoubleNumQStrat{R} <: AbstractQuadStrat
    outer_rule::R 
    inner_rule::R
end

"""
momintegrals!(biop, tshs, bshs, tcell, bcell, interactions, strat)

Function for the computation of moment integrals using simple double quadrature.
"""
function momintegrals!(biop::MWDoubleLayer3D,
    tshs, bshs, tcell, bcell, z, strat::SIMDDoubleQuadRule)
    momintegrals!(CompDoubleKern(Dot(),BEAST.HH3DGradGreen(biop.gamma),Cross()), tshs, bshs, tcell, bcell, z, strat::SIMDDoubleQuadRule)
end

function momintegrals!(biop::MWSingleLayer3D,
    tshs, bshs, tcell, bcell, z, strat::SIMDDoubleQuadRule)

    igd = Integrand(biop, tshs, bshs, tcell, bcell)

    womps = strat.outer_quad_points
    wimps = strat.inner_quad_points

    αim = imag(igd.operator.α)
    αre = real(igd.operator.α)
    βim = imag(igd.operator.β)
    βre = real(igd.operator.β)
    γim = imag(igd.operator.gamma)
    γre = real(igd.operator.gamma)
    @turbo for wimpid in eachindex(wimps.weights)
            bgeoxi = wimps.points[1,wimpid]
            bgeoyi = wimps.points[2,wimpid]
            bgeozi = wimps.points[3,wimpid]
           for wompid in eachindex(womps.weights)
                tgeoxi = womps.points[1,wompid]
                tgeoyi = womps.points[2,wompid]
                tgeozi = womps.points[3,wompid]
                rx = tgeoxi-bgeoxi; ry = tgeoyi-bgeoyi; rz = tgeozi-bgeozi
                R = sqrt(rx*rx + ry*ry + rz*rz)
                iR = 1 / R
                e = exp(-γre*R)
                s, c = sincos(γim*R)
                scale = e * (i4pi * iR)
                greenre = c * scale
                greenim = -s * scale

                strat.g_cache[1][wompid, wimpid] = greenre 
                strat.g_cache[2][wompid, wimpid] = greenim
            end
        end
    for l in 1:3 
        for k in 1:3
            zre = 0.0
            zim = 0.0
            @turbo for wimpid in eachindex(wimps.weights)
                    bvalsxi = wimps.values[1,wimpid,l]
                    bvalsyi = wimps.values[2,wimpid,l]
                    bvalszi = wimps.values[3,wimpid,l]
                    bdiv = wimps.derivatives[wimpid,l]
                    jy = wimps.weights[wimpid]
            
            for wompid in eachindex(womps.weights)
                tvalsxi = womps.values[1,wompid,k]
                tvalsyi = womps.values[2,wompid,k]
                tvalszi = womps.values[3,wompid,k]
                tdiv = womps.derivatives[wompid,k]
                jx = womps.weights[wompid]
                    j = jx * jy
                    greenre = strat.g_cache[1][wompid, wimpid]
                    greenim = strat.g_cache[2][wompid, wimpid]
                    dotprod = tvalsxi*bvalsxi + tvalsyi*bvalsyi + tvalszi*bvalszi
                    tdivbdiv = tdiv * bdiv
                    A_re = αre*dotprod + βre*tdivbdiv
                    A_im = αim*dotprod + βim*tdivbdiv

                    zre += j * (greenre*A_re - greenim*A_im)
                    zim += j * (greenre*A_im + greenim*A_re)
                end
            end
 
            z[k,l] += zre + im*zim


        end
    end
   
    return z
end

TransposedStrat(a::SIMDDoubleQuadRule) = a


@generated function momintegrals!(biop::CompDoubleKern{T,Op1,Op2,Kern},
    tshs::TSpace, bshs::BSpace, tcell, bcell, z, strat::SIMDDoubleQuadRule) where {T,Op1,Kern,Op2,TSpace,BSpace}
    M = numfunctions(tshs, domain(tcell))
    N = numfunctions(bshs, domain(bcell))
    # Core.println("M = $M, N = $N")
    # load variables
    vars = _load_variables(Kern)

    kernexp = _kernexp(Kern)

    assign_green_realexp = assign_green_real(Kern)
    assign_green_imagexp = assign_green_imag(Kern)



    get_greenrealexp = get_greenreal(Kern)
    get_greenimagexp = get_greenimag(Kern)


    op1 = operationexp(Op1)
    op2 = operationexp(Op2)
    trialspacetovarmapexp = trialspacetovarmap(BSpace)
    testspacetovarmapexp = testspacetovarmap(TSpace)
    loadvaluestestexp = loadvalues_test(TSpace)
    loadvaluestrialexp = loadvalues_trial(BSpace)
    ex = quote 
        womps = strat.outer_quad_points
        wimps = strat.inner_quad_points
        $vars
       @turbo for wimpid in eachindex(wimps.weights)
            bgeoxi = wimps.points[1,wimpid]
            bgeoyi = wimps.points[2,wimpid]
            bgeozi = wimps.points[3,wimpid]
           for wompid in eachindex(womps.weights)
                tgeoxi = womps.points[1,wompid]
                tgeoyi = womps.points[2,wompid]
                tgeozi = womps.points[3,wompid]

                $kernexp # evaluates the kernel between the given points
                $assign_green_realexp # assigns the real part of the kernel to the cache
                $assign_green_imagexp # assigns the imaginary part of the kernel to the cache

           end
        end
        for li in 1:$N 
            for ki in 1:$M
                zre = 0.0
                zim = 0.0
            for wimpid in eachindex(wimps.weights)
                            $loadvaluestrialexp # loads the trial space values at the given quadrature point into variables
                            jy = wimps.weights[wimpid]
                    
                            for wompid in eachindex(womps.weights)
                                $loadvaluestestexp # loads the test space values at the given quadrature point into variables
                                jx = womps.weights[wompid]
                                j = jx * jy

                                $get_greenrealexp # retrieves the real part of the kernel from the cache
                                $trialspacetovarmapexp # maps the trial space values to the working variables used in the next operation
                                $op2 # evaluates the second operation (e.g. multiplying by the trial space values)
                                $testspacetovarmapexp # maps the test space values to the working variables used in the next operation
                                $op1 # evaluates the first operation (e.g. multiplying by the test space values)
                                zre += j * r

                                $get_greenimagexp # retrieves the imaginary part of the kernel from the cache
                                $trialspacetovarmapexp # maps the trial space values to the working variables used in the next operation
                                $op2 # evaluates the second operation (e.g. multiplying by the trial space values)
                                $testspacetovarmapexp # maps the test space values to the working variables used in the next operation
                                $op1 # evaluates the first operation (e.g. multiplying by the test space values)
                                zim += j * r

                            end
                        end
                z[ki,li] += zre + im*zim
            end

        end

            return z
    end
    # display(ex)
    return ex
end
valuetype(::Type{RTRefSpace{T}}) where {T} = SVector{3,T}
CompScienceMeshes.domain(a::Type{CompScienceMeshes.Simplex{A,B,C,D,T}}) where {A,B,C,D,T} = CompScienceMeshes.ReferenceSimplex{B} 

function _load_variables(::Type{U}) where {U <: Union{HH3DGreen, HH3DGradGreen}}
    quote
        γim = imag(biop.kernel.gamma)
        γre = real(biop.kernel.gamma)
    end
end
function loadvalues_test(::Type{U}) where {U}
    if valuetype(U) <: Number
        return quote
            tvals = womps.values[1, wompid,ki]
        end
    end
    if valuetype(U) <: SVector{3}
        return quote
            tvalsxi = womps.values[1,wompid,ki]
            tvalsyi = womps.values[2,wompid,ki]
            tvalszi = womps.values[3,wompid,ki]
        end
    end
end
function loadvalues_trial(::Type{U}) where {U}
    if valuetype(U) <: Number
        return quote
            bvals = wimps.values[1,wimpid,li]
        end
    end
    if valuetype(U) <: SVector{3}
        return quote
            bvalsxi = wimps.values[1,wimpid,li]
            bvalsyi = wimps.values[2,wimpid,li]
            bvalszi = wimps.values[3,wimpid,li]
        end
    end
end
function _kernexp(::Type{<:HH3DGreen})
    quote
        rx = tgeoxi-bgeoxi; ry = tgeoyi-bgeoyi; rz = tgeozi-bgeozi
        R = sqrt(rx*rx + ry*ry + rz*rz)
        iR = 1 / R
        e = exp(-γre*R)
        s, c = sincos(γim*R)
        scale = e * (i4pi * iR)
        greenre = c * scale
        greenim = -s * scale
    end
end
function _kernexp(::Type{<:HH3DGradGreen})
    quote
        rx = tgeoxi-bgeoxi; ry = tgeoyi-bgeoyi; rz = tgeozi-bgeozi
        R = sqrt(rx*rx + ry*ry + rz*rz)
        iR = 1 / R
        e = exp(-γre*R)
        s, c = sincos(γim*R)
        scale = e * (i4pi * iR)
        greenre = c*scale
        greenim = -s*scale

        gradgreenre = iR*((-γre-iR)*greenre + γim*greenim)
        gradgreenim = iR*((-γre-iR)*greenim - γim*greenre)

       
        greenrex = gradgreenre * rx
        greenrey = gradgreenre * ry
        greenrez = gradgreenre * rz
        greenimx = gradgreenim * rx
        greenimy = gradgreenim * ry
        greenimz = gradgreenim * rz
    end
end
operationexp(::Type{Cross}) = _cross_exp()
operationexp(::Type{Dot}) = _dot_exp()
operationexp(::Type{STimesS}) = _scalar_times_scalar_exp()
operationexp(::Type{VTimesS}) = _vector_times_scalar_exp()
operationexp(::Type{STimesV}) = _scalar_times_vector_exp()

function trialspacetovarmap(a::Type{T}) where {T <: RefSpace}

    if valuetype(T) <: SVector{3}
        ex = quote
            rx = bvalsxi; ry = bvalsyi; rz = bvalszi
        end
    elseif valuetype(T) <: Number
        ex = quote
            r = bvals
        end
    else
        # k = _size(valuetype(T))
        # Core.println(k)
        error("unsupported trial space value type")
    end
    return ex
end
function testspacetovarmap(a::Type{T}) where {T <: RefSpace}

    if valuetype(T) <: SVector{3}
        ex = quote
            lx = tvalsxi; ly = tvalsyi; lz = tvalszi
        end
    elseif valuetype(T) <: Number
        ex = quote
            l = tvals
        end
    else
        error("unsupported trial space value type")
    end
    return ex
end
""" 
    valuetype(::Type{<:Space})

return the type of the values of the given space evaluated in a meshpoint
"""

valuetype(::Type{_LocalBasisDot{T}}) where {T} = T
function valuetype(::Type{_LocalBasisTimes{T,U,V}}) where {T,U,V}

    valuetype(U) <: Number && valuetype(V) <: Number && return promote_type(U,V)
    return SVector{3,T}
end
valuetype(::Type{_LocalBasisCross{T,U,V}}) where {T,U,V} = SVector{3,T}
valuetype(::Type{NormalVector}) = SVector{3,T}

# _size(a::Type{<:SVector}) = size(a)
# _size(a::Type{<:Number}) = 1
function _scalar_times_scalar_exp()
    quote
        sol = l*r
        r = sol
    end
end
function _scalar_times_vector_exp()
    quote
        solx = l*rx
        soly = l*ry
        solz = l*rz
        rx = solx
        ry = soly
        rz = solz
    end
end
function _vector_times_scalar_exp()
    quote
        solx = lx*r
        soly = ly*r
        solz = lz*r
        rx = solx
        ry = soly
        rz = solz
    end
end
function _dot_exp()
    quote
        r = rx*lx + ry*ly + rz*lz
    end
end
function _cross_exp()
    quote
        solx = -ry*lz + rz*ly
        soly = -rz*lx + rx*lz
        solz = -rx*ly + ry*lx
        rx = solx
        ry = soly
        rz = solz
    end
end

function get_greenreal(kern::Type{<:HH3DGreen})
    quote
        l = strat.g_cache[1][wompid, wimpid]
    end
end
function get_greenimag(kern::Type{<:HH3DGreen})
    quote
        l = strat.g_cache[2][wompid, wimpid]
    end
end
function get_greenreal(kern::Type{<:HH3DGradGreen})
    quote
        lx = strat.g_cache[1][wompid, wimpid]
        ly = strat.g_cache[2][wompid, wimpid]
        lz = strat.g_cache[3][wompid, wimpid]
    end
end
function get_greenimag(kern::Type{<:HH3DGradGreen})
    quote
        lx = strat.g_cache[4][wompid, wimpid]
        ly = strat.g_cache[5][wompid, wimpid]
        lz = strat.g_cache[6][wompid, wimpid]
    end 
end
function assign_green_real(kern::Type{<:HH3DGreen})
    quote
        strat.g_cache[1][wompid, wimpid] = greenre
    end
end
function assign_green_imag(kern::Type{<:HH3DGreen})
    quote
        strat.g_cache[2][wompid, wimpid] = greenim
    end
end
function assign_green_real(kern::Type{<:HH3DGradGreen})
    quote
        strat.g_cache[1][wompid, wimpid] = greenrex
        strat.g_cache[2][wompid, wimpid] = greenrey
        strat.g_cache[3][wompid, wimpid] = greenrez
    end
end
function assign_green_imag(kern::Type{<:HH3DGradGreen})
    quote
        strat.g_cache[4][wompid, wimpid] = greenimx
        strat.g_cache[5][wompid, wimpid] = greenimy
        strat.g_cache[6][wompid, wimpid] = greenimz
    end
end


function quadrule(operator::IntegralOperator,
    local_test_basis, local_trial_basis,
    test_id, test_element, trial_id, trial_element, qd, qs::SIMDDoubleNumQStrat)
    qr = SIMDDoubleQuadRule(qd.tpoints[1, test_id], qd.bpoints[1, trial_id], qd.cache)
end
function quaddata(operator::IntegralOperator,
    local_test_basis, local_trial_basis,
    test_elements, trial_elements, qs::SIMDDoubleNumQStrat) 
    derivative_needed = requires_derivative(operator)
    (tpoints = quadpoints_simd(local_test_basis, test_elements, (qs.outer_rule,),derivative_needed), 
     bpoints = quadpoints_simd(local_trial_basis, trial_elements, (qs.inner_rule,),derivative_needed))
end

function quadcache(operator::IntegralOperator,
    local_test_basis, local_trial_basis,
    test_elements, trial_elements, qd, qs::SIMDDoubleNumQStrat) 
    g_cache = generate_cacheSIMDDoubleNum(operator, qd.tpoints, qd.bpoints)
    return (tpoints = qd.tpoints, bpoints = qd.bpoints, cache = g_cache,)
end

# frequently used quadstrats
SIMDDoubleNumSauterQStrat(a,b,c,d,e,f) = SauterQStrat(SIMDDoubleNumQStrat(a,b),c,d,e,f)

@testitem "SIMD DoubleNumSauter Quadrature" begin
    using CompScienceMeshes

    Γ = meshcuboid(1.0,1.0,1.0,0.3)
    X = raviartthomas(Γ);

    ops = [BEAST.CompDoubleInt(x->x, BEAST.Dot(),BEAST.HH3DGreen(1.0*im),BEAST.STimesV(), x->x),
        BEAST.CompDoubleInt(x->x, BEAST.Dot(),BEAST.HH3DGradGreen(1.0*im),BEAST.Cross(), x->x),
        Maxwell3D.singlelayer(wavenumber = 1.0),
        Maxwell3D.doublelayer(wavenumber = 1.0)]
    for op2 in ops
        M3 = assemble(op2,X,X,quadstrat = BEAST.DoubleNumSauterQstrat(6,6,6,6,6,6));
        M4 = assemble(op2,X,X,quadstrat = BEAST.SauterQStrat(BEAST.SIMDDoubleNumQStrat(6,6),6,6,6,6));
        @test isapprox(M3,M4, atol = 1e-5)
    end
end