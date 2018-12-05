

using Compat

mutable struct MonomialBasis{T,Degree,NF} <: RefSpace{T,NF} end

valuetype(::MonomialBasis{T}) where {T} = T
degree(::MonomialBasis{T,D}) where {T,D} = D
numfunctions(x::MonomialBasis) = degree(x)+1
function (x::MonomialBasis)(p)
    t = cartesian(p)[1]
    valuetype(x)[t^d for d in 0 : degree(x)]
end

abstract type AbstractTimeBasisFunction end

"""
    TimeBasisFunction{N,D}

T: the value type of the time basis function
N: the number of intervals in the support (this included the semi infinite interval
    stretching to +∞)
D1: the degree of the TBF restricted to each of the intervals **plus one**
"""
mutable struct TimeBasisFunction{T,N,D1,D} <: AbstractTimeBasisFunction
    timestep::T
    numfunctions::Int
    polys::SVector{N,Polynomial{D1,T}}
end


function derive(tbf::TimeBasisFunction{T,N,D1,D}) where {T,N,D1,D}
    dpolys = Polynomial{D1-1,T}[ derive(p) for p in tbf.polys ]
    TimeBasisFunction{T,N,D1-1,D-1}(
        tbf.timestep,
        tbf.numfunctions,
        SVector{N,Polynomial{D1-1,T}}(dpolys...)
    )
end


function integrate(tbf::TimeBasisFunction{T,NI,D1,D}) where {T,NI,D1,D}

    timestep = tbf.timestep
    numfunctions = tbf.numfunctions

    polys = Vector{Polynomial{D1+1,T}}(undef, length(tbf.polys))
    c = zero(timestep)
    for i in 1:length(polys)
        t0 = (i-2) * timestep
        t1 = (i-1) * timestep
        polys[i] = integrate(tbf.polys[i],t0,c)
        c = polys[i](t1)
    end

    TimeBasisFunction{T,NI,D1+1,D+1}(
        tbf.timestep,
        tbf.numfunctions,
        SVector(polys...)
    )
end


mutable struct TimeBasisDelta{T} <: AbstractTimeBasisFunction
    timestep::T
    numfunctions::Int
end

mutable struct DiracBoundary{T} <: RefSpace{T,1} end
numfunctions(x::DiracBoundary) = 1

scalartype(x::TimeBasisDelta{T}) where {T} = T
numfunctions(x::TimeBasisDelta) = x.numfunctions
refspace(x::TimeBasisDelta{T}) where {T} = DiracBoundary{T}()
timestep(x::TimeBasisDelta) = x.timestep
numintervals(x::TimeBasisDelta) = 1


timebasisdelta(dt, Nt)= TimeBasisDelta(dt, Nt)

numintervals(tbf::TimeBasisFunction{T,N,D1,D}) where {T,N,D1,D} = N
degree(tbf::TimeBasisFunction{T,N,D1,D}) where {T,N,D1,D} = D
timestep(tbf::TimeBasisFunction) = tbf.timestep
scalartype(tbf::TimeBasisFunction{T,N,D1,D}) where {T,N,D1,D} = T

numfunctions(t::TimeBasisFunction) = t.numfunctions
refspace(t::TimeBasisFunction{T,N,D1,D}) where {T,N,D1,D} = MonomialBasis{T,D,D1}()




geometry(t::AbstractTimeBasisFunction) = [SegmentedAxis(timestep(t), numfunctions(t)+numintervals(t)-3)]


"""
    timebasisc0d1(type, timestep, numfunctions)

Build the space of continuous, piecewise linear time basis functions. The DoFs
are the time steps. `numfunctions` basis functions will be built in total.
"""
function timebasisc0d1(timestep, numfunctions, T::Type=Float64)
    i, z = one(T), zero(T)
    polys = SVector(
        Polynomial(SVector(i,  i/timestep)), # 1 + t
        Polynomial(SVector(i, -i/timestep)), # 1 - t
        Polynomial(SVector(z,  z))
    )
    TimeBasisFunction{T,3,2,1}(timestep, numfunctions, polys)
end

"""
    timebasiscxd0(timestep, numfunctions, T::Type=Float64)

Create a temporal basis based on shifted copies of the nodal continuous, piecewise
linear interpolant.
"""
function timebasiscxd0(timestep, numfunctions, T::Type=Float64)
    i, z = one(T), zero(T)
    polys = SVector(
        Polynomial(i), # 1
        Polynomial(z), # 0
    )
    TimeBasisFunction{T,2,1,0}(timestep, numfunctions, polys)
end


"""
    timebasisspline2(timestep, numfunctions, T::Type=Float64)

Create a temporal basis based on shifted copies of the quadratic spline. The
spline is the convolution of a cxd0 and a c0d1 basis function.
"""
function timebasisspline2(dt, numfunctions, T::Type=Float64)
    i, z = one(T), zero(T)
    polys = SVector{4,Polynomial{3,T}}(
        Polynomial(SVector(i/2, i/dt, i/2/dt/dt)),
        Polynomial(SVector(i/2, i/dt, -i/dt/dt)),
        Polynomial(SVector(2*i, -2*i/dt, i/2/dt/dt)),
        Polynomial(SVector(z, z, z))
    )
    TimeBasisFunction{T,4,3,2}(dt, numfunctions, polys)
end

function timebasisshiftedlagrange(dt, numfunctions, degree, T::Type=typeof(dt))
    z, i = zero(T), one(T)
    c = Polynomial(SVector(i))
    t = Polynomial(SVector(z,i/dt))
    polys = Polynomial{degree+1,T}[]
    for k = 0:degree
        f = c
        for i in 1:k
            f = (1/i) * f * (i - t)
        end
        g = c
        for i in 1:(degree-k)
            g = (1/i) * g * (i + t)
        end
        push!(polys, f*g)
    end
    push!(polys, 0*Base.power_by_squaring(t, degree))
    @assert length(polys) == degree+2
    polys = SVector{degree+2,Polynomial{degree+1,T}}(polys...)
    TimeBasisFunction{T,degree+2,degree+1,degree}(dt, numfunctions, polys)
end

function (f::TimeBasisFunction)(t::Real)
    dt = timestep(f)
    ni = numintervals(f)
    t < -dt && return zero(t)
    i = floor(Int, (t/dt)) + 2
    p = i <= ni ? f.polys[i] : f.polys[end]
    return p(t)
end





function assemblydata(tbf::TimeBasisFunction)

    T = scalartype(tbf)
    Δt = timestep(tbf)
    z = (0, zero(T))

    t = Polynomial(zero(T), one(T))

    num_cells = numfunctions(tbf)
    num_refs  = degree(tbf)+1

    max_num_funcs = numintervals(tbf)
    numfuncs = zeros(Int, num_cells, num_refs)
    data = fill(z, max_num_funcs, num_refs, num_cells)

    els = [ simplex(point((i-1)*Δt),point(i*Δt)) for i in 1:num_cells ]

    for k in 1 : numfunctions(tbf)
        tk = (k-1) * Δt
        for i in 1 : numintervals(tbf)
            # Focus on interval [(i-2)Δt,(i-1)Δt]
            p = tbf.polys[i]
            q = substitute(p,t-tk)

            c = k + i - 2
            1 <= c <= num_cells || continue
            for d = 0 : degree(q)
                r = d + 1
                w = q[d]

                j = (numfuncs[c,r] += 1)
                data[j,r,c] = (k,w)
            end
        end
    end

    return els, AssemblyData(data)
end



function temporalassemblydata(tbf; kmax=typemax(Int))

    T = scalartype(tbf)
    Δt = timestep(tbf)

    t = Polynomial(zero(T), one(T))

    num_cells = numfunctions(tbf)
    num_refs  = degree(tbf)+1

    has_zero_tail = all(tbf.polys[end].data .== 0)
    @show has_zero_tail

    if has_zero_tail
        max_num_funcs = numintervals(tbf) # 3
    else
        max_num_funcs = min(kmax, numfunctions(tbf))
    end

    numfuncs = zeros(Int, num_cells, num_refs)
    data = fill((0,zero(T)), max_num_funcs, num_refs, num_cells)
    for k in 1 : min(kmax, numfunctions(tbf))
        tk = (k-1) * Δt
        for i in 1 : numintervals(tbf)-1
            p = tbf.polys[i]
            q = substitute(p,t+tk)

            c = k - i + 1
            c < 1 && continue
            for d = 0 : degree(q)
                r = d + 1
                w = q[d]

                j = (numfuncs[c,r] += 1)
                data[j,r,c] = (k,w)
            end
        end

        # Treat the half-open interval separately
        if !has_zero_tail
            p = tbf.polys[numintervals(tbf)]
            q = substitute(p, t+tk)
            c = k - numintervals(tbf) + 1
            while c > 0
                for d = 0 : degree(q)
                    r = d+1
                    w = q[d]
                    j = (numfuncs[c,r] += 1)
                    data[j,r,c] = (k,w)
                end
                c = c-1
            end
        end
    end

    return AssemblyData(data)
end


function assemblydata(tbf::TimeBasisDelta)

    T = scalartype(tbf)
    Δt = timestep(tbf)

    z = zero(scalartype(tbf))
    w = one(scalartype(tbf))

    num_cells = numfunctions(tbf)
    num_refs  = 1

    max_num_funcs = 1
    num_funcs = zeros(Int, num_cells, num_refs)
    data = fill((0,z), max_num_funcs, num_refs, num_cells)

    els = [ simplex(point((i-0)*Δt),point((i+1)*Δt)) for i in 1:num_cells ]

    for k in 1 : numfunctions(tbf)-1
        data[1,1,k] = (k+1,w)
    end

    return els, AssemblyData(data)
end


function convolve(f::TimeBasisFunction, g::TimeBasisFunction)
    @info "BEAST.convolve: computing quadratic spline out of a pulse and a hat"
    dt = timestep(f)
    fg = timebasisspline2(dt, numfunctions(f), scalartype(f))
    fg.polys = dt * fg.polys
    return fg
end

# function convolve(f::TimeBasisFunction{T,2,1,0}, g::TimeBasisFunction{T,2,1,0}) where {T}
#     @info "Convolving two pulses into a hat"
#     dt = timestep(f)
#     fg = timebasisshiftedlagrange(dt, numfunctions(f), 1, scalartype(f))
#     fg.polys = dt * fg.polys
#     return fg
# end

function convolve(δ::TimeBasisDelta, g::TimeBasisFunction)
    @info "Computing the trivial convolution between a delta and something else"
    return g
end

function convolve(f::TimeBasisFunction{T,2,1,0}, g::TimeBasisFunction{T,2,1,0}) where {T}
    @info "Convolving two pulses into a hat"
    dt = BEAST.timestep(f)
    fg = timebasisshiftedlagrange(dt, numfunctions(f), 1, scalartype(f))
    fg.polys = dt * fg.polys
    return fg
end

function BEAST.convolve(
    f::TimeBasisFunction{T,2,1,0},
    g::TimeBasisFunction{T,N,D1,D}) where {T,N,D1,D}

    @info "Convolving a left-of-origin pulse with a piecewise poly"

    @assert degree(f.polys[1]) == 0
    @assert timestep(f) ≈ timestep(g)

    # The result of convolving with a LoO pulse extends the support
    # by Δt to the right.
    Δt = timestep(f)
    P = Polynomial{D1+1,T}
    polys = Vector{P}(undef, N+1)
    t = Polynomial{2,T}(SVector{2,T}(0, 1))

    # Treat the leftmost interval separately
    polys[1] = integrate(g.polys[1], -Δt, 0)

    for i in 2:N
        q1 = -1 * substitute(integrate(g.polys[i-1], (i-2)*Δt, 0), t-Δt)
        q2 = integrate(g.polys[i], (i-2)*Δt, 0)
        polys[i] = q1 + q2
    end

    q1 = -1 * substitute(integrate(g.polys[N], (N-1)*Δt, 0), t-Δt)
    q2 = integrate(g.polys[N], (N-1)*Δt, 0)
    polys[N+1] = q1 + q2

    TimeBasisFunction{T,N+1,D1+1,D+1}(
        g.timestep,
        g.numfunctions,
        polys)
end

"""
    TimeBasisDeltaShifted{T}

Represents a TimeBasisDelta{T} retarded by a fraction of the time step.
"""
struct TimeBasisDeltaShifted{T} <: AbstractTimeBasisFunction
	tbf   :: TimeBasisDelta{T}
	shift :: T
end
scalartype(x::TimeBasisDeltaShifted) = scalartype(x.tbf)
numfunctions(x::TimeBasisDeltaShifted) = numfunctions(x.tbf)
refspace(x::TimeBasisDeltaShifted) = refspace(x.tbf)
timestep(x::TimeBasisDeltaShifted) = timestep(x.tbf)
numintervals(x::TimeBasisDeltaShifted) = numintervals(x.tbf)

function assemblydata(tbds::TimeBasisDeltaShifted)
	tbf = tbds.tbf

    T = scalartype(tbf)
    Δt = timestep(tbf)

    z = zero(scalartype(tbf))
    w = one(scalartype(tbf))

    num_cells = numfunctions(tbf)
    num_refs  = 1

    max_num_funcs = 1
    num_funcs = zeros(Int, num_cells, num_refs)
    data = fill((0,z), max_num_funcs, num_refs, num_cells)

    els = [ simplex(point((i-0+tbds.shift)*Δt),point((i+1+tbds.shift)*Δt)) for i in 1:num_cells ]

    for k in 1 : numfunctions(tbf)
        data[1,1,k] = (k,w)
    end

    return els, AssemblyData(data)
end
