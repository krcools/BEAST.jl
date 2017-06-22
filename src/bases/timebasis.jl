export timebasisc0d1
export timebasiscxd0
export timebasisdelta
export timebasisshiftedlagrange

using Compat

type MonomialBasis{T,Degree,NF} <: RefSpace{T,NF} end

valuetype{T}(::MonomialBasis{T}) = T
degree{T,D}(::MonomialBasis{T,D}) = D
numfunctions(x::MonomialBasis) = degree(x)+1
@compat function (x::MonomialBasis)(p)
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
type TimeBasisFunction{T,N,D1,D} <: AbstractTimeBasisFunction
    timestep::T
    numfunctions::Int
    polys::SVector{N,Polynomial{D1,T}}
end


function derive{T,N,D1,D}(tbf::TimeBasisFunction{T,N,D1,D})
    dpolys = Polynomial{D1-1,T}[ derive(p) for p in tbf.polys ]
    TimeBasisFunction{T,N,D1-1,D-1}(
        tbf.timestep,
        tbf.numfunctions,
        SVector{N,Polynomial{D1-1,T}}(dpolys...)
    )
end


type TimeBasisDelta{T} <: AbstractTimeBasisFunction
    timestep::T
    numfunctions::Int
end

type DiracBoundary{T} <: RefSpace{T,1} end
numfunctions(x::DiracBoundary) = 1

scalartype{T}(x::TimeBasisDelta{T}) = T
numfunctions(x::TimeBasisDelta) = x.numfunctions
refspace{T}(x::TimeBasisDelta{T}) = DiracBoundary{T}()
timestep(x::TimeBasisDelta) = x.timestep
numintervals(x::TimeBasisDelta) = 1


timebasisdelta(dt, Nt)= TimeBasisDelta(dt, Nt)

numintervals{T,N,D1,D}(tbf::TimeBasisFunction{T,N,D1,D}) = N
degree{T,N,D1,D}(tbf::TimeBasisFunction{T,N,D1,D}) = D
timestep(tbf::TimeBasisFunction) = tbf.timestep
scalartype{T,N,D1,D}(tbf::TimeBasisFunction{T,N,D1,D}) = T

numfunctions(t::TimeBasisFunction) = t.numfunctions
refspace{T,N,D1,D}(t::TimeBasisFunction{T,N,D1,D}) = MonomialBasis{T,D,D1}()




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

function timebasisshiftedlagrange(dt, numfunctions, degree, T::Type=Float64)
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
    push!(polys, 0*(t^degree))
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

    for k in 1 : numfunctions(tbf)
        data[1,1,k] = (k,w)
    end

    return els, AssemblyData(data)
end


function convolve(f::TimeBasisFunction, g::TimeBasisFunction)
    warn("BEAST.convolve only returns correct result for constant * continuous,linear")
    dt = timestep(f)
    fg = timebasisspline2(dt, numfunctions(f), scalartype(f))
    fg.polys = dt * fg.polys
    return fg
end

function convolve(δ::TimeBasisDelta, g::TimeBasisFunction)
    return g
end
