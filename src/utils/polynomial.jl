using StaticArrays

struct Polynomial{N,T}
    data::SVector{N,T}
end
Polynomial(x::T...) where {T} = Polynomial{length(x),T}(SVector{length(x),T}(x))
Polynomial(a::T) where {T<:Number} = Polynomial{1,T}(SVector{1,T}(a))

Base.copy(p::Polynomial) = p
Base.one(p::Polynomial{N,T}) where{N,T} = Polynomial(@SVector[one(T)])

Base.length(p::Polynomial) = length(p.data)
Base.eltype(p::Polynomial) = eltype(p.data)
Base.getindex(p::Polynomial,i::Int) = p.data[i+1]
Base.setindex!(p::Polynomial,v,i::Int) = (p.data[i+1] = v)

degree(p::Polynomial) = length(p.data)-1

function derive(p)
    x = eltype(p)[d*p[d] for d in 1:degree(p)]
    Polynomial(x...)
end

function integrate(p::Polynomial{N,T},x0,y0) where {N,T}

    coeffs = similar(p.data, length(p.data)+1)
    fill!(coeffs, 0)
    for i in 2:length(coeffs)
        coeffs[i] = p.data[i-1]/(i-1)
    end
    temp = Polynomial(SVector(coeffs...))
    coeffs[1] = y0 - temp(x0)
    # TODO: Remove high order terms with zero coefficient
    return Polynomial(SVector(coeffs...))
end

function evaluate(p::Polynomial, t)
    d = p.data
    T = promote_type(typeof(t), eltype(d))
    r = zero(T)
    for i in 1:length(d)
        r += d[i] * t^(i-1)
    end
    return r
end

(p::Polynomial)(t) = evaluate(p,t)


function integrate(p)
    x = zeros(eltype(p), degree(p)+2)
    for i in 1:degree(p)+1
        x[i+1] = p[i-1]/i
    end
    Polynomial(x...)
end

import Base: *,+,-,promote_rule

#The sizes are static (ntuple over Val) so that the arithmetic below is type stable and does
#not allocate: temporalassemblydata does a lot of it for every block of a time-domain matrix.
function *(p::Polynomial{N},q::Polynomial{M}) where {N,M}
    T = promote_type(eltype(p), eltype(q))
    D = N+M-1
    x = ntuple(Val(D)) do i
        d = i-1
        s = zero(T)
        for k in max(0,d-M+1) : min(d,N-1)
            s += p[k] * q[d-k]
        end
        s
    end
    Polynomial{D,T}(SVector{D,T}(x))
end

# function *(p::Polynomial{1}, q::Polynomial)
#     Polynomial(p.data[1] * q.data)
# end
#
# function *(p::Polynomial, q::Polynomial{1})
#     Polynomial(q.data[1] * p.data)
# end


function *(a::Number, q::Polynomial)
    T = promote_type(typeof(a), eltype(q))
    Polynomial{length(q),T}(a*q.data)
end

function +(p::Polynomial{N}, q::Polynomial{M}) where {N,M}
    T = promote_type(eltype(p), eltype(q))
    D = max(N, M)
    x = ntuple(Val(D)) do i
        s = zero(T)
        i <= N && (s += p.data[i])
        i <= M && (s += q.data[i])
        s
    end
    Polynomial{D,T}(SVector{D,T}(x))
end

+(p::Polynomial, a::Number) = p + Polynomial(a)
-(p::Polynomial, a::Number) = p - Polynomial(a)
+(a::Number, p::Polynomial) = p + a
-(a::Number, p::Polynomial) = Polynomial(a) - p


-(p::Polynomial, q::Polynomial) = p + (-one(eltype(q))*q)

#Horner's scheme; the degree of r grows in every step, so the loop is unrolled at compile
#time to keep it type stable (a recursion over Val is not inferred beyond degree 2)
@generated function substitute(p::Polynomial{N}, q::Polynomial) where {N}
    steps = [:(r = q*r + p[$d]) for d in N-2:-1:0]
    return quote
        r = Polynomial(p[$(N-1)])
        $(steps...)
        return r
    end
end

# mutable struct PieceWisePolynomial
# end

# function Bernstein(n,u)
#     basis = zeros(n+1,1)
#     for i = 0:n
#         basis[i+1] = binomial(n,i) * u^i * (1.0-u)^(n-i)
#     end
#     basis
# end
