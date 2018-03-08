using StaticArrays

struct Polynomial{N,T}
    data::SVector{N,T}
end
Polynomial(x::T...) where {T} = Polynomial{length(x),T}(SVector{length(x),T}(x))
Polynomial(a::T) where {T<:Number} = Polynomial{1,T}(SVector{1,T}(a))

Base.length(p::Polynomial) = length(p.data)
Base.eltype(p::Polynomial) = eltype(p.data)
Base.getindex(p::Polynomial,i::Int) = p.data[i+1]
Base.setindex!(p::Polynomial,v,i::Int) = (p.data[i+1] = v)

degree(p::Polynomial) = length(p.data)-1

function derive(p)
    x = eltype(p)[d*p[d] for d in 1:degree(p)]
    Polynomial(x...)
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

function *(p::Polynomial,q::Polynomial)
    T = promote_type(eltype(p), eltype(q))
    D = degree(p)+degree(q)+1
    x = zeros(T, D)
    for d in 0 : D
        for k in 0 : D
            k <= degree(p) || continue
            0 <= d-k <= degree(q) || continue
            x[d+1] += p[k] * q[d-k]
        end
    end
    Polynomial(x...)
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

function +(p::Polynomial, q::Polynomial)
    T = promote_type(eltype(p), eltype(q))
    D = max(length(p), length(q))
    x = zeros(T,D)
    for i in 1 : length(p); x[i] += p.data[i]; end
    for i in 1 : length(q); x[i] += q.data[i]; end
    Polynomial(x...)
end

+(p::Polynomial, a::Number) = p + Polynomial(a)
-(p::Polynomial, a::Number) = p - Polynomial(a)
+(a::Number, p::Polynomial) = p + a
-(a::Number, p::Polynomial) = Polynomial(a) - p


-(p::Polynomial, q::Polynomial) = p + (-one(eltype(q))*q)

function substitute(p::Polynomial, q::Polynomial)
    T = promote_type(eltype(p), eltype(q))
    r = Polynomial(p[degree(p)])
    for d = degree(p)-1 : -1 : 0
        r = q*r + p[d]
    end
    return r
end

mutable struct PieceWisePolynomial
end

function Bernstein(n,u)
    basis = zeros(n+1,1)
    for i = 0:n
        basis[i+1] = binomial(n,i) * u^i * (1.0-u)^(n-i)
    end
    basis
end
