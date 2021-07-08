# Utils to lift LinearMaps to LinearMaps acting on larger vector spaces
import LinearMaps
import LinearAlgebra

struct LiftedMap{T,TI,TJ} <: LinearMap{T}
    A::LinearMap{T}
    I::TI
    J::TJ
end

function LinearAlgebra.mul!(y::AbstractVector, L::LiftedMap, x::AbstractVector, α::Number, β::Number)

    yI = view(y, L.I)
    xJ = view(x, L.J)
    AIJ = L.A
    LinearAlgebra.mul!(yI, AIJ, xJ, α, β)
end