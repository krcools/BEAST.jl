import WiltonInts84

mutable struct WiltonSEStrategy{P,Q} <: SingularityExtractionStrategy
    outer_quad_points::P
    regularpart_quadrule::DoubleQuadStrategy{P,Q}
end

function innerintegrals!(op::MWSingleLayer3DSng, p, g, f, t, s, z,
        strat::WiltonSEStrategy, dx)

    γ = op.gamma
    x = cartesian(p)
    n = cross(s[1]-s[3],s[2]-s[3])
    n /= norm(n)
    ρ = x - ((x-s[1]) ⋅ n) * n

    scal, vec = WiltonInts84.wiltonints(s[1], s[2], s[3], x, Val{1})

    # \int \frac{1}{4 \pi R}
    ∫G = (scal[2] - γ*scal[3] + 0.5*γ^2*scal[4]) / (4π)
    # \int \frac{y}{4 \pi R}
    ∫Gy = SVector((
        (vec[2][1] + scal[2]*ρ[1] - γ*(vec[3][1]+scal[3]*ρ[1]) + 0.5*γ^2*(vec[4][1]+scal[4]*ρ[1]))/(4π),
        (vec[2][2] + scal[2]*ρ[2] - γ*(vec[3][2]+scal[3]*ρ[2]) + 0.5*γ^2*(vec[4][2]+scal[4]*ρ[2]))/(4π),
        (vec[2][3] + scal[2]*ρ[3] - γ*(vec[3][3]+scal[3]*ρ[3]) + 0.5*γ^2*(vec[4][3]+scal[4]*ρ[3]))/(4π),
    ))

    c₁ = op.α
    c₂ = op.β

    α = 1 / volume(t) / volume(s) / 4
    for i in 1 : numfunctions(g)
        a = t[i]
        g = x - a
        dg = 2

        for j in 1 : numfunctions(f)
            b = s[j]

            ∫Gf = SVector(∫Gy[1]-∫G*b[1], ∫Gy[2]-∫G*b[2], ∫Gy[3]-∫G*b[3])
            ∫Gdf = 2 * ∫G
            dg∫Gf = g[1]*∫Gf[1] + g[2]*∫Gf[2] + g[3]*∫Gf[3]
            z[i,j] += ( α*c₁*dg∫Gf + α*c₂*dg*∫Gdf ) * dx

        end # next j
    end #

end


function innerintegrals!(op::MWDoubleLayer3DSng, p, g, f, t, s, z, strat::WiltonSEStrategy, dx)

    γ = op.gamma
    x = cartesian(p)
    n = cross(s[1]-s[3],s[2]-s[3])
    n /= norm(n)
    ρ = x - ((x-s[1]) ⋅ n) * n

    scal, vec, grad = WiltonInts84.wiltonints(s[1], s[2], s[3], x, Val{1})

    # \int \nabla G_s with G_s = \nabla (1/R + 0.5*γ^2*R) / (4\pi)
    ∫∇G = (-grad[1] - 0.5*γ^2*grad[3]) / (4π)

    α = 1 / volume(t) / volume(s) / 4
    for i in 1 : numfunctions(g)
        a = t[i]
        g = (x - a)

        for j in 1 : numfunctions(f)
            b = s[j]

            z[i,j] += ( α * ( (x-b) × g ) ⋅ ∫∇G ) * dx

        end # next j
    end #

end
