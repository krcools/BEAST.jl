import BEAST; BE = BEAST;
using CompScienceMeshes
using Test

T = Float64
P = SVector{3,T}
tol = eps(T)*10^3

struct IntegrationDomain{T,P}
    signed_height::T
    dist_to_plane::T
    center::P
    plane_center::P
    unit_normal::P
    segments::Vector{SVector{2,P}}
end


function IntegrationDomain(p1, p2, p3, center)

    t1 = p2 - p1
    t2 = p3 - p1

    unit_normal = t1 × t2
    unit_normal /= norm(unit_normal)
    signed_height = dot(center-p1, unit_normal)
    dist_to_plane = abs(signed_height)
    plane_center = center - signed_height * unit_normal

    s1 = SVector((p1, p2,))
    s2 = SVector((p2, p3,))
    s3 = SVector((p3, p1,))
    segments = [s1, s2, s3]

    #T = eltype(p1)
    IntegrationDomain(
        signed_height,
        dist_to_plane,
        center,
        plane_center,
        unit_normal,
        segments)

end

wiltonints(a,b,c,center) = wiltonints(IntegrationDomain(a,b,c,center))

# type WiltonIntsWS{T}
#     segments_through_c::Vector{T}
# end

# type WiltonIntsRV{T,N,M}
#     scalars::SVector{N,T}
#     gradients::SVector{M,T}
# end

# This function implements the Wilton algorithm for Integration of
# linear potentials on triangular domains
function wiltonints(domain::IntegrationDomain{T,Q}, tol=eps(T)*10^3) where {T,Q}

    dim = 3

    d = domain.signed_height
    D = domain.dist_to_plane
    c = domain.plane_center
    n = domain.unit_normal

    num_segments = length(domain.segments)

    # outside the domain of integration, or whether it is on the boundary.
    # keep track of whether the projected singularity is inside or
    # If it is on the boundary, make a note of on which part of the boundary
    # it is located.
    c_inside = true
    num_segments_through_c = 0
    segments_through_c = zeros(Int, num_segments)

    s = zeros(T,6) # TODO alloc
    v = zeros(Q, 6) # TODO alloc

    for i in 1 : num_segments

        a = domain.segments[i][1]
        b = domain.segments[i][2]

        ab = b - a     # segment tangent
        l = norm(ab)   # segment length
        t = ab / l     # unit tangent
        m = cross(t,n) # unit binormal

        ca = a - c
        cb = b - c

        # coordinates of a and b w.r.t. the projection
        # of c on the line [a,b]
        xa = dot(ca,t)
        xb = dot(cb,t)

        # signed and absolute distance of plane center to line
        p = dot(ca,m)
        P = abs(p)

        if p < 0
          c_inside = false
        end

        w2 = p^2 + d^2
        ra = sqrt(xa^2 + w2)
        rb = sqrt(xb^2 + w2)
        w = sqrt(w2)

        # if the plane center is NOT on the line through this segment
        # use the general formula; no singularities are encountered in
        # the computation
        if P > tol
            lg = log((xb+rb) / (xa+ra))
            at = atan(d*xb/(p*rb)) - atan(d*xa/(p*ra))

            s[1] += -at
            s[2] += p*lg + d*at
            s[3] += 1/2 * p * (xb-xa)
            pt1 = 1/6*p*(xb*rb-xa*ra)
            pt2 = (1/2*d*d+1/6*p*p)*p*lg
            pt3 = 1/3*d*d*d*at
            s[4] += pt1 + pt2 + pt3
            pt1 = 1/2*xb*p*(1/6*xb*xb+1/2*p*p+d*d)
            pt2 = 1/2*xa*p*(1/6*xa*xa+1/2*p*p+d*d)
            s[5] += pt1 - pt2
            pt1b = 1/20*p*xb*rb*(xb*xb+1/2*(9*d*d+5*p*p))
            pt1a = 1/20*p*xa*ra*(xa*xa+1/2*(9*d*d+5*p*p))
            pt2 = 1/40*p*(15*d*d*d*d+10*d*d*p*p+3*p*p*p*p)*lg + 0.2*d*d*d*d*d*at;
            s[6] += (pt1b - pt1a) + pt2

            v[1] += -lg * m
            v[2] += 0.5 *((xb*rb-xa*ra) + w*w*lg) * m
            pt1b = 0.5 * xb * (1/3*xb*xb + p*p)
            pt1a = 0.5 * xa * (1/3*xa*xa + p*p)
            v[3] += (pt1b - pt1a) * m
            pt1 = 1/12*(xb*rb*rb*rb-xa*ra*ra*ra) + 1/8*w*w*(xb*rb-xa*ra)
            pt2 = 1/8*w*w*w*w*lg
            v[4] += (pt1 + pt2) * m
            pt1b = 1/20*xb*xb*xb*xb*xb + 1/6*p*p*xb*xb*xb + 0.25*p*p*p*p*xb + 0.5*d*d*(1/3*xb*xb*xb + p*p*xb)
            pt1a = 1/20*xa*xa*xa*xa*xa + 1/6*p*p*xa*xa*xa + 0.25*p*p*p*p*xa + 0.5*d*d*(1/3*xa*xa*xa + p*p*xa)
            v[5] += (pt1b - pt1a) * m


        else # plane center lies on the line ab

            if xa*xb < tol^2
                # projection c on ab lies in between a and b
                num_segments_through_c += 1
                segments_through_c[num_segments_through_c] = i
            end

            if D > tol
                # singularity not in the plane of the triangle and thus not on ab

                lg = log((xb+rb)/(xa+ra))

                v[1] += -lg * m
                pt1 = 0.5*(xb*rb-xa*ra)
                pt2 = 0.5*w*w*lg
                v[2] += (pt1 + pt2) * m
                pt1b = 1/6*xb*xb*xb
                pt1a = 1/6*xa*xa*xa
                v[3] += (pt1b - pt1a) * m
                pt1 = 1/12*(xb*rb*rb*rb-xa*ra*ra*ra) + 1/8*w*w*(xb*rb-xa*ra)
                pt2 = 1/8*w*w*w*w*lg
                v[4] += (pt1 + pt2) * m
                pt1b = 1/20*xb*xb*xb*xb*xb + 1/6*p*p*xb*xb*xb + 0.25*p*p*p*p*xb + 0.5*d*d*(1/3*xb*xb*xb + p*p*xb)
                pt1a = 1/20*xa*xa*xa*xa*xa + 1/6*p*p*xa*xa*xa + 0.25*p*p*p*p*xa + 0.5*d*d*(1/3*xa*xa*xa + p*p*xa)
                v[5] += (pt1b - pt1a) * m

            else # D <= tol

                # the projected singularity lies on ab. It is assumed the singularity
                # itself does not lie on [a,b] however, since this would lead to some
                # of the integrals to become infinite.

                lg = (xb > 0) ? log(xb/xa) : log(xa/xb)

                v[1] += -lg * m
                v[2] += 0.5*(xb*rb-xa*ra) * m
                pt1b = 1/6*xb*xb*xb
                pt1a = 1/6*xa*xa*xa
                v[3] += (pt1b - pt1a) * m
                pt1 = 1/12*(xb*rb*rb*rb-xa*ra*ra*ra) + 1/8*w*w*(xb*rb-xa*ra)
                pt2 = 1/8*w*w*w*w*lg
                v[4] += (pt1 + pt2) * m
                pt1b = 1/20*xb*xb*xb*xb*xb + 1/6*p*p*xb*xb*xb + 0.25*p*p*p*p*xb + 0.5*d*d*(1/3*xb*xb*xb + p*p*xb)
                pt1a = 1/20*xa*xa*xa*xa*xa + 1/6*p*p*xa*xa*xa + 0.25*p*p*p*p*xa + 0.5*d*d*(1/3*xa*xa*xa + p*p*xa)
                v[5] += (pt1b - pt1a) * m

            end # if D < tol

        end # if p < tol

    end # next segment

    # contributions from intersection of epsilon cirlce around c with triangle
    sgn = D < tol ? zero(d) : sign(d)
    alpha = zero(d)
    if num_segments_through_c == 0 && c_inside
        alpha = 2pi
    elseif num_segments_through_c == 1
        alpha = pi
    elseif num_segments_through_c == 2
        i = segments_through_c[1]
        t1 = domain.segments[i][2] - domain.segments[i][1]
        i = segments_through_c[2]
        t2 = domain.segments[i][2] - domain.segments[i][1]
        cosalpha = dot(-t1,t2)/(norm(t1)*norm(t2))
        alpha = acos(clamp(cosalpha,-1,+1))
    elseif num_segments_through_c > 2
        error("num_intersections must be either 0, 1, or 2.")
    end

    s[1] += alpha * sgn
    s[2] += -alpha * D
    s[4] += -1/3 * alpha * D*D*D
    s[6] += -1/5 * alpha * D*D*D*D*D

    # build gradients
    g = Vector{Q}(undef,4)
    g[1] = -1( v[1] - s[1] * n )
    g[2] = +1( v[2] - d*s[2]*n )
    g[3] = +2( v[3] - d*s[3]*n )
    g[4] = +3( v[4] - d*s[4]*n )
    return s, v, g

end

# s1 = int d/R^3
# s2 = int 1/R
# s3 = int 1
# s4 = int R
# s5 = int R^2
# s6 = int R^3
# v1 = int (rho'-rho)/R^3
# v2 = int (rho'-rho)/R
# v3 = int (rho'-rho)
# v4 = int (rho'-rho)*R
# v5 = int (rho'-rho)*R
# g1 = int grad' (1/R)
# g2 = int grad' R
# g3 = int grad' R^2
# g4 = int grad' R^3
function withquadrules(triangle, r, n)

	u, w = BE.trgauss(n)

    n = (triangle[1] - triangle[3]) × (triangle[2] - triangle[3])
    n /= norm(n)
    d = (r - triangle[1]) ⋅ n
    ρ = r - d*n

    scalar = zeros(T,6)
    vector = zeros(P,6)
    gradgr = zeros(P,4)

	m = simplex(triangle[1], triangle[2], triangle[3])

    for i in 1 : length(w)
        q = neighborhood(m, u[:,i])
        s = cartesian(q)
        dq = w[i] * jacobian(q)

        R = norm(s - r)
        scalar[1] += d / R^3 * dq
        scalar[2] += 1 / R * dq
        scalar[3] += 1 * dq
        scalar[4] += R * dq
        scalar[5] += R^2 * dq
        scalar[6] += R^3 * dq
        vector[1] += (s - ρ) / R^3 * dq
        vector[2] += (s - ρ) / R * dq
        vector[3] += (s - ρ) * dq
        vector[4] += (s - ρ) * R * dq
        vector[5] += (s - ρ) * R^2 * dq
    end

    gradgr[1] = -1( vector[1] - scalar[1] * n )
    gradgr[2] = +1( vector[2] - d*scalar[2]*n )
    gradgr[3] = +2( vector[3] - d*scalar[3]*n )
    gradgr[4] = +3( vector[4] - d*scalar[4]*n )

    return scalar, vector, gradgr
end

tr = [
	point(0.0, 0.0, 0.0),
	point(2.0, 0.0, 0.0),
	point(0.0, 3.0, 0.0)]

pts = [
    (tr[1]+2*tr[2]-1*tr[3]),
    (tr[1]+tr[2])/2 + point(0.0, 0.0, 0.3),
    point(3.0, 0.0, 0.0),
    (tr[1]+tr[2]+tr[3])/3 + point(0.0, 0.0, 0.3),
    tr[2] + point(0.0, 0.0, 0.3),
    (tr[1]+tr[2]+tr[3])/3 + point(0.0, 0.0, 100.0)
]

for _p in pts
    s1, v1, g1 = wiltonints(tr[1],tr[2],tr[3],_p)
    s2, v2, g2 = withquadrules(tr,_p,13)

    for i in length(s1)
        @test norm(s1[i]-s2[i]) < 1.0e-5
    end

    for i in length(v1)
        @test norm(v1[i]-v2[i]) < 1.0e-5
    end

    for i in 1:length(g1)
        for _j in 1:3
            @test norm(g1[i][_j]-g2[i][_j]) < 1.0e-5
        end
    end
end


Γ = readmesh(joinpath(dirname(@__FILE__),"assets","sphere2.in"))
nc = numcells(Γ)
t = chart(Γ, first(cells(Γ)))
s = chart(Γ, last(cells(Γ)))

X = BEAST.raviartthomas(Γ)
x = BEAST.refspace(X)

κ = 1.0
op = BEAST.MWSingleLayer3D(κ)

n = BE.numfunctions(x)
z1 = zeros(ComplexF64, n, n)
z2 = zeros(ComplexF64, n, n)

BE = BEAST

tqd = BE.quadpoints(x, [t], (12,))
bqd = BE.quadpoints(x, [s], (13,))

DQ_strategy = BE.DoubleQuadStrategy(tqd[1,1], bqd[1,1])
BEAST.momintegrals!(op, x, x, t, s, z1, DQ_strategy)

SE_strategy = BE.WiltonSEStrategy(
  tqd[1,1],
  BE.DoubleQuadStrategy(
	tqd[1,1],
	bqd[1,1],
  ),
)
BEAST.momintegrals!(op, x, x, t, s, z2, SE_strategy)

@test norm(z1-z2) < 1.0e-7
