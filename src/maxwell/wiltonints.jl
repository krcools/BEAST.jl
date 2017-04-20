export wiltonints
export WiltonSEStrategy


immutable IntegrationDomain{T,P}
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

type WiltonSEStrategy{P,Q} <: SingularityExtractionStrategy
    outer_quad_points::P
    regularpart_quadrule::DoubleQuadStrategy{P,Q}
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
# g1 = int grad' (1/R) = -int (r' - r) / R^3
# g2 = int grad' R     = int (r' -r) / R
# g3 = int grad' R^2   = 2 * int (r' - r)
# g4 = int grad' R^3   = 3 * int (r' - r) * R
wiltonints(a,b,c,center) = wiltonints(IntegrationDomain(a,b,c,center))

type WiltonIntsWS{T}
    segments_through_c::Vector{T}
end

type WiltonIntsRV{T,N,M}
    scalars::SVector{N,T}
    gradients::SVector{M,T}
end

# This function implements the Wilton algorithm for Integration of
# linear potentials on triangular domains
function wiltonints{T,Q}(domain::IntegrationDomain{T,Q}, tol=eps(T)*10^3)

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
            s[3] += 1./2. * p * (xb-xa)
            pt1 = 1./6.*p*(xb*rb-xa*ra)
            pt2 = (1./2.*d*d+1./6.*p*p)*p*lg
            pt3 = 1./3.*d*d*d*at
            s[4] += pt1 + pt2 + pt3
            pt1 = 1./2.*xb*p*(1./6.*xb*xb+1./2.*p*p+d*d)
            pt2 = 1./2.*xa*p*(1./6.*xa*xa+1./2.*p*p+d*d)
            s[5] += pt1 - pt2
            pt1b = 1./20.*p*xb*rb*(xb*xb+1./2.*(9*d*d+5*p*p))
            pt1a = 1./20.*p*xa*ra*(xa*xa+1./2.*(9*d*d+5*p*p))
            pt2 = 1./40.*p*(15*d*d*d*d+10*d*d*p*p+3*p*p*p*p)*lg + 0.2*d*d*d*d*d*at;
            s[6] += (pt1b - pt1a) + pt2

            v[1] += -lg * m
            v[2] += 0.5 *((xb*rb-xa*ra) + w*w*lg) * m
            pt1b = 0.5 * xb * (1./3.*xb*xb + p*p)
            pt1a = 0.5 * xa * (1./3.*xa*xa + p*p)
            v[3] += (pt1b - pt1a) * m
            pt1 = 1./12.*(xb*rb*rb*rb-xa*ra*ra*ra) + 1./8.*w*w*(xb*rb-xa*ra)
            pt2 = 1./8.*w*w*w*w*lg
            v[4] += (pt1 + pt2) * m
            pt1b = 1./20.*xb*xb*xb*xb*xb + 1./6.*p*p*xb*xb*xb + 0.25*p*p*p*p*xb + 0.5*d*d*(1./3.*xb*xb*xb + p*p*xb)
            pt1a = 1./20.*xa*xa*xa*xa*xa + 1./6.*p*p*xa*xa*xa + 0.25*p*p*p*p*xa + 0.5*d*d*(1./3.*xa*xa*xa + p*p*xa)
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
                pt1b = 1./6.*xb*xb*xb
                pt1a = 1./6.*xa*xa*xa
                v[3] += (pt1b - pt1a) * m
                pt1 = 1./12.*(xb*rb*rb*rb-xa*ra*ra*ra) + 1./8.*w*w*(xb*rb-xa*ra)
                pt2 = 1./8.*w*w*w*w*lg
                v[4] += (pt1 + pt2) * m
                pt1b = 1./20.*xb*xb*xb*xb*xb + 1./6.*p*p*xb*xb*xb + 0.25*p*p*p*p*xb + 0.5*d*d*(1./3.*xb*xb*xb + p*p*xb)
                pt1a = 1./20.*xa*xa*xa*xa*xa + 1./6.*p*p*xa*xa*xa + 0.25*p*p*p*p*xa + 0.5*d*d*(1./3.*xa*xa*xa + p*p*xa)
                v[5] += (pt1b - pt1a) * m

            else # D <= tol

                # the projected singularity lies on ab. It is assumed the singularity
                # itself does not lie on [a,b] however, since this would lead to some
                # of the integrals to become infinite.

                lg = (xb > 0) ? log(xb/xa) : log(xa/xb)

                v[1] += -lg * m
                v[2] += 0.5*(xb*rb-xa*ra) * m
                pt1b = 1./6.*xb*xb*xb
                pt1a = 1./6.*xa*xa*xa
                v[3] += (pt1b - pt1a) * m
                pt1 = 1./12.*(xb*rb*rb*rb-xa*ra*ra*ra) + 1./8.*w*w*(xb*rb-xa*ra)
                pt2 = 1./8.*w*w*w*w*lg
                v[4] += (pt1 + pt2) * m
                pt1b = 1./20.*xb*xb*xb*xb*xb + 1./6.*p*p*xb*xb*xb + 0.25*p*p*p*p*xb + 0.5*d*d*(1./3.*xb*xb*xb + p*p*xb)
                pt1a = 1./20.*xa*xa*xa*xa*xa + 1./6.*p*p*xa*xa*xa + 0.25*p*p*p*p*xa + 0.5*d*d*(1./3.*xa*xa*xa + p*p*xa)
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
    s[4] += -1./3. * alpha * D*D*D
    s[6] += -1./5. * alpha * D*D*D*D*D

    # build gradients
    g = Vector{Q}(4)
    g[1] = -1( v[1] - s[1] * n )
    g[2] = +1( v[2] - d*s[2]*n )
    g[3] = +2( v[3] - d*s[3]*n )
    g[4] = +3( v[4] - d*s[4]*n )
    return s, v, g

end
