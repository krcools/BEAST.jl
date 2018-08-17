struct BDMRefSpace{T} <: RefSpace{T,6} end

function (f::BDMRefSpace)(p)

    u,v = parametric(p)

    tu = tangents(p,1)
    tv = tangents(p,2)

    j = jacobian(p)
    d = 1/j

    return @SVector[
        (value=(u+v-1)*tu/j,   divergence=d),
        (value=(-v*tu+v*tv)/j, divergence=d),
        (value=(u+v-1)*tv/j,   divergence=d),
        (value=(u*tu-u*tv)/j,  divergence=d),
        (value=u*tu/j,         divergence=d),
        (value=v*tv/j,         divergence=d),]
end

# Default quadrule for use with these shape functions. Overload if
# better options are available for the specific operator kernels used.

function quadrule(op, g::BDMRefSpace, f::BDMRefSpace, i, τ, j, σ, qd)
    # defines coincidence of points
    dtol = 1.0e3 * eps(eltype(eltype(τ.vertices)))

    # decides on whether to use singularity extraction
    xtol = 0.2
    k = norm(op.gamma)

    hits = 0
    xmin = xtol
    for t in τ.vertices
      for s in σ.vertices
        d = norm(t-s)
        xmin = min(xmin, k*d)
        if d < dtol
          hits +=1
          break
        end
      end
    end

  hits == 3   && return SauterSchwabQuadrature.CommonFace(qd.gausslegendre[3])
  hits == 2   && return SauterSchwabQuadrature.CommonEdge(qd.gausslegendre[2])
  hits == 1   && return SauterSchwabQuadrature.CommonVertex(qd.gausslegendre[1])
  xmin < xtol && return DoubleQuadStrategy(
      qd.tpoints[2,i],
      qd.bpoints[2,j],)
  return DoubleQuadStrategy(
    qd.tpoints[1,i],
    qd.bpoints[1,j],)
end
