struct RTRefSpace{T} <: DivRefSpace{T} end

# valuetype(ref::RTRefSpace{T}, charttype) where {T} = SVector{3,Tuple{SVector{universedimension(charttype),T},T}}
function valuetype(ref::RTRefSpace{T}, charttype::Type) where {T}
    SVector{universedimension(charttype),T}
end

function (ϕ::RTRefSpace)(mp)

    u, v = parametric(mp)
    j = jacobian(mp)

    tu = tangents(mp,1)
    tv = tangents(mp,2)

    inv_j = 1/j
    d = 2 * inv_j

    u_tu = u*tu
    v_tv = v*tv

    return SVector((
        (value=(u_tu-tu + v_tv    )*inv_j, divergence=d),
        (value=(u_tu     + v_tv-tv)*inv_j, divergence=d),
        (value=(u_tu     + v_tv    )*inv_j, divergence=d)
    ))
end

numfunctions(x::RTRefSpace, dom::CompScienceMeshes.ReferenceSimplex{2}) = 3

divergence(ref::RTRefSpace, sh, el) = [Shape(sh.cellid, 1, sh.coeff/volume(el))]

"""
    ntrace(refspace, element, localindex, face)

Compute the normal trace of all local shape functions on `elements` belonging to
`refspace` on `face`. This function returns a matrix expressing the traces of local
shape functions in `refspace` as linear combinations of functions in the local
trace space. Cf. `restrict`. `localindex` is the index of `face` in the enumeration
of faces of `elements`. In many special cases knowing this index allows for highly
optimised implementations.
"""
function ntrace(x::RTRefSpace, el, q, fc)
    t = zeros(scalartype(x),1,3)
    t[q] = 1 / volume(fc)
    return t
end


const _vert_perms_rt = [
    (1,2,3),
    (2,3,1),
    (3,1,2),
    (2,1,3),
    (1,3,2),
    (3,2,1),
]
const _dof_perms_rt = [
    (1,2,3),
    (3,1,2),
    (2,3,1),
    (2,1,3),
    (1,3,2),
    (3,2,1),
]

function dof_permutation(::RTRefSpace, vert_permutation)
    i = something(findfirst(==(tuple(vert_permutation...)), _vert_perms_rt),0)
    @assert i != 0
    return _dof_perms_rt[i]
end

function dof_perm_matrix(::RTRefSpace, vert_permutation)
    i = findfirst(==(tuple(vert_permutation...)), _vert_perms_rt)
    @assert i != nothing
    return _dof_rtperm_matrix[i]
end

"""
    interpolate(interpolant::RefSpace, chart1, interpolee::RefSpace, chart2)

Computes by interpolation approximations of the local shape functions for
`interpolee` on `chart2` in terms of the local shape functions for `interpolant`
on `chart1`. The returned value is a matrix `Q` such that

```math
\\phi_i \\approx \\sum_j Q_{ij} \\psi_j
```

with ``\\phi_i`` the i-th local shape function for `interpolee` and ``\\psi_j`` the
j-th local shape function for `interpolant`.
"""
function interpolate(interpolant::RefSpace, chart1, interpolee::RefSpace, chart2)
    function fields(p)
        x = cartesian(p)
        v = carttobary(chart2, x)
        r = neighborhood(chart2, v)
        fieldvals = [f.value for f in interpolee(r)]
    end

    interpolate(fields, interpolant, chart1)
end


function interpolate(interpolant::RefSpace, chart1, interpolee::RefSpace, chart2, ch1toch2)
    function fields(p1)
        u1 = parametric(p1)
        u2 = cartesian(ch1toch2, u1)
        p2 = neighborhood(chart2, u2)
        fieldvals = [f.value for f in interpolee(p2)]
    end

    interpolate(fields, interpolant, chart1)
end


function interpolate(fields, interpolant::RTRefSpace, chart)
    Q = map(faces(chart)) do face
        p = center(face)
        x = cartesian(p)
        u = carttobary(chart, x)
        q = neighborhood(chart, u)
        n = normal(q)

        # minus because in CSM the tangent points towards vertex[1]
        t = -tangents(p,1)
        m = cross(t,n)

        fieldvals = fields(q)
        q = [dot(fv,m) for fv in fieldvals]
    end

    return hcat(Q...)
end


function restrict(ϕ::RefSpace, dom1, dom2)
    interpolate(ϕ, dom2, ϕ, dom1)
end

function restrict(ϕ::RefSpace, dom1, dom2, dom2todom1)
    interpolate(ϕ, dom2, ϕ, dom1, dom2todom1)
end

const _dof_rtperm_matrix = [
    @SMatrix[1 0 0;         # 1. {1,2,3}
                0 1 0;
                0 0 1],

    @SMatrix[0 0 1;         # 2. {2,3,1}
                1 0 0;
                0 1 0],

    @SMatrix[0 1 0;         # 3. {3,1,2}
                0 0 1;
                1 0 0],

    @SMatrix[0 1 0;         # 4. {2,1,3}
                1 0 0;
                0 0 1],

    @SMatrix[1 0 0;         # 5. {1,3,2}
                0 0 1;
                0 1 0],

    @SMatrix[0 0 1;         # 6. {3,2,1}
                0 1 0;
                1 0 0]
]


# Support for zeroth order elements on quadrilaterals
