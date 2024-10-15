struct NCrossBDMRefSpace{T} <: RefSpace{T} end

function (f::NCrossBDMRefSpace{T})(p) where T

    u,v = parametric(p)
    n = normal(p)
    tu = tangents(p,1)
    tv = tangents(p,2)

    j = jacobian(p)
    d = 1/j

    return @SVector[
        (value= n × (-v*tu+v*tv)/j, curl=d),
        (value= n × ((u+v-1)*tu)/j,   curl=d),
        (value= n × ((u+v-1)*tv) /j,   curl=d),
        (value= n × (u*tu-u*tv)/j,  curl=d),
        (value= n × (u*tu)/j,         curl=d),
        (value= n × (v*tv)/j,         curl=d),]
end

const _vert_perms_ncrossbdm = [
    (1,2,3),
    (2,3,1),
    (3,1,2),
    (2,1,3),
    (1,3,2),
    (3,2,1),
]
const _dof_perms_ncrossbdm = [
    (1,2,3,4,5,6),
    (5,6,1,2,3,4),
    (3,4,5,6,1,2),
    (4,3,2,1,6,5),
    (2,1,6,5,4,3),
    (6,5,4,3,2,1),
]

function dof_permutation(::NCrossBDMRefSpace, vert_permutation)
    i = findfirst(==(tuple(vert_permutation...)), _vert_perms_ncrossbdm)
    return _dof_perms_ncrossbdm[i]
end