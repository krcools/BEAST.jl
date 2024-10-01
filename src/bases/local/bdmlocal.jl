struct BDMRefSpace{T} <: RefSpace{T} end

function (f::BDMRefSpace)(p)

    u,v = parametric(p)

    tu = tangents(p,1)
    tv = tangents(p,2)

    j = jacobian(p)
    d = 1/j

    return @SVector[
        (value=(-v*tu+v*tv)/j, divergence=d),
        (value=(u+v-1)*tu/j,   divergence=d),
        (value=(u+v-1)*tv/j,   divergence=d),
        (value=(u*tu-u*tv)/j,  divergence=d),
        (value=u*tu/j,         divergence=d),
        (value=v*tv/j,         divergence=d),]
end

divergence(ref::BDMRefSpace, sh, el) = Shape(sh.cellid, 1, sh.coeff/(2*volume(el)))

const _vert_perms_bdm = [
    (1,2,3),
    (2,3,1),
    (3,1,2),
    (2,1,3),
    (1,3,2),
    (3,2,1),
]
const _dof_perms_bdm = [
    (1,2,3,4,5,6),
    (5,6,1,2,3,4),
    (3,4,5,6,1,2),
    (4,3,2,1,6,5),
    (2,1,6,5,4,3),
    (6,5,4,3,2,1),
]

function dof_permutation(::BDMRefSpace, vert_permutation)
    i = findfirst(==(tuple(vert_permutation...)), _vert_perms_bdm)
    return _dof_perms_bdm[i]
end

function dof_perm_matrix(::BDMRefSpace, vert_permutation)
    i = findfirst(==(tuple(vert_permutation...)), _vert_perms_bdm)
    return _dof_bdmperm_matrix[i]
end

dimtype(::BDMRefSpace, ::CompScienceMeshes.Simplex{U,2}) where {U} = Val{6}

const _dof_bdmperm_matrix = [
    
    @SMatrix[1 0 0 0 0 0;  #1. {1,2,3}
                0 1 0 0 0 0;
                0 0 1 0 0 0;
                0 0 0 1 0 0;
                0 0 0 0 1 0;
                0 0 0 0 0 1],

    @SMatrix[0 0 0 0 1 0;  #2. {2,3,1}
                0 0 0 0 0 1;
                1 0 0 0 0 0;
                0 1 0 0 0 0;
                0 0 1 0 0 0;
                0 0 0 1 0 0],

    @SMatrix[0 0 1 0 0 0;  #3. {3,1,2}
                0 0 0 1 0 0;
                0 0 0 0 1 0;
                0 0 0 0 0 1;
                1 0 0 0 0 0;
                0 1 0 0 0 0],

    @SMatrix[0 0 0 1 0 0;  #4. {2,1,3}
                0 0 1 0 0 0;
                0 1 0 0 0 0;
                1 0 0 0 0 0;
                0 0 0 0 0 1;
                0 0 0 0 1 0],

    @SMatrix[0 1 0 0 0 0;  #5. {1,3,2}
                1 0 0 0 0 0;
                0 0 0 0 0 1;
                0 0 0 0 1 0;
                0 0 0 1 0 0;
                0 0 1 0 0 0],

    @SMatrix[0 0 0 0 0 1;  #6. {3,2,1}
                0 0 0 0 1 0;
                0 0 0 1 0 0;
                0 0 1 0 0 0;
                0 1 0 0 0 0;
                1 0 0 0 0 0]

]