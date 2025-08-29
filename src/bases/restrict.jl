

"""
    restrict(refspace, element1, element2)

Computes the restriction of a set of local shape functions on `element1` as linear
combinations of the set of local shape functions on `element2`. More precisely `restrict`
returns an `NxM` matrix `P` such that the `i`-th local shape ``g_i`` function on element2
can be written as:

``g_i = sum_{j=1}^{M} P_{ij} f_j``
"""
function restrict end

function restrict(sp::Space, submesh::CompScienceMeshes.AbstractMesh)

    fns = similar(sp.fns)
    S = eltype(eltype(fns))
    for (i,fn) in pairs(sp.fns)
        fns[i] = Vector{S}()
        for sh in fn
            cell = CompScienceMeshes.restrict(submesh, sh.cellid)
            cell == 0 && continue
            push!(fns[i], S(cell, sh.refid, sh.coeff))
        end
    end

    similar(sp, submesh, fns, sp.pos)
end


"""
    extend(submesh, parentmesh)

Extend functions defined on a submesh to their extension by zero on
the parent mesh.
"""
function extend(space::S, submesh::CompScienceMeshes.SubMesh, supermesh) where {S <: Space}

    @assert CompScienceMeshes.parent(submesh) == supermesh

    fns = similar(space.fns)
    for i in 1:numfunctions(space)
        fn = space.fns[i]
        fns[i] = [Shape(CompScienceMeshes.extend(submesh, shape.cellid),
            shape.refid, shape.coeff) for shape in fn]
    end

    return similar(space, supermesh, fns, space.pos)
end


function extend(space::S, supermesh) where {S <: Space}
    extend(space, geometry(space), supermesh)
end


@testitem "extend by zero" begin
    using CompScienceMeshes

    m1 = meshrectangle(1.0, 1.0, 1.0, 3)
    m2 = m1

    m, s = union([m1,m2], Val{:submeshes})


    X1 = raviartthomas(s[1])
    X = BEAST.extend(X1, m)

    @test numfunctions(X1) == 1
    @test numfunctions(X) == 1

    Id = Identity()

    for i in 1:2
        G1 = assemble(Id, X1, X1)
        G2 = assemble(Id, X, X)
        @test G1 ≈ G2
    end
end

function extend(space::F, submesh, supermesh) where {F <: Space}

    S = copy(transpose(CompScienceMeshes.embedding(submesh, supermesh)))

    R = rowvals(S)
    V = nonzeros(S)

    fns = similar(space.fns)
    for i in 1:numfunctions(space)
        fn = space.fns[i]
        fns[i] = map(fn) do sh
            k = nzrange(S, sh.cellid)[1]
            cellid = R[k]
            sign = V[k]
            @assert sign == +1
            Shape(cellid, sh.refid, sh.coeff)
        end
    end

    return similar(space, supermesh, fns, space.pos)
end


@testitem "extend by zero" begin
    using CompScienceMeshes

    m1 = meshrectangle(1.0, 1.0, 1.0, 3)
    m2 = CompScienceMeshes.translate(m1, point(1.0, 0.0, 0.0))

    m = CompScienceMeshes.weld(m1, m2)


    X1 = raviartthomas(m1)
    X = BEAST.extend(X1, m)

    @test numfunctions(X1) == 1
    @test numfunctions(X) == 1

    Id = Identity()

    for i in 1:2
        G1 = assemble(Id, X1, X1)
        G2 = assemble(Id, X, X)
        @test G1 ≈ G2
    end
end