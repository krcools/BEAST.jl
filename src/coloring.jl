"""
    function conflicts(
        space::Space;
        addata=assemblydata(space),
        kwargs...,
    )

Computes conflict indices for a Space. Two elements are in conflict if they are both
part of the support of the same basis function.

# Arguments

  - `space::Space`: The space.
  - `addata=assemblydata(space)`: The assembly data for the space (default: computed using `assemblydata`).
  - `kwargs...`: Additional keyword arguments.

# Returns

A tuple containing:

  - `eachindex(elements)`: The indices of the elements.
  - `ConflictFunctor(conflictindices)`: A functor that maps element indices to conflict indices.
  - `Base.OneTo(numfunctions(space))`: The indices of the conflicts.
"""
function GraphsColoring.conflicts(
    space::Space;
    addata=assemblydata(space),
    kwargs...,
)
    elements, ad, _ = addata

    conflictindices = Vector{Int}[Int[] for _ in eachindex(elements)]

    for elementid in eachindex(elements)
        for i in 1:length(ad[elementid])
            for (functionid, a) in ad[elementid, i]
                iszero(a) && continue
                push!(conflictindices[elementid], functionid)
            end
        end
    end

    for i in eachindex(conflictindices)
        conflictindices[i] = unique(conflictindices[i])
    end

    return eachindex(elements),
    GraphsColoring.ConflictFunctor(conflictindices),
    Base.OneTo(numfunctions(space))
end

function color(basis, ::Any; addata=assemblydata(basis), kwargs...)
    return GraphsColoring.color(GraphsColoring.conflictmatrix(basis; addata=addata, kwargs...), WorkstreamDSATUR).colors
end

# for serial scheduling, just return all elements in one color
function color(basis, ::SerialScheduler; addata=assemblydata(basis), kwargs...)
    return [collect(1:length(first(addata)))]
end
