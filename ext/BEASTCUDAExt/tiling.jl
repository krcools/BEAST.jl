abstract type Tiling end

struct EqualTiling <: Tiling
    num_tiles::Int

    function EqualTiling(num_tiles)
        num_tiles > 0 || throw(ArgumentError("number of tiles must be positive"))
        return new(num_tiles)
    end
end

struct WorksizeTiling <: Tiling
    tile_size::Int

    function WorksizeTiling(tile_size)
        tile_size > 0 || throw(ArgumentError("tile size must be positive"))
        return new(tile_size)
    end
end

function tileindices(workload_size::Int, tiling::EqualTiling)
    num_tiles = min(tiling.num_tiles, workload_size)
    boundaries = round.(Int, range(0, stop=workload_size, length=num_tiles + 1))
    return [collect(boundaries[i]+1:boundaries[i+1]) for i in 1:num_tiles]
end

function tileindices(workload_size::Int, tiling::WorksizeTiling)
    return [collect(start:min(start + tiling.tile_size - 1, workload_size))
        for start in 1:tiling.tile_size:workload_size]
end

struct TilingStrategy{A<:Tiling,B<:Tiling}
    test::A
    trial::B
end

function Base.getindex(strategy::TilingStrategy, i::Int)
    i == 1 && return strategy.test
    i == 2 && return strategy.trial
    throw(BoundsError(strategy, i))
end
