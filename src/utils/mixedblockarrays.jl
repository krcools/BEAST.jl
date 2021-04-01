struct ZeroBlockInitializer end
const zero_blocks = ZeroBlockInitializer()

import Base.Cartesian: @nloops, @ntuple
import BlockArrays: BlockArray, undef_blocks, blockaxes
import FillArrays: Zeros, Fill

"""
Initialise a BlockArray where each block can have a different type.

    BlockArray{T}(zero_blocks, blocksize...)
"""
@generated function BlockArray{T}(::ZeroBlockInitializer, blocksizes::Vararg{AbstractVector{Int},N}) where {T,N}
    return quote
        B = BlockArray{T,N,Array{AbstractArray{T,N},N}}(undef_blocks, blocksizes...)
        axs = axes(B)
        @nloops $N block dim->blockaxes(B,dim) begin
            block_indices = @ntuple $N block
            indices = getindex.(axs, block_indices)
            block_size = length.(indices)
            B[block_indices...] = Zeros{T}(block_size...)
        end
        return B
    end
end

@generated function BlockArray{T}(::ZeroBlockInitializer, blocksizes::Vararg{AbstractVector{Int},N}) where {T <: Array,N}
    return quote
        B = BlockArray{T,N,Array{AbstractArray{T,N},N}}(undef_blocks, blocksizes...)
        axs = axes(B)
        @nloops $N block dim->blockaxes(B,dim) begin
            block_indices = @ntuple $N block
            indices = getindex.(axs, block_indices)
            block_size = length.(indices)
            B[block_indices...] = Fill{T}(T(), block_size...)
        end
        return B
    end
end
