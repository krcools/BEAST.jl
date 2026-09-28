shapetype(::RTRefSpace{T}) where T =
    @NamedTuple{value::SVector{3,T},divergence::T}

shapetype(::LagrangeRefSpace{T,0,3}) where T =
    @NamedTuple{value::T,derivative::T}

shapetype(::LagrangeRefSpace{T,D,3}) where {T,D} =
    @NamedTuple{value::T,curl::SVector{3,T}}
