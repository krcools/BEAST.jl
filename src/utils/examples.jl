function example(name)
    fn = joinpath(dirname(@__FILE__),"..","..","examples",name*".jl")
    include(fn)
end
