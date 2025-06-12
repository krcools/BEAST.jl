module SauterSchwabQuadrature1D

# -------- exportet parts
# types
export SauterSchwabStrategy1D
export CommonEdge, CommonVertex

# functions
export sauterschwab_parameterized1D, _MRWrules, reorder

# -------- included files
include("doublesauterschwabint.jl")

end