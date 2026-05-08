using Literate

Literate.markdown("./examples/composedoperator.jl","./docs/src/composedoperator"; execute=true)
Literate.markdown("./examples/pmchwt_theta.jl","./docs/src/projectors"; execute=true)
Literate.markdown("./examples/pmchwt_lf.jl","./docs/src/projectors"; execute=true)
Literate.markdown("./examples/OSRC_efie.jl","./docs/src/operators"; execute=true)

