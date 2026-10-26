using Literate

Literate.markdown("./examples/composedoperator.jl","./docs/src/composedoperator"; execute=true)
Literate.markdown("./examples/pmchwt_theta.jl","./docs/src/projectors"; execute=true)
Literate.markdown("./examples/pmchwt_lf.jl","./docs/src/projectors"; execute=true)

Literate.markdown("./examples/efie_geop2.jl","./docs/src/tutorials"; execute=true)
