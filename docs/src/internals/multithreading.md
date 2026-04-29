# Multi-threading and assemble chain of custody

The essentials of the assembly process are explained in [assemble](assemble.md). To optimally use the computational capacity offered by multi-threaded systems and to account for usage scenarios as they appear in e.g. hierarchical matrix based compression algorithms, the work required for assembly is split over a number of functions, each responsible for part of the total work.

- `assemble(op, tfs, bfs; kwargs...)`:  Builds an AbstracMatrix representation of the discrete operator `(op, tfs, bfs)`. Is responsible for allocating storage for the result, building the quaddata cache and multi-threading.

- `assemble!(op, tfs, bfs, store, threading; quadstrat, scheduler)`: Builds an AbstractMatrix representation of the discrete operator `(op, tfs, bfs)` into preallocated storage. Is responsible for building the quaddata cache and multi-threading (dofspltting).

- `assemblechunk!(op, tfs, bfs, store; quadstrat)`: Builds an AbstractMatrix representation of the discrete operator `(op, tfs, bfs)`. Is responsible for building the quaddata cache and multi-threading (cellcoloring). Its responsibilities coincide with those of `assemble!` but this layer in the assembly stack is nevertheless needed to allow for *dof-splitting*, multi-threading based on the distribution of matrix rows and columns over different, potentially concurrent, tasks.

```julia
function assemblechunk_body!(op, tfs, tels,
    tad, telptrs, bfs,
    bels, bad, belptrs,
    qd, zlocal, store; quadstrat)
```

Runs on a single thread. Takes a pre-allocated zlocal for the storage of element-element interactions and a pre-built quaddata cache.

```julia
function assembleblock_body!(biop::IntegralOperator,
        tfs, test_ids, test_elements, test_assembly_data, active_test_el_ids,
        bfs, trial_ids, bsis_elements, trial_assembly_data, active_trial_el_ids,
        quadrature_data, zlocals, store; quadstrat)
```

*Goal*: consolidate assemblechunk_body and assembleblock_body
*Strategy*: new assembleblock_body:

- construct the Vector of ptrs to active cells
- build a reduced assemblydata structure
- call assemblychunk_body


New function signature:
```julia
    assemblechunk_body!(op, tfs, bfs,
        tels, tptrs, tad, active_tels,
        bels, bptrs, bad, active_bels,
        quaddata, store; quadstrat)
```

**Arguments**:
- `tels::Vector{<:Chart}`: a cached vector of charts that covers an area that contains the support that will be considered in the assembly process. In general this covered area is larger than the area considered in this call but smaller than the total mesh `geometry(tfs)` on which the test functions are defined. 

- `tad`: assembly data of size `num_active_tfs x (num_active_tels, num_tshapes)`. The row indices refer to the storage index for the corresponding test function. This row index does not necessarily coincide to the index of this test function within `tfs`. This generality enables using `assemblechunk_body!` to build subblocks of the total BEM interaction matrix (in e.g. H-matrix methods). The column index `p` refers to the index within `active_tels`. The index `P = active_tels[p]` gives the index within `tels` of the corresponding chart. This index `P` should also be passeed to `quadrule` to allow for efficient retrieval of any cached quadrature data relating to this element.

- `active_tels::Vector{Int}`: a vector of indices into `tels` such that the `tels[active_tels]` is the support that will be considered in the assembly process. This generality allows for `assemblechunk_body!` to be used with lock-free threading approaches based on coloring schemes. In simple single-threaded scenarios it holds that `length(active_tels)` equals the number of charts in the support of all the active test functions.

*Problems*:

- `tels`, `tel_ptrs`, and `qd` are to be considered caches, that are computed for an area of the mesh that in general is larger than the support of `tfs` but smaller than the complete mesh `geometry(tfs)`. This is where all the confusion and complexity stems from.
