# BEAST CUDA extension

The extension provides CUDA specific implementations not yet GPU agnostic, we should also consider [KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl) to make AMD available in future. 

```julia
using BEAST, CUDA
ext = Base.get_extension(BEAST, :BEASTCUDAExt)
```

## Available functionality

- `assemble(...; threading=:gpu)` assembles a dense matrix on one or more GPUs.
  `devices`, `nstreams`, and `tilingstrat` control device selection,
  concurrency, and peak memory usage.
- `ext.gpu_blockassembler(...)` creates a reusable assembler for individual
  dense matrix blocks.
- `ext.gpu_batched_blockassemble(...)` assembles several blocks together. Its
  `budget` keyword limits the intermediate staging memory.
- `ext.gpu_blockassemblers(...; devices)` creates the per-device assemblers for
  multi-GPU batched block assembly.

The GPU results are tested against CPU assembly for Helmholtz and Maxwell
single-layer operators with Lagrange, Raviart-Thomas, and
Buffa-Christiansen spaces. Supported quadrature strategies are
`DoubleNumSauterQstrat` and `DoubleNumQStrat`.

See [`examples/gpu_assemble.jl`](../../examples/gpu_assemble.jl) and
[`examples/gpu_blockassembler.jl`](../../examples/gpu_blockassembler.jl).

## Current limitations

- Only three-dimensional triangular surface meshes are supported.
- Singular block assembly does not yet support test and trial meshes in a
  refinement relation. Separate nonintersecting meshes are supported.
- The default Wilton-Sauter quadrature strategy is not implemented on the GPU,
  so a supported strategy must be selected explicitly.
- CUDA is the only GPU backend. 
- Additional operators and reference spaces need tests before they can be considered supported.

## Code layout

- `gpu_assemble_integralop.jl`: dense assembly, quadrature kernels, tiling,
  streams, and multi-GPU scheduling.
- `gpu_blockassembler.jl`: reusable single-block assembly and GPU-resident
  space data.
- `gpu_batched_blockassembler.jl`: batched projection, chunking, and
  multi-GPU block distribution.
- `gpu_basis.jl` and `gpu_integrals.jl`: GPU shape and integrand support.
- `tiling.jl` and `gpu_utils.jl`: tiling types and kernel launch helpers.

GPU regression tests are located in `test/gpu`.
