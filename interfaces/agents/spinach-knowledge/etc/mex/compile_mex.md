# etc/mex/compile_mex.m

- MATLAB implementation: [etc/mex/compile_mex.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/mex/compile_mex.m)

## Call and effect

Run `compile_mex()` in MATLAB from the Spinach installation that contains this function; it takes no arguments and returns no MATLAB output. It locates its installation root from `mfilename('fullpath')`, then calls `mex` and `mexcuda` sequentially to rebuild five C++ MEX targets and the CUDA sparse-product MEX, placing each binary beside its source:

| Source | Destination directory | Kernel |
|---|---|---|
| `kernel/line_shapes/lorentzcon.cpp` | `kernel/line_shapes` | Lorentzian convolution |
| `kernel/line_shapes/gausscon.cpp` | `kernel/line_shapes` | Gaussian convolution |
| `kernel/eigenfields/cubic_roots.cpp` | `kernel/eigenfields` | Cubic polynomial roots |
| `kernel/indexing/spsortrows.cpp` | `kernel/indexing` | Sparse-double row sorting |
| `kernel/indexing/spunicols.cpp` | `kernel/indexing` | Sparse-double unique columns |
| `kernel/overloads/cuda_sparse_by_sparse_mex.cu` | `kernel/overloads` | GPU sparse matrix product |

The C++ build calls select the `-R2018a` MEX API, optimisation (`-O`), and `-DNDEBUG`. The first three also pass the current `COMPFLAGS` and `LINKFLAGS` through to MATLAB's compiler invocation. The target list and output paths are fixed in the function; this is an in-place rebuild, not a configurable target selector or installer.

The CUDA build uses `-R2018a`, `-O`, and `-lcuda` to link the CUDA driver API. On Windows, it searches `CUDA_PATH/lib/x64` for the driver library. A working MATLAB `mex` compiler configuration, `mexcuda`, a supported CUDA toolkit, and the source tree are prerequisites; Windows requires `CUDA_PATH` to point to the toolkit. Platform support also depends on the CUDA source: its current POSIX host-memory probe is not Windows-compatible. The function contains no separate post-build test or result report, so successful return alone is not evidence that these kernels pass runtime tests.

Source: [compile_mex.m in Spinach](https://spindynamics.org/wiki/index.php?title=compile_mex.m).
