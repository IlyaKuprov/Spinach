# etc/mex/compile_mex.m

- MATLAB implementation: [etc/mex/compile_mex.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/mex/compile_mex.m)

## Call and effect

Run `compile_mex()` in MATLAB from the Spinach installation that contains this function; it takes no arguments and returns no MATLAB output. It locates its installation root from `mfilename('fullpath')`, then calls `mex` sequentially to rebuild five C++ MEX targets, placing each binary beside its source:

| Source | Destination directory | Kernel |
|---|---|---|
| `kernel/line_shapes/lorentzcon.cpp` | `kernel/line_shapes` | Lorentzian convolution |
| `kernel/line_shapes/gausscon.cpp` | `kernel/line_shapes` | Gaussian convolution |
| `kernel/eigenfields/cubic_roots.cpp` | `kernel/eigenfields` | Cubic polynomial roots |
| `kernel/indexing/spsortrows.cpp` | `kernel/indexing` | Sparse-double row sorting |
| `kernel/indexing/spunicols.cpp` | `kernel/indexing` | Sparse-double unique columns |

The build calls select the `-R2018a` MEX API, optimisation (`-O`), and `-DNDEBUG`. The first three also pass the current `COMPFLAGS` and `LINKFLAGS` through to MATLAB's compiler invocation. The target list and output paths are fixed in the function; this is an in-place rebuild, not a configurable target selector or installer.

A working MATLAB `mex` compiler configuration and the source tree are prerequisites. The function contains no separate post-build test or result report, so successful return alone is not evidence that these kernels pass runtime tests.

Source: [compile_mex.m in Spinach](https://spindynamics.org/wiki/index.php?title=compile_mex.m).
