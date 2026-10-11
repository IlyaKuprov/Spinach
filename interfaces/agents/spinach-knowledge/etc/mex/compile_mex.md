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

The C++ build calls select the `-R2018a` MEX API, optimisation (`-O`), and `-DNDEBUG`. The first three also pass the current `COMPFLAGS` and `LINKFLAGS` through to MATLAB's compiler invocation. MATLAB releases before R2026b retain the five C++ builds and skip the CUDA target, whose internal layout contract is R2026b. Windows, Linux, and macOS all reach the CUDA compiler on newer releases. The target list and output paths are fixed in the function; this is an in-place rebuild, not a configurable target selector or installer.

The CUDA build uses `-R2018a`, `-O`, and `-lcuda` to link the CUDA driver API. The CUDA target requires MATLAB R2026b or later, `mexcuda`, a supported host C++ compiler, and the CUDA driver library. `mexcuda` can use the NVIDIA compiler bundled with MATLAB; a standalone CUDA toolkit is not required ([MathWorks documentation](https://www.mathworks.com/help/parallel-computing/mexcuda.html)). The source tree is required for every build. Windows adds `CUDA_PATH/lib/x64` before `-lcuda` when `CUDA_PATH` is set, as required by [MATLAB library-link option ordering](https://www.mathworks.com/help/matlab/ref/mex.html); otherwise compiler library search remains unchanged. No operating system is excluded before calling the compiler. On macOS, the specific `parallel:gpu:device:MacOSGPUNotSupported` error is reported without aborting the completed C++ builds; every other compiler error propagates. Modern MATLAB itself has no CUDA GPU support on macOS ([MathWorks requirements](https://www.mathworks.com/support/requirements/matlab-mac.html)); retaining a compilation route does not provide that missing platform capability. The CUDA source probes host pointers through `ReadProcessMemory` on Windows and POSIX calls on Linux/macOS. The function contains no separate post-build test or result report, so successful return alone is not evidence that these kernels pass runtime tests.

Source: [compile_mex.m in Spinach](https://spindynamics.org/wiki/index.php?title=compile_mex.m).
