# Integrity, include, and compiler regression invariants

`test_dynamic_integrity_includes_mex_suite()` returns a standard regression
result. Its tests exercise infrastructure without rebuilding shipped binaries,
rewriting the production integrity database, or running destructive cleanup.

The checked host defaults are 32 process workers on ALAUNDO and 12 on TALOS
when GPU arithmetic is requested. Host-specific parallelisation must preserve
user overrides; unknown hosts and
CPU-only routes must not acquire unsolicited GPU worker settings. GPU guard
includes must remove and restore only the requested GPU flag, preserve other
flags, and reject an unmatched end guard. Environment, graphics defaults, and
any temporary pool created by the test are restored on exit.

Parallel profiling must return finite communication and non-negative timing
measurements. Serial and asynchronous Redfield includes must agree with the
analytical one-dimensional exponential integral and release their large
intermediates. Direct include dispatch must preserve these same behaviours.
For the scalar rank-one Redfield fixture, the analytical reference is
`R=-(1-exp(-T))/3`, with `T=1.5*log(1/rlx_integration)` and
`rlx_integration=1e-6`; absolute and relative comparison tolerances are `1e-10`.
The detailed internal profiler path is inspected statically, not executed.

Startup checks must use canonical production utilities. Integrity tools must
reject invalid modes/actions and empty selections; a freshly rearmed isolated
tree must pass its hash check. Production sniffing is read-only. Destructive
cleanup and continuous patrol are deliberately inspected rather than run.

The compiler helper must confine outputs to their source directories, retain
the interleaved-complex API and CUDA driver link, and avoid OpenMP dependence.
Windows, Linux, and macOS on R2026b or later include the CUDA gateway; older
MATLAB releases retain only the five C++ targets. Windows driver-library
lookup is checked both with an external CUDA path and with compiler defaults.
The specific macOS CUDA-unavailable error must be reported without aborting
the C++ builds; other macOS errors and all Windows errors must propagate.
Compiler invocations are intercepted in an isolated fixture, so these checks never rebuild production
binaries. This validates routing and compiler arguments, not kernel arithmetic;
`kernel/cuda_product` supplies the independent CPU-reference numerical checks.

Run `run_test('kernel/dynamic_integrity_includes_mex')` from the main suite.
