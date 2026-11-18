# examples/cce/gcce_highspin_siv.m

A spin-3/2 vacancy-like centre with carbon-13 and silicon-29 spin-half bath nuclei undergoes free evolution at generalised CCE order two.

## Interface

`[answer,diagnostics]=gcce_highspin_siv()` returns numeric rows of column-major observed density matrices. `diagnostics.parameters.times` gives elapsed times; `diagnostics.reduction` contains order-resolved answers, validity masks, and numerical error bounds. No clipping or positivity repair is applied.

The example uses standard sys/inter/bas/parameters blocks, a product-ms transition selector, and `clexpn(spin_system,@free_evolution,parameters,'labframe')`. The selected manifolds define the equal coherent preparation; no state or operator matrix is assembled in the example. The plot shows observed populations.

The independent dense Hamiltonian, preparation, propagation, and reduced-density comparison live in `tests/kernel/test_cce_examples.m`.

The example takes no inputs and uses times in seconds and coupling tensors in Hz. It creates the spin system without a full-bath basis and uses a local one-worker pool; do not put Spinach and EasySpin on the same MATLAB path. This deliberately small illustrative model is not a publication-scale material reproduction. Quantitative comparisons require material structures, tensor frames, defect provenance, and source licences; the finite bath does not establish a thermodynamic T2 or microscopic T1.

[MATLAB source](../../../../../examples/cce/gcce_highspin_siv.m) · [Independent reference tests](../../../../../tests/kernel/test_cce_examples.m)
