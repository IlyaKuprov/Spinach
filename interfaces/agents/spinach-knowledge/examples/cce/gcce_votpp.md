# examples/cce/gcce_votpp.m

## Interface

`[answer,diagnostics]=gcce_votpp()` returns numeric rows of column-major observed density matrices. `diagnostics.parameters.times` gives elapsed times; `diagnostics.reduction` contains order-resolved answers and numerical bounds.

The example uses standard sys/inter/bas/parameters blocks, a product-ms transition selector, and `clexpn(spin_system,@free_evolution,parameters,'labframe')`. The selected manifolds define the equal coherent preparation; no state or operator matrix is assembled in the example. The plot shows observed populations. Original tensors, time samples, retained order, and numerical tolerances are unchanged.

The independent dense Hamiltonian, preparation, propagation, and reduced-density comparison live in `tests/kernel/test_cce_examples.m`. This is an illustrative finite model, not a material-specific quantitative reproduction.
