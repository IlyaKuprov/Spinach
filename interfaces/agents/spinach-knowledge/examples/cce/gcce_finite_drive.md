# examples/cce/gcce_finite_drive.m

Syntax: `[answer,diagnostics]=gcce_finite_drive()`.

Small electron–carbon model demonstrating a 20-microsecond rectangular observed
spin drive at 40000 rad/s, including the coupled nuclear drift, followed by free evolution. Canonical sys/inter,
transition preparation, a scalar parameters.drive descriptor, and generalised
clexpn with @finite_drive avoid manually assembled states and event arrays.
The returned answer is the column-vectorised observed density, one row per time;
`diagnostics.parameters.times` gives elapsed times; `diagnostics.reduction` contains order-resolved answers, validity masks, and numerical error bounds. No clipping or positivity repair is applied. The plot shows populations.
The independent two-spin expm reference lives in test_cce_examples.

The example takes no inputs and uses times in seconds and coupling tensors in Hz. It creates the spin system without a full-bath basis and uses a local one-worker pool; do not put Spinach and EasySpin on the same MATLAB path. This deliberately small illustrative model is not a publication-scale material reproduction. Quantitative comparisons require material structures, tensor frames, defect provenance, and source licences; the finite bath does not establish a thermodynamic T2 or microscopic T1.

[MATLAB source](../../../../../examples/cce/gcce_finite_drive.m) · [Independent reference tests](../../../../../tests/kernel/test_cce_examples.m)
