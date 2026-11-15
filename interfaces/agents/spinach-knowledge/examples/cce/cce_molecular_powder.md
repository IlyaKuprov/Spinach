# examples/cce/cce_molecular_powder.m

- Signature: `result=cce_molecular_powder()`

## Purpose

Orientation-averaged molecular Hahn echo with explicit quadrature weights. The example reports its finite-cluster result without claiming a published material fit.

## Spin system and sequence

Standard Spinach sys/inter inputs are passed to create; the CCE driver compiles metadata without creating a full-bath basis. The example is self-contained, takes no inputs, uses SI times and explicit Hz coupling tensors, and opens a plot. It uses a local one-worker pool. Spinach and EasySpin must not share the MATLAB path.

## Observable and limits

The returned structure is the public cce result with any printed comparison metrics added by the example. Density examples return all observed density elements; the two-centre example checks trace/positivity before computing spin correlations. Other examples plot complex signal or noise. Validity masks and conservative numerical bounds remain visible; no clipping or positivity repair is applied.

Except for the explicit pinned NV fixture and analytic identities, these are deliberately small illustrative systems, not publication-scale material reproductions. Material structures, tensor frames, defect origins and source licences must be resolved before quantitative experimental comparisons. A finite bath does not establish a thermodynamic T2 or microscopic T1. Read the header and printed source gate in [cce_molecular_powder.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/cce/cce_molecular_powder.m).

Uses generalized CCE with an explicit observed superposition: powder rotation creates transverse observed terms, so fixed projected manifolds are not assumed.

Response aliases, by-order values and masks, numerical bounds, and raw propagation records are now aggregated over the same positive quadrature weights; numerical diagnostics and budget flags are recomputed for the powder average. Orientation-resolved absolute bounds are retained in `orientation_abs_error_bound`.
