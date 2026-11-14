# examples/cce/cce_polarized_bath.m

- Signature: `result=cce_polarized_bath()`

## Purpose

Nuclear polarization and exhaustive versus sampled bath estimators. The example reports its finite-cluster result without claiming a published material fit.

## Spin system and sequence

Standard Spinach sys/inter inputs are passed to create; the CCE driver compiles metadata without creating a full-bath basis. The example is self-contained, takes no inputs, uses SI times and explicit Hz coupling tensors, and opens a plot. It uses a local one-worker pool. Spinach and EasySpin must not share the MATLAB path.

## Observable and limits

The returned structure is the public cce result with any printed comparison metrics added by the example. Density examples return all central density elements; the two-centre example checks trace/positivity before computing entanglement. Other examples plot complex coherence or noise. Validity masks and conservative numerical bounds remain visible; no clipping or positivity repair is applied.

Except for the explicit pinned NV fixture and analytic identities, these are deliberately small illustrative systems, not publication-scale material reproductions. Material structures, tensor frames, defect origins and source licences must be resolved before quantitative experimental comparisons. A finite bath does not establish a thermodynamic T2 or microscopic T1. Read the header and printed source gate in [cce_polarized_bath.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/cce/cce_polarized_bath.m).
