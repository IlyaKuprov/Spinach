# examples/cce/gcce_two_centres.m

- Signature: `result=gcce_two_centres()`

## Purpose

Two electron centres: Bell-state decay and a physicality-gated negativity. The example prints a comparison against a separately implemented dense Kronecker/expm calculation.

## Spin system and sequence

Standard Spinach sys/inter inputs are passed to create; the CCE driver compiles metadata without creating a full-bath basis. The example is self-contained, takes no inputs, uses SI times and explicit Hz coupling tensors, and opens a plot. It uses a local one-worker pool. Spinach and EasySpin must not share the MATLAB path.

## Observable and limits

The returned structure is the public cce result with any printed comparison metrics added by the example. Density examples return all central density elements; the two-centre example checks trace/positivity before computing entanglement. Other examples plot complex coherence or noise. Validity masks and conservative numerical bounds remain visible; no clipping or positivity repair is applied.

Except for the explicit pinned NV fixture and analytic identities, these are deliberately small illustrative systems, not publication-scale material reproductions. Material structures, tensor frames, defect origins and source licences must be resolved before quantitative experimental comparisons. A finite bath does not establish a thermodynamic T2 or microscopic T1. Read the header and printed source gate in [gcce_two_centres.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/cce/gcce_two_centres.m).
