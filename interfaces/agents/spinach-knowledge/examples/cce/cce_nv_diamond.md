# examples/cce/cce_nv_diamond.m

- Signature: `result=cce_nv_diamond()`

## Purpose

Conditional NV Hahn echo for the pinned PyCCE four-nucleus numerical fixture. The example prints the maximum complex discrepancy against the pinned PyCCE four-nucleus CCE2 fixture at commit d74d3cb683bf3caa07231654314e554f440833b4. The PyCCE off-diagonal trace is explicitly conjugated to Spinach’s convention. The source qubit is mS=-1/+1 at 500 G, not the planned 65 G ensemble application.

## Spin system and sequence

Standard Spinach sys/inter inputs are passed to create; the CCE driver compiles metadata without creating a full-bath basis. The example is self-contained, takes no inputs, uses SI times and explicit Hz coupling tensors, and opens a plot. It uses a local one-worker pool. Spinach and EasySpin must not share the MATLAB path.

## Observable and limits

The returned structure is the public cce result with any printed comparison metrics added by the example. Density examples return all central density elements; the two-centre example checks trace/positivity before computing entanglement. Other examples plot complex coherence or noise. Validity masks and conservative numerical bounds remain visible; no clipping or positivity repair is applied.

Except for the explicit pinned NV fixture and analytic identities, these are deliberately small illustrative systems, not publication-scale material reproductions. Material structures, tensor frames, defect origins and source licences must be resolved before quantitative experimental comparisons. A finite bath does not establish a thermodynamic T2 or microscopic T1. Read the header and printed source gate in [cce_nv_diamond.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/cce/cce_nv_diamond.m).
