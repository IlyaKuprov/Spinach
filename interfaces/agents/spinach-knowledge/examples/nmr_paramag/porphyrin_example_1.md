# examples/nmr_paramag/porphyrin_example_1.m

- Signature: `porphyrin_example_1()`
- Source: [examples/nmr_paramag/porphyrin_example_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_paramag/porphyrin_example_1.m)
- Manual cited by the source: [Pseudocontact shift analysis](http://spindynamics.org/wiki/index.php?title=Pseudocontact_shift_analysis)

## Point-model calculation

The source defines 12 porphyrin-ring proton coordinates and places the metal at `mxyz=[0 0 0]`. It supplies diagonal g-tensors `diag([3.0 3.0 2.0])` for Co(II) and `diag([2.0 2.0 2.2])` for Cu(II), then calls `g2chi(g,298,1/2)` to form each Curie susceptibility tensor. The source does not annotate the unit of the `298` argument or the units of the coordinate values.

`ppcs(nxyz,mxyz,chi)` computes point-model PCS for each ion. The example displays the two result columns in Co, Cu order and labels the PCS output in ppm.

## Scope and omissions

This is a basic Cu(II)/Co(II) porphyrin comparison, not a carbonic-anhydrase calculation: there is no protein, residue, or distinct metal site. It does not use a distributed-density model, fit parameters, experimental spectra, a specified magnetic field, or a field/temperature sweep. The only temperature-like input is the unlabelled value 298 passed to `g2chi`.
