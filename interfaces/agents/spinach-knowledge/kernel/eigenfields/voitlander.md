# kernel/eigenfields/voitlander.m

- Signature: `spec=voitlander(spin_system,parameters,triangle,Ic,Iz,Qc,Qz,Hmw)`

## Purpose

Integrates field-swept EPR transition contributions over a spherical triangle. The triangle vertices carry transition data evaluated at their orientations; the routine recursively subdivides the triangle until its integral estimate meets `parameters.int_tol`.

## Physical / mathematical content

At each triangle vertex, the transition contribution uses the transition moment, scaled field-sweep Jacobian, and energy-level population difference. Corresponding transitions are matched across vertices by transition identities when available, with field continuation used to match roots. The finite vertex-wise products `tm*tj*pd` are averaged, the three transition fields are supplied to the Lorentzian-triangle convolution, and the contribution is weighted by the spherical triangle area.

## Numerical / algorithmic content

The routine bisects the spherical edges to obtain three midpoint orientations, evaluates their eigensets with `eigenfields`, and forms four child triangles. It compares the direct triangle integral `spec_dir` with the sum over the four children `spec_sub`, then applies the Simpson-Richardson estimate `spec_sim = (4*spec_sub - spec_dir)/3`. If `norm(spec_dir - spec_sub, 2) > parameters.int_tol`, it recursively integrates the four children and sums their results; the child calls are submitted with `parfeval`. Otherwise it returns `spec_sim`.

## Parameters / inputs

- `triangle(1:3).xyz` - Cartesian coordinates of the spherical-triangle corners, unit column vectors.
- `triangle(1:3).tf` - transition fields at the corners, real column vectors, one element per transition.
- `triangle(1:3).tm` - transition moments at the corners, positive column vectors, one element per transition.
- `triangle(1:3).tw` - transition widths at the corners, positive column vectors, one element per transition.
- `triangle(1:3).pd` - energy-level population differences at the corners, real column vectors, one element per transition.
- `triangle(1:3).ti` - transition identity arrays at the triangle corners, one row per transition.
- `triangle(1:3).tj` - scaled field-sweep Jacobians at the triangle corners, real column vectors, one element per transition.
- `Ic` - isotropic part of the coupling Hamiltonian, a Hermitian matrix (set retention to `'couplings'` in `assume.m` and then call `hamiltonian.m`).
- `Qc` - irreducible components of the anisotropic part of the coupling Hamiltonian, a cell array returned by `hamiltonian.m`.
- `Iz` - isotropic part of the Zeeman Hamiltonian, a Hermitian matrix (set retention to `'zeeman'` in `assume.m` and then call `hamiltonian.m`), normalised to 1 Tesla.
- `Qz` - irreducible components of the anisotropic part of the Zeeman Hamiltonian, a cell array returned by `hamiltonian.m`, normalised to 1 Tesla.
- `Hmw` - perturbation operator, a Hermitian matrix.
- `parameters.b_axis` - vector of magnetic-field values, in Tesla.
- `parameters.int_tol` - integration accuracy tolerance; the implementation requires a positive real scalar.
- The implementation also validates `parameters.mw_freq`, `parameters.window`, `parameters.pp_tol`, `parameters.tm_tol`, and `parameters.fwhm`; for the `zeeman-hilb` formalism it also requires `parameters.rspt_order`.

## Outputs

- `spec` - ESR spectrum integral over the triangle, an array with the same dimensions as `parameters.b_axis`.

## Reference

[voitlander.m](https://spindynamics.org/wiki/index.php?title=voitlander.m)