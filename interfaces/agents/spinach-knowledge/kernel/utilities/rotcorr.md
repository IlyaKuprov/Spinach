# rotcorr: solvent-aware rigid-body rotational diffusion

`[tau,D,axes_len]=rotcorr(atom_symbols,xyz,solvent,temperature)` estimates rotational diffusion from one rigid molecular geometry in dilute pure solvent. Supply an N-by-1 cell array of element symbols and an N-by-3 real floating-point coordinate array in Angstrom, a character string `'water'` or `'chloroform'`, and temperature in Kelvin. Supported elements are H, B, C, N, O, F, Si, P, S, Cl, Se, Br, and I. Unsupported symbols and coincident atom centres are rejected. Use a complete atom list and the correct oligomeric state.

For example, a complete ideal tetrahedral methane geometry can be supplied as:

```matlab
xyz=[0 0 0;1 1 1;1 -1 -1;-1 1 -1;-1 -1 1]*1.09/sqrt(3);
[tau,D,axes_len]=rotcorr({'C';'H';'H';'H';'H'},xyz,'chloroform',298.15);
```

This illustrates the input convention, not validated molecular-scale hydrodynamic accuracy.

## Fixed solvent models

The elemental van der Waals radii come from Mantina et al., DOI [10.1021/jp8111556](https://doi.org/10.1021/jp8111556). Water adds a 2.2 Å shell and uses a 1.4 Å probe, following ROTDIF3's ELM settings. Chloroform uses an explicitly unhydrated van der Waals envelope with zero shell and zero probe; this is not a calibrated chloroform solvation model. Temperature-dependent viscosities represent ordinary H2O and CHCl3, not D2O/CDCl3, mixtures, or pressure-dependent fluids. Accepted temperatures are 273.16–373 K for water and 210–334 K for chloroform. The shipped `etc/rotcorr_sources.md` contains the complete radius table, viscosity equations, source links, and the distinctions between ELM, ROTDIF3, and HYDROPRO conventions.

## Calculation and interpretation

Exposed contact patches on atomic spheres are sampled with area weights. Their covariance defines an equivalent ellipsoid, and Perrin integrals give rotational diffusion. The design follows the surface-based ELM idea of Ryabov et al., JACS 128 (2006), 15432–15444, DOI [10.1021/ja062715t](https://doi.org/10.1021/ja062715t), but does not reproduce SURF tessellation or include re-entrant patches. The declared elemental radii also differ from the original protein model.

The routine doubles surface directions per atom from 256 until two consecutive refinements change `tau`, each sorted semi-axis, and diffusion along every direction by less than 1%. It raises an error if convergence is not reached at 1,048,576 directions per atom; it never silently returns the last unconverged estimate. This is numerical self-consistency rather than an absolute error bound. Rigid rotations change the finite-grid sampling error, not the underlying physical model.

`tau=1/(2*trace(D))` is an isotropic-equivalent rank-2 time in seconds. It is not the correlation time of every interaction under anisotropic rotation. `D` is the 3-by-3 rotational diffusion tensor in the input Cartesian frame, in inverse seconds; `axes_len` contains the three descending ellipsoid semi-axes in Angstrom.

## Applicability

The water shell is protein-derived, not a universal small-solute hydration calibration. Chloroform predictions are uncalibrated continuum estimates. Missing atoms, flexibility, specific binding, aggregation, and highly non-globular geometry can dominate quadrature error. The model retains anisotropy and exposed shape, but this extra capability is not a universal accuracy guarantee. For higher-fidelity shape-resolved hydrodynamics, consider validated bead/shell calculations such as HYDROPRO: Ortega et al., DOI [10.1016/j.bpj.2011.06.046](https://doi.org/10.1016/j.bpj.2011.06.046). Its effective-radius calibration is distinct from an ELM additive hydration layer.
