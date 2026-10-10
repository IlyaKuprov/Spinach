# rotcorr solvent and atomic parameters

`rotcorr(atom_symbols,xyz,solvent,temperature)` is a rigid-body, dilute-solution, stick-boundary continuum estimate, not a solute-specific solvation calculation. Coordinates are in Angstrom, temperature in Kelvin, and viscosity in Pa s. The fixed parameters below are embedded in `kernel/utilities/rotcorr.m`; no download is required.

## Elemental van der Waals radii

Mantina et al., *Consistent van der Waals Radii for the Whole Main Group*, J. Phys. Chem. A 113 (2009), 5806–5812, DOI [10.1021/jp8111556](https://doi.org/10.1021/jp8111556), [Table 12](https://pmc.ncbi.nlm.nih.gov/articles/PMC3658832/#T12):

| Element | H | B | C | N | O | F | Si | P | S | Cl | Se | Br | I |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Radius / Å | 1.10 | 1.92 | 1.70 | 1.55 | 1.52 | 1.47 | 2.10 | 1.80 | 1.80 | 1.75 | 1.90 | 1.83 | 1.98 |

These are elemental geometric parameters, not solvent-dependent hydrodynamic radii. Hydrogen is 1.10 Å in this table, not Bondi's older 1.20 Å convention. Unsupported elements and coincident atom centres are rejected. Supply the complete molecular atom list, including hydrogens when available; this is not a united-atom force-field model.

## Water surface

Water uses an additive 2.2 Å shell and a 1.4 Å rolling probe, the explicit defaults in ROTDIF3/ARMOR's [`ElmPredictor.java`](https://api.bitbucket.org/2.0/repositories/kberlin/armor/src/b540d7389b04a1934a84792198b0ea41d5019b50/src/main/java/edu/umd/umiacs/armor/nmr/relax/ElmPredictor.java). Its [`SolventAccessibleSurfaceImpl.java`](https://api.bitbucket.org/2.0/repositories/kberlin/armor/src/b540d7389b04a1934a84792198b0ea41d5019b50/src/main/java/edu/umd/umiacs/armor/molecule/SolventAccessibleSurfaceImpl.java) distinguishes the atomic surface radius plus shell from the probe-centre radius. In `rotcorr`, the probe tests accessibility but is not added a second time to the contact-surface radius.

These settings originate in protein hydrodynamics. Combining them with the declared elemental table and sampled contact patches is an ELM-inspired approximation, **not a reproduction or recalibration of ROTDIF3 or the original ELM model**. Ryabov et al., JACS 128 (2006), 15432–15444, DOI [10.1021/ja062715t](https://doi.org/10.1021/ja062715t), use SURF triangulation and discuss a 2.8 Å shell with different atomic radii (e.g. C 1.85, N 1.75, O 1.60 Å). Surface weighting and re-entrant patches differ. Neither protein shell value establishes hydration of an arbitrary small molecule.

HYDROPRO is a distinct bead/shell method: Ortega, Amorós, and García de la Torre, Biophys. J. 101 (2011), 892–898, DOI [10.1016/j.bpj.2011.06.046](https://doi.org/10.1016/j.bpj.2011.06.046), [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC3175065/). Its 2.9 Å atomic parameter is a **total effective radius per nonhydrogen atom**, not an additive hydration layer. It is therefore not substituted for the ELM layer here.

## Chloroform surface

Chloroform uses zero additive shell and zero probe: exposed patches of the union of elemental van der Waals spheres. This is an explicitly **unhydrated envelope approximation**, not a sourced chloroform hydration thickness or a calibrated solvent-accessible surface. No water probe or protein hydration layer is transferred to chloroform. The hydrodynamic accuracy of this approximation for small solutes has not been established; specific solvation and molecular-scale slip can dominate the result. Numerical convergence does not remove this limitation.

## Temperature-dependent dynamic viscosity

Only ordinary pure H2O and CHCl3 at approximately ambient pressure are represented. D2O, CDCl3, mixtures, dissolved salts, and pressure-dependent viscosities are not represented.

For water, let `x=T/300`:

```
eta = 1e-6*(280.68*x^(-1.9) + 511.45*x^(-7.7)
           + 61.131*x^(-19.6) + 0.45903*x^(-40))
```

Source: Assael et al., *Reference Values and Reference Correlations for the Thermal Conductivity and Viscosity of Fluids*, J. Phys. Chem. Ref. Data 47 (2018), DOI [10.1063/1.5036625](https://doi.org/10.1063/1.5036625), [equation (8)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6463310/#FD8). The published 0.1 MPa range is 253.15–383.15 K, including metastable regions; expanded relative uncertainty is 1.5% at 95% confidence. `rotcorr` deliberately restricts water to **273.16–373 K**. This compact correlation is not the full density-dependent IAPWS-2008 equation. It gives 1.001567 mPa s at 293.15 K and 0.889997 mPa s at 298.15 K.

For chloroform, with natural logarithm:

```
eta = exp(-14.109 + 1049.2/T + 0.5377*log(T))
```

Source: the `chemicals` package's [Perry Table 2-313 viscosity data](https://raw.githubusercontent.com/CalebBell/chemicals/master/chemicals/Viscosity/Table%202-313%20Viscosity%20of%20Inorganic%20and%20Organic%20Liquids.tsv), CAS 67-66-3, coefficients `[-14.109,1049.2,0.5377,0,0]`, interpreted with [DIPPR equation 101](https://raw.githubusercontent.com/CalebBell/chemicals/master/chemicals/dippr.py). The package attributes the data to Green and Perry, *Perry's Chemical Engineers' Handbook*, 8th ed. (2007). The [thermo implementation](https://raw.githubusercontent.com/CalebBell/thermo/master/thermo/viscosity.py) identifies the units as Pa s. This is package provenance to Perry, not an independently refitted experimental dataset. Its fit range is 209.63–353.2 K and no uncertainty is supplied. `rotcorr` restricts it to **210–334 K**, below the normal boiling point. It gives 0.566836 mPa s at 293.15 K and 0.538691 mPa s at 298.15 K.

## Numerical and physical applicability

The ellipsoid comes from area-weighted sampled contact-surface covariance, followed by Perrin rotational friction integrals. Re-entrant surface patches are absent. Surface points per atom double from 256; two successive refinements must change the scalar time, every sorted semi-axis, and diffusion along every direction by less than 1%. The tensor criterion is the largest absolute generalised eigenvalue of `(D_new-D_old,D_new)`, so it also checks slow rotational modes without dividing by zero off-diagonal entries. At 1,048,576 points per atom the routine raises an error rather than returning an unconverged answer. This finite cap limits quadrature memory; it does not guarantee quick execution for large molecules.

The 1% criterion measures successive-grid self-consistency, not a rigorous bound on discretisation error and certainly not 1% experimental accuracy. Use one rigid body of known oligomeric state. Flexible or highly non-globular structures need a more complete hydrodynamic model. `tau=1/(2*trace(D))` is an isotropic-equivalent rank-2 time, not the correlation time of every interaction in an anisotropically rotating molecule. The full tensor is returned in the input Cartesian frame.
