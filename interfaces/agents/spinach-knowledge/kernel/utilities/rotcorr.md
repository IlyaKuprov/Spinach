# rotcorr: rigid-body surface-ellipsoid rotational diffusion

`rotcorr(xyz,radii,layer,probe,temp,visc,npoints)` estimates an isotropic-equivalent rank-2 correlation time and the full rotational diffusion tensor from a rigid atom geometry. Coordinates and radii are in Angstrom, temperature is in Kelvin, and dynamic viscosity is in Pa s. Every solvent and surface-model parameter is explicit: atom coordinates cannot determine hydration or viscosity.

The estimator samples exposed contact patches on hydration-inflated atom spheres, forms an area-weighted surface covariance ellipsoid, and applies Perrin's rotational friction integrals. Its design follows the surface-based ELM idea of Ryabov, Geraghty, Varshney, and Fushman, JACS 128 (2006), 15432–15444, DOI [10.1021/ja062715t](https://doi.org/10.1021/ja062715t). It is not an exact implementation of that paper's SURF tessellation: re-entrant patches are omitted and the contact patches use numerical spherical sampling.

The returned scalar is `1/(2*trace(D))`. For anisotropic rotation it is not the correlation time of every interaction, and it must not replace an orientation-dependent diffusion calculation without an isotropic approximation. The tensor is returned in the input coordinate frame. Increase `npoints` to establish convergence; rigid rotation changes finite-grid sampling error but not the converged physical result.

## Choosing the model

The method improves on a uniform-atom-radius, volume-equivalent sphere by retaining atom-specific radii, solvent accessibility, hydration, and rotational anisotropy. This is greater model capability, not a universal accuracy guarantee. For strongly irregular shapes or higher-fidelity hydrodynamics, use a validated bead/shell solver such as HYDROPRO; see Ortega, Amorós, and García de la Torre, Biophysical Journal 101 (2011), 892–898, DOI [10.1016/j.bpj.2011.06.046](https://doi.org/10.1016/j.bpj.2011.06.046). Protein hydration calibration should not be transferred uncritically to small molecules or nonaqueous solvents.

ROTDIF3's ELM source uses a 2.2 Angstrom shell and 1.4 Angstrom water probe, whereas the original paper discusses a 2.8 Angstrom hydration layer with SURF. These are distinct model settings, not interchangeable fitted constants. Supply radii and hydration appropriate to the desired model, and document them. Missing atoms, flexible termini, oligomeric state, and non-globular geometry can dominate numerical quadrature error.
