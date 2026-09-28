# etc/diamond_defects/diamond_p1_13c.m

- Signature: [sys,inter]=diamond_p1_13c(parameters)

## Purpose

Builds a P1-centre spin system with one nitrogen nucleus and eighteen 13C neighbours. The electron and nitrogen parameters match diamond_p1.m. Carbon hyperfine tensors and site assignments are based on R. C. Barklie and J. Guven, *J. Phys. C: Solid State Phys.* **14**, 3621–3631 (1981), doi:10.1088/0022-3719/14/25/009; A. Cox, M. E. Newton, and J. M. Baker, *J. Phys.: Condens. Matter* **6**, 551–563 (1994), doi:10.1088/0953-8984/6/2/012; and C. V. Peaker et al., *Diamond Relat. Mater.* **70**, 118–123 (2016), doi:10.1016/j.diamond.2016.10.013.

Five sites have literature-based labels and assignments: G1_C1, G2_C3, G4_C5, G8_C2, and G14_C4. G1_C1 and G2_C3 use the Cox et al. Table 2 assignments; the sign for G2_C3 and the assignments of G4_C5, G8_C2, and G14_C4 follow Peaker et al. The remaining thirteen sites carry _calc labels and use Peaker et al.'s Table 4 values.

## Physical / mathematical content

Each 13C nucleus is coupled to the electron by its own anisotropic hyperfine tensor. For 14N, the nitrogen hyperfine tensor and quadrupole interaction match diamond_p1.m; for 15N, the opposite-sign hyperfine tensor is used and no quadrupole term is added. The electron coordinate is deliberately left unspecified so that point-dipolar electron–nuclear couplings are not added on top of the measured or calculated hyperfine tensors.

## Numerical / algorithmic content

The source reconstructs representative carbon coordinates from ideal diamond-lattice directions and Peaker et al.'s tabulated distances from the midpoint of the broken N–C bond, expressed in units of a0. Relaxed Cartesian coordinates are not reported in that paper. The code sets a0 = 3.567 Å, places nitrogen at the origin, and leaves the electron coordinate empty. It constructs each carbon tensor from its principal values and polar/azimuthal axes, symmetrises it, and rotates it consistently with the requested orientation. These coordinates are representative reconstructions, not reported relaxed positions.

## Parameters / inputs

- parameters.orientation: '111', '110', or '100'; the corresponding crystal-plane normal is aligned with the magnetic field.
- parameters.nitrogen: '14N' or '15N'. Both fields are required.

## Outputs

- sys: Spinach system specification structure, with the electron, nitrogen, and eighteen 13C nuclei and their labels.
- inter: Spinach interaction specification structure, including coordinates for nitrogen and carbon; the electron coordinate is empty.

## Implementation structure

The routine validates its two required fields, constructs the orientation rotation and site list, sets isotope and label arrays, builds the nitrogen interactions, reconstructs nuclear coordinates, and adds the eighteen rotated carbon hyperfine tensors to Spinach's coupling matrix.
