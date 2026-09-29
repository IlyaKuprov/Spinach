# etc/diamond_defects/diamond_p1_13c.m

- MATLAB implementation: [etc/diamond_defects/diamond_p1_13c.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/diamond_defects/diamond_p1_13c.m)

**Call:** `[sys,inter] = diamond_p1_13c(parameters)`

Builds a P1-centre spin model containing one electron, one nitrogen nucleus and eighteen `13C` neighbours.

## Inputs

- `parameters.orientation`: `'111'`, `'110'` or `'100'`; the selected crystal direction is aligned with the magnetic field.
- `parameters.nitrogen`: `'14N'` or `'15N'`.

Both fields are required; the source checks these character values exactly. The electron/nitrogen parameter convention is documented as following `diamond_p1.m`.

## Model details

The carbon sites are labelled `G1_C1` through `G18`; the source identifies five measured/assigned sites (G1_C1, G2_C3, G4_C5, G8_C2 and G14_C4), while the other thirteen carry `_calc` labels and use Peaker et al.'s calculated Table 4 values. For example, G1_C1 has principal hyperfine values 139.531, 139.531 and 338.171 MHz in the source table; the code scales the table by `1e6` when constructing the tensor.

The tensors and site assignments draw on Barklie & Guven (1981), Cox et al. (1994), and Peaker et al. (2016). The coordinate model is deliberately representative, not a set of relaxed Cartesian coordinates: it uses ideal diamond-lattice directions, Peaker's tabulated distances `d/a0`, and `a0 = 3.567` Å. Nitrogen is at the origin. The electron coordinate is left empty so Spinach does not add point-dipolar electron–nuclear couplings on top of the explicit measured/calculated hyperfine tensors. Each carbon tensor's principal axes are rotated into the chosen field orientation before insertion into `inter.coupling.matrix`.

## Outputs and scope

- `sys`: isotope and label arrays for the electron, selected nitrogen isotope and 18 carbons, with the reconstructed nuclear coordinates.
- `inter`: Zeeman and coupling data, including the carbon hyperfine tensors.

The carbon count and site parameters are fixed by the routine; it does not expose a subset-size parameter. Use the returned tensors as the source model provides them rather than interpreting the representative coordinates as relaxed structure.

## References

- R. C. Barklie and J. Guven, *J. Phys. C: Solid State Phys.* **14**, 3621–3631 (1981), [doi:10.1088/0022-3719/14/25/009](https://doi.org/10.1088/0022-3719/14/25/009).
- A. Cox, M. E. Newton and J. M. Baker, *J. Phys.: Condens. Matter* **6**, 551–563 (1994), [doi:10.1088/0953-8984/6/2/012](https://doi.org/10.1088/0953-8984/6/2/012).
- C. V. Peaker et al., *Diamond and Related Materials* **70**, 118–123 (2016), [doi:10.1016/j.diamond.2016.10.013](https://doi.org/10.1016/j.diamond.2016.10.013).

[Spin Dynamics Wiki source page](https://spindynamics.org/wiki/index.php?title=diamond_p1_13c.m).
