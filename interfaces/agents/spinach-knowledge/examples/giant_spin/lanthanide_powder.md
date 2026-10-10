# examples/giant_spin/lanthanide_powder.m

- MATLAB implementation: [examples/giant_spin/lanthanide_powder.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/giant_spin/lanthanide_powder.m)

- Source: [examples/giant_spin/lanthanide_powder.m](../../../../../examples/giant_spin/lanthanide_powder.m)
- Signature: `lanthanide_powder()` (no input arguments)

## Purpose

Calculates and plots a powder spectrum for Gd(III) using a giant-spin Hamiltonian with zero-field splitting (ZFS) through spherical rank 4. The source describes a sweepable 400 MHz NMR magnet and 263.2 GHz microwaves, and attributes the Stevens parameters to Gd(III) in tetragonal BaTiO3. It gives the reference DOI [10.1103/PhysRev.127.702](https://doi.org/10.1103/PhysRev.127.702) (Rimai and deMars).

## Spin model and basis

The source sets `sys.isotopes={'E8'}`, `inter.zeeman.scalar={1.9918}`, and `sys.magnet=1`. Its giant-spin coefficient arrays include a rank-2 component `[0, 0, -4.65e8, 0, 0]` and rank-4 components `[-2.00e5, 0, 0, 0, 3.34e6, 0, 0, 0, -2.00e5]`; all listed giant-spin Euler angles are zero. The source comments identify the rank-4 terms as converted from `b40=4e-4 cm^-1` and `b44=-2e-4 cm^-1`, and explain that odd ranks are zero by time-reversal symmetry. It builds a `zeeman-hilb` basis with `approximation='none'`.

## Powder field sweep

The high-temperature calculation calls `fieldsweep` with `spins={'E8'}`, grid `rep_2ang_100pts_sph`, microwave frequency `263.2e9` Hz, `fwhm=2e-4`, `int_tol=10.0`, `tm_tol=0.1`, field window `[9.32 9.56]` T, `npoints=4096`, and `rspt_order=Inf`. The window is in tesla as indicated by the plot's field-axis label; the source does not specify units for `fwhm`, `int_tol`, or `tm_tol`. The initial state is `-state(spin_system,'Lz','E8')`. The returned spectrum is plotted against `parameters.b_axis` with intensity labelled in arbitrary units. These are script settings, not a spectrum measured or reproduced for this note.
