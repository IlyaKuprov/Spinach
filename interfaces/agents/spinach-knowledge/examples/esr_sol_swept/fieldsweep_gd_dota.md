# examples/esr_sol_swept/fieldsweep_gd_dota.m

[Source code](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_swept/fieldsweep_gd_dota.m) · Function: `fieldsweep_gd_dota()`

This example calculates a powder-averaged, W-band field-swept ESR spectrum for a Gd(III) DOTA complex. The source specifies an `E8` electron, scalar Zeeman parameter `1.9918`, and an axial self-coupling tensor with principal values `[0.57e9, 0.57e9, -2*0.57e9]/3` and zero Euler angles. No nuclear spins, hyperfine terms, or anisotropic Zeeman tensor are specified in this model.

## Calculation and scan

The code builds a `zeeman-hilb` basis with no basis approximation and calls `fieldsweep`; its source header describes the calculation as exact diagonalisation. This is a field sweep, not a pulse sequence. The powder average uses `rep_2ang_100pts_sph`. The initial state is `-state(spin_system,'Lz','E8')`, identified in the source as the high-temperature approximation.

The scan settings are `mw_freq=90e9` (the source labels the example W-band but does not annotate the parameter's unit), field window `[3.05 3.4]`, and `npoints=4096`; the plotted field axis is labelled in tesla. The file sets `fwhm=2e-4`, `int_tol=10.0`, `tm_tol=0.1`, and `rspt_order=Inf`; units for `fwhm` and the tolerance parameters are not stated in this source. `sys.magnet=1` is explicitly commented “must be 1”; it is not the scan window.

## Output and scope

`fieldsweep` returns `spec` and updated parameters. The example plots `spec` against `parameters.b_axis`, with magnetic field in tesla and intensity in arbitrary units. The source estimates seconds for the calculation. Treat the result as the spectrum of the stated one-electron Hamiltonian, not a complete molecular model of DOTA; adding omitted nuclei or interactions would change the model.
