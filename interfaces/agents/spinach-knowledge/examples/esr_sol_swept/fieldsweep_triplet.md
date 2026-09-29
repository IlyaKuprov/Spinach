# examples/esr_sol_swept/fieldsweep_triplet.m

[Source code](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_swept/fieldsweep_triplet.m) · Function: `fieldsweep_triplet()`

This example calculates a powder-averaged X-band field-swept ESR spectrum for a photo-generated pentacene triplet. Its spin system is `E3`, with an isotropic Zeeman tensor `diag([2.0 2.0 2.0])`. The source sets zero-field-splitting parameters `D=1360.1*1e6 Hz` and `E=-47.2*1e6 Hz` and converts them with `zfs2mat`. These correspond to `1360.1 MHz` and `-47.2 MHz`, respectively. The source estimates seconds for the calculation.

## Calculation and scan

The code uses a `zeeman-hilb` basis with no approximation and samples `rep_2ang_100pts_sph` orientations. It converts the Zeeman tensor to Hz/T and supplies an orientation- and field-dependent initial state through `zftrip`, using the rotated ZFS and Zeeman tensors and the vector `[0.56 0.31 0.13]`. It then calls `fieldsweep`; there is no pulse sequence in this field-swept calculation.

The explicitly unit-labelled settings are `mw_freq=9e9 Hz`, `fwhm=5e-4 Tesla`, and `window=[0.25 0.40] Tesla`. The scan has `npoints=512`, `int_tol=0.1`, `tm_tol=0.1`, and `rspt_order=Inf`. As in the source comment, `sys.magnet=1` is required.

## Output and scope

`fieldsweep` returns `spec`, plotted against `parameters.b_axis` with magnetic field in tesla and intensity in arbitrary units. The modeled system is the specified triplet electron and its orientation-dependent state; this example does not include nuclear spins or a broader pentacene host/defect environment.
