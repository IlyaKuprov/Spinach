# examples/esr_sol_swept/fieldsweep_nitroxide.m

[Source code](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_swept/fieldsweep_nitroxide.m) · Function: `fieldsweep_nitroxide()`

This example computes a field-swept EPR spectrum of a nitroxide by finding resonance fields and transition moments; the source estimates seconds for the calculation. The model contains an electron (`E`) and one `14N` nucleus. The electron Zeeman tensor is diagonal with entries `[2.01045, 2.00641, 2.00211]`. The electron–nitrogen coupling matrix is `[1.2356 0 0.6322; 0 1.1266 0; 0.6322 0 8.2230]*1e7`; the source gives no unit annotation for these entries.

## Calculation and scan

The example uses the `zeeman-hilb` formalism with `bas.approximation='none'`, samples orientations on `rep_2ang_100pts_sph`, and calls `fieldsweep`. It is a field-swept spectrum calculation rather than a pulse sequence. The initial state `-state(spin_system,'Lz','E')` is described in the source as the high-temperature approximation.

The supplied scan values are `mw_freq=9e9` (unit not annotated in this file), field window `[0.316 0.326]`, and `npoints=1024`. The plotted field axis is labelled tesla, so the window is a magnetic-field interval in the plotted axis. Other inputs are `fwhm=1e-5` (unit not stated), `int_tol=10.0`, `tm_tol=0.1`, and `rspt_order=Inf`. The source comments that `sys.magnet=1` is required.

## Output and scope

`fieldsweep` returns the spectrum `spec` and scan parameters; the example plots intensity in arbitrary units against `parameters.b_axis` in tesla. This is the two-spin model encoded in the example, not a general nitroxide conformational or solvent model. To represent a different nitroxide, update its interaction tensors and scan inputs rather than interpreting the plotted line shape as universal.
