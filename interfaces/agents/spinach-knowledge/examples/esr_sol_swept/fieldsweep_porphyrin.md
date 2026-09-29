# examples/esr_sol_swept/fieldsweep_porphyrin.m

[Source code](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_swept/fieldsweep_porphyrin.m) · Function: `fieldsweep_porphyrin()`

This example computes a field-swept EPR spectrum of a copper porphyrin complex by finding resonance fields and transition moments; the source estimates minutes for the calculation. The specified spins are four `14N` nuclei, an electron (`E`), and `63Cu`. The electron Zeeman principal values are `[2.0509, 2.0509, 2.1801]`. The electron–copper coupling principal values are `[-70.9257, -70.9257, -575.0219]*1e6`, with zero Euler angles; the four electron–nitrogen scalar couplings are each `46.0345*1e6`. The source supplies these numerical interaction values and scale factors without explicit unit annotations.

## Calculation and scan

The basis uses `zeeman-hilb` with no approximation. `S4` symmetry is assigned to spins `[1 2 3 4]`, the four nitrogens. The source sets the high-temperature initial state to `-state(spin_system,'Lz','E')` and calls `fieldsweep`, which returns a swept spectrum; this is not a pulse sequence.

Scan inputs are `mw_freq=9.39e9` (unit not annotated in this source), `window=[0.27 0.35]`, and `npoints=512`. The plotted magnetic-field axis is labelled tesla. The code also sets `fwhm=5e-4` (unit not stated), `int_tol=1e3`, `tm_tol=0.1`, `rspt_order=Inf`, and the powder orientation grid `rep_2ang_100pts_sph`. `sys.magnet=1` is set without an additional comment here.

## Output and scope

The example plots `spec` intensity in arbitrary units versus `parameters.b_axis` in tesla. It is a compact six-spin model with the interactions listed above and the imposed nitrogen permutation symmetry; it does not specify additional molecular or environmental degrees of freedom. The source's runtime note is minutes.
