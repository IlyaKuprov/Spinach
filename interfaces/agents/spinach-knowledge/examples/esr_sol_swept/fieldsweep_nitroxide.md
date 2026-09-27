# examples/esr_sol_swept/fieldsweep_nitroxide.m

- Signature: `fieldsweep_nitroxide()`

## Purpose

Compute a field-swept nitroxide EPR spectrum from resonance fields and transition moments. The source notes a calculation time of seconds.

## Physical / mathematical content

- The spin system contains an electron (`E`) and a nitrogen-14 nucleus (`14N`).
- The electron has an anisotropic Zeeman tensor with diagonal elements 2.01045, 2.00641, and 2.00211. The electron–nitrogen coupling tensor, scaled by `1e7`, has diagonal elements 1.2356, 1.1266, and 8.2230 and symmetric x–z elements of 0.6322.
- The calculation uses the `zeeman-hilb` formalism without a basis approximation and starts from `-state(spin_system,'Lz','E')` for the high-temperature approximation.

## Numerical / algorithmic content

- `fieldsweep` evaluates the spectrum using the `rep_2ang_100pts_sph` orientation grid, a 9 GHz microwave frequency, and `rspt_order=Inf`.
- The magnetic-field window is 0.316–0.326 T with 1024 points. The specified linewidth and tolerances are `fwhm=1e-5`, `int_tol=10.0`, and `tm_tol=0.1`.

## Implementation structure

- Set the isotopes and `sys.magnet=1`, then populate the Zeeman and coupling tensors.
- Create the spin system with `create`, apply the basis with `basis`, set the experiment parameters and initial state, and call `[spec,parameters]=fieldsweep(spin_system,parameters)`.
- Plot `spec` against `parameters.b_axis`, labeling magnetic field in tesla and intensity in arbitrary units.