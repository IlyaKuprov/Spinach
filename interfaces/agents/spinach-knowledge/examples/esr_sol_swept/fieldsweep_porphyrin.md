# examples/esr_sol_swept/fieldsweep_porphyrin.m

- Signature: `fieldsweep_porphyrin()`

## Purpose

Simulate a field-swept EPR spectrum of a copper porphyrin complex by finding resonance fields and transition moments. The source notes a calculation time of minutes.

## Physical / mathematical content

- The spin system contains one electron (`E`), one `63Cu` nucleus, and four equivalent `14N` nuclei.
- The electron has an axial Zeeman tensor with principal values `[2.0509 2.0509 2.1801]`. Its hyperfine tensor with `63Cu` has principal values `[-70.9257 -70.9257 -575.0219]` MHz; each `14N` has a scalar electron hyperfine coupling of `46.0345` MHz.
- The calculation uses a 9.39 GHz microwave frequency and samples orientations on the `rep_2ang_100pts_sph` grid. The initial state is `-state(spin_system,'Lz','E')` for the high-temperature approximation.

## Numerical / algorithmic content

- The calculation uses the `zeeman-hilb` formalism with no basis approximation and applies `S4` symmetry to the four nitrogen spins.
- `fieldsweep` calculates 512 spectrum points over a 0.27–0.35 T field window. Parameters include a `5e-4` T linewidth, intensity tolerance `1e3`, transition-moment tolerance `0.1`, and `rspt_order=Inf`.

## Implementation structure

The function preallocates interaction cells, specifies the Zeeman and hyperfine interactions, and constructs the Spinach spin system with `create` and `basis`. It sets the experiment parameters and initial state, calls `fieldsweep`, then plots intensity against `parameters.b_axis` with the magnetic field labeled in tesla.