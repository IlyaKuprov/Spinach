# examples/relaxation_theory/dd_relaxation_2.m

- MATLAB implementation: [examples/relaxation_theory/dd_relaxation_2.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/dd_relaxation_2.m)

- Signature: `dd_relaxation_2()`.
- Returns: no MATLAB output arguments; prints the relaxation superoperator in full (dense) form.

## Purpose

Build and display a complete Bloch–Redfield–Wangsness relaxation superoperator for three dipole-coupled spins. The dipolar couplings are derived from Cartesian coordinates. The source comment states that the result should be independent of the coordinate rotation; the script constructs one rotated geometry but does not itself compare multiple rotations or assert that invariance.

## Model and calculation

The example uses three `1H` spins, `sys.magnet=14.1`, and Zeeman scalar entries `{1.0,2.0,3.0}`. It rotates three fixed coordinate rows by `euler2dcm([pi/3 pi/4 pi/5])`; the coordinates are multiplied by the resulting matrix, and their units are not specified in this file. Relaxation is Redfield with `inter.equilibrium='zero'`, `inter.rlx_keep='labframe'`, one correlation-time entry `1e-9`, and integration tolerance `sys.tols.rlx_integration=1e-5`. It builds an untruncated `sphten-liouv` basis, creates the spin system, then evaluates `relaxation(spin_system)`; `full(...)` converts the result for display with `disp`.

The function takes no parameters, so changing the field, offsets, coordinates, relaxation settings, or basis requires editing the script. The header describes the run as taking seconds, but actual runtime depends on MATLAB and the host. The dense display is intended for this small three-spin example, not a scalable way to inspect large superoperators. No FID or spectrum is calculated.
