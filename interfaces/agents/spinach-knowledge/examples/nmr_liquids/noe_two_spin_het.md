# examples/nmr_liquids/noe_two_spin_het.m

- Signature: `noe_two_spin_het()`

## Purpose

Nuclear Overhauser effect in a heteronuclear two-spin system in the short correlation time case. Calculation time: seconds.

## Physical / mathematical content

- The model contains one proton and one carbon-13 spin, separated by 1.03 Å, with zero isotropic Zeeman offsets. Redfield relaxation uses a 100 ps correlation time, 298 K, and the Di Bari equilibrium convention.
- The initial density operator is thermal equilibrium with the proton spin inverted. The calculation tracks the longitudinal magnetization of both spins, showing their NOE relaxation response.

## Numerical / algorithmic content

- The calculation uses the full spherical-tensor Liouville basis (no basis approximation) and retains the Redfield relaxation terms with the kite selection.
- The relaxation-only evolution is sampled every 0.01 s for 400 intervals, covering 0–4 s.

## Implementation structure

- Create the H-1/C-13 system at 14.1 T, construct the basis and Redfield superoperator, and calculate thermal equilibrium.
- Invert the proton's `Lz` component, then call multichannel evolution with both spins' `Lz` operators as detection channels; plot and label the proton and carbon longitudinal signals.
