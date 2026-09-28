# kernel/utilities/rlx_modes.m

- Signature: `R=rlx_modes(spin_system)`

## Purpose

Builds the Liouville-space relaxation superoperator for dissipative bosonic modes described in `spin_system.inter.modes`. It includes thermal amplitude damping and number-operator pure dephasing for modes of types `C`, `V`, and `T`.

## Physical / mathematical content

For each mode, the implementation uses cooling and heating terms with coefficients `kappa*(1+nbar)` and `kappa*nbar`, respectively, plus pure dephasing with coefficient `2*gamma_phi`. Here `nbar` is the Bose-Einstein occupation at the mode's physical frequency. At zero temperature it is set to zero; the spin relaxation convention equating zero temperature with the high-temperature limit does not apply to these bosonic modes.

## Numerical / algorithmic content

For a damped mode at positive temperature, the physical frequency is the declared carrier plus the mode frequency when the carrier is positive, and otherwise the absolute value of the declared mode frequency. A frequency below `2*pi*spin_system.tols.inter_cutoff` causes an error because its thermal occupation is undefined. The superoperator is accumulated from left- and right-action ladder-operator superoperators; modes with both rates zero are skipped.

## Parameters / inputs

- `spin_system` - Spinach system structure containing bosonic mode information in `inter.modes` (including the damping and dephasing rates), relaxation temperature in `rlx.temperature`, and tolerances. The required formalism is `zeeman-liouv` or `sphten-liouv`; damping and dephasing rates are in s^-1.

## Outputs

- `R` - bosonic-mode dissipation superoperator. The implemented terms are `kappa*(1+nbar)*D[a]`, `kappa*nbar*D[c]`, and `2*gamma_phi*D[n]`, with ladder operators truncated to each mode's declared level count.

## Implementation structure

The function validates the mode information, Liouville-space formalism, and temperature field, preallocates `R`, and visits the bosonic modes. For each active mode it reports the damping and dephasing rates and the occupation number, constructs the needed left/right ladder-operator superoperators, and adds the cooling, heating (when `nbar>0`), and dephasing (when `gamma_phi>0`) contributions.

## Reference

[Spin Dynamics Wiki: rlx_modes.m](https://spindynamics.org/wiki/index.php?title=rlx_modes.m)
