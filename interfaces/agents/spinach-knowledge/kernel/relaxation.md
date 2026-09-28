# kernel/relaxation.m

- Signature: `R=relaxation(spin_system,euler_angles)`

## Purpose

Builds the relaxation superoperator from the theories selected in `spin_system.rlx.theories`. Orientation-dependent contributions use `euler_angles` only when the selected theory supports relaxation-rate anisotropy; dissipative bosonic modes declared in `spin_system.inter.modes` contribute thermalised GKSL terms in Liouville-space formalisms.

## Physical / mathematical content

- Redfield is one available relaxation theory. Its rotationally averaged rates are not changed by `euler_angles`; orientation matters for theories with relaxation-rate anisotropy.
- Diagonal retention is not implemented for `zeeman-liouv` and is rejected explicitly. Spherical-tensor diagonal retention and full Zeeman (`labframe`) retention are available.
- SRSK is supported in `sphten-liouv`; a request in `zeeman-liouv` is rejected. Its recursive call sets equilibrium to `'zero'`, so IME or DiBari-Levitt thermalisation is applied once to the accumulated spin relaxation, not to the recursive SRSK contribution. Mode damping and dephasing are disabled in that recursive copy and left to the outer call.

## Numerical / algorithmic content

The selected spin-relaxation contributions are accumulated in `R`. Spinach context functions add relaxation and kinetics superoperators to the total Liouvillian automatically. A manually assembled Liouvillian must include this dissipative term with the factor `1i`: use `L=H+1i*R+1i*K`, not `H+R+K`.

## Parameters / inputs

- `euler_angles` - three Euler angles in the active ZYZ convention, in radians, specifying system orientation relative to the input orientation. They are required by theories that support relaxation-rate anisotropy and have no effect on theories such as Redfield that use rotational averaging.

## Outputs

- `R` - relaxation superoperator. In Liouville-space formalisms, dissipative bosonic modes declared in `inter.modes` receive thermalised GKSL dissipators from `rlx_modes.m`; `euler_angles` refers to the spin subsystem and does not affect mode terms.

## Source

- <https://spindynamics.org/wiki/index.php?title=relaxation.m>
