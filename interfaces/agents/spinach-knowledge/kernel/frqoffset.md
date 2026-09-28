# kernel/frqoffset.m

- Signature: `H=frqoffset(spin_system,H,parameters)`

## Purpose

Adds Larmor-frequency offsets to selected spins in a Hamiltonian or commutation superoperator; this is useful in liquid-state NMR.

## Physical / mathematical content

- For each selected spin, adds `2*pi*offset*Lz(spin)` to `H`; `offset` is given in Hz.
- The transformation is approximate. Use `rotframe.m` or `intrep.m` for a rigorous treatment of second-order effects.

## Parameters / inputs

- `parameters.spins` — cell array of spin labels to which offsets are applied (for example, `{'1H','13C'}`).
- `parameters.offset` — vector of offsets in Hz, one for each listed spin. If multiple channels refer to the same spin, their offsets must agree.

## Outputs

- `H` — the Hamiltonian operator or commutation superoperator with the offsets added.

<https://spindynamics.org/wiki/index.php?title=frqoffset.m>
