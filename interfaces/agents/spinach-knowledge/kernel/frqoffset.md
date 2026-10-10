# kernel/frqoffset.m

- Signature: `H=frqoffset(spin_system,H,parameters)`

## Purpose

Adds frequency-offset terms for selected spins to a Hamiltonian operator or commutation superoperator.

## Frequency-offset rule

For each nonzero offset, the function adds `2*pi*offset*Lz(spin)` to `H`. An offset is specified in Hz; multiplying by `2*pi` converts its coefficient to angular frequency in radians per second. The function constructs an offset Hamiltonian term; it does not itself propagate a state or evolve a signal. Matrix addition leaves the returned `H` at the input dimensions and retains its operator or superoperator representation.

If multiple entries in `parameters.spins` name the same spin, their offset values must agree; the routine applies that spin offset once, rather than combining different channel values. Zero offsets contribute no term.

## Parameters / inputs

- `spin_system` — Spinach spin-system structure used to resolve the spin operators.
- `H` — Hamiltonian operator or commutation superoperator to which the offset terms are added.
- `parameters.spins` — non-empty cell array of character spin labels present in `spin_system.comp.isotopes`, for example `{'1H','13C'}`.
- `parameters.offset` — non-empty real numeric vector in Hz, with one value per spin label. The implementation checks the vector length and real-numeric form; it does not explicitly require finite values.

## Output

- `H` — the input operator or superoperator with the selected offset terms added.

This is the documented approximate offset transformation; the source recommends `rotframe.m` or `intrep.m` when a rigorous treatment of second-order effects is required.

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/frqoffset.m)
<https://spindynamics.org/wiki/index.php?title=frqoffset.m>
