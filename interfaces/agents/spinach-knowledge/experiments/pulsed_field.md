# experiments/pulsed_field.m

- Signature: `answer=pulsed_field(spin_system,parameters,H,R,K) %#ok<INUSD>`

## Purpose

Simulates magnetisation dynamics in a time-dependent laboratory-frame Z field with spin-phonon relaxation, as used in pulsed-field magnetometry of molecular magnets. The field is treated as a staircase, and the spin-phonon dissipator is rebuilt in the eigenbasis of the Hamiltonian on each stair.

## Numerical / algorithmic content

The initial state is thermal equilibrium of the field-free Hamiltonian at the phonon-bath temperature. Each stair uses the field profile at its midpoint. In the stair-Hamiltonian eigenbasis, propagation is a symmetric split: exact coherent half-stair phases, a second-order dissipative step, then the coherent phases again. The dissipator-times-stair-width must be small; the coherent part is treated exactly. The dissipator is applied using Hilbert-space matrix products, with cubic cost in Hilbert-space dimension. Recorded field values are evaluated at the ends of the recorded stairs.

## Parameters / inputs

- `parameters.field_prof` — function handle returning the field in Tesla at a time in seconds.
- `parameters.hzeeman` — Zeeman operator per Tesla in rad/s/T, supplied by the context when `zeeman_op` is requested in `parameters.needs`.
- `parameters.timestep` — stair width in seconds.
- `parameters.nsteps` — number of stairs.
- `parameters.coil` — Hermitian Hilbert-space observable operator or a cell array of such operators.
- `parameters.phonon_x` — spin-phonon coupling operator (see `rlx_phonon.m`).
- `parameters.phonon_i0` — phonon spectral-density prefactor (see `rlx_phonon.m`).
- `parameters.phonon_alpha` — phonon spectral-density exponent, at least 1 (see `rlx_phonon.m`).
- `parameters.nout` — number of stairs between recorded observable values.
- `H` — Hilbert-space Hamiltonian from the context, containing the Zeeman term at `sys.magnet=1` Tesla; the routine removes that term and adds the field for each stair.
- `R` — context relaxation superoperator; ignored because the spin-phonon dissipator is constructed at each stair.
- `K` — context kinetics superoperator; ignored.

## Outputs

- `answer.t` — column of recording times in seconds.
- `answer.field` — column of field values at those times in Tesla.
- `answer.obs` — matrix of observable expectation values, one column per coil, at the recording times.

## Context requirements and limitations

Use `zeeman-hilb` formalism and call the context with the `labframe` assumption. Crystal and powder contexts assemble the anisotropic Hamiltonian term; the liquid context drops that term, including a giant-spin crystal field. For a powder context, set `parameters.sum_up=false` because the result is a structure. Additional rotating frames (`parameters.rframes`) and frequency offsets (`parameters.offset`) are unsupported because the field operator is added in the laboratory frame. `sys.magnet` must be 1 Tesla so `parameters.hzeeman` is per Tesla. The bath temperature is `inter.temperature`.
