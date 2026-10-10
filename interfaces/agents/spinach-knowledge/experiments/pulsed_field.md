# experiments/pulsed_field.m

Source: [experiments/pulsed_field.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/pulsed_field.m)
Spinach Wiki: [pulsed_field.m](https://spindynamics.org/wiki/index.php?title=pulsed_field.m)

- Signature: `answer=pulsed_field(spin_system,parameters,H,R,K)`

## Purpose

Simulates magnetisation dynamics under a time-dependent laboratory-frame Z field with spin-phonon relaxation, as used in pulsed-field magnetometry of molecular magnets. It is a field-sweep propagator, not a generic pulse-sequence builder: it replaces the field profile by constant-field stairs and rebuilds the phonon dissipator for each stair Hamiltonian.

## Inputs and frame

- `parameters.field_prof` is a function handle returning a real finite field in Tesla for a time in seconds. The Hamiltonian for stair `n` uses the profile at its midpoint, `(n-0.5)*timestep`; recorded fields are evaluated at the stair endpoints.
- `parameters.hzeeman` is the Hilbert-space Zeeman operator per Tesla, in rad/s/T. `H` is the context Hamiltonian containing the Zeeman term at `sys.magnet=1` Tesla; the routine subtracts `hzeeman` and then adds the requested field on each stair. For `crystal` or `powder` contexts, set `sys.magnet=1` and `parameters.needs={'zeeman_op'}`: the contexts supply `hzeeman` only when requested, and the `pulsed_field` grumbler requires both this operator and the 1 T reference field.
- `parameters.timestep` is the stair width in seconds, `parameters.nsteps` is the number of stairs, and `parameters.nout` is the number of stairs between records.
- `parameters.coil` is one Hermitian Hilbert-space observable or a cell array of them. `parameters.phonon_x`, `parameters.phonon_i0`, and `parameters.phonon_alpha` provide the spin-phonon coupling operator and spectral-density parameters; the exponent must be at least 1. The phonon-bath temperature is `spin_system.rlx.temperature`.
- The signature also accepts `R` and `K`, but this implementation does not use them; it constructs its own spin-phonon dissipator.

The calculation requires the `zeeman-hilb` formalism and a context call with the `labframe` assumption. Crystal and powder contexts supply the anisotropic Hamiltonian; the liquid context omits it. Use `parameters.sum_up=false` for the powder context so the result remains a structure. Additional rotating frames and frequency offsets are unsupported.

## Propagation and output

The initial density matrix is thermal equilibrium of the field-free Hamiltonian obtained by removing the 1-T Zeeman term, with populations proportional to `exp(-hbar*(E-min(E))/(kbol*T))`.  On each stair the routine diagonalises the instantaneous Hamiltonian, builds the phonon relaxation superoperator in that eigenbasis, applies an exact coherent half-stair phase, a second-order Taylor step for the dissipator, and the second coherent half-stair phase. The dissipator-times-stair-width must remain small; the coherent part is treated exactly for any stair width. The implementation uses Hilbert-space matrix products, with cubic cost in Hilbert-space dimension.

The returned `answer` structure contains `t` (record times in seconds), `field` (profile values in Tesla at those times), and `obs` (real observable expectations, one column per coil). Records are made every `nout` stairs, so the allocated number is `floor(nsteps/nout)`.
