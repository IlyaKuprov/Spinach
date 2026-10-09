# kernel/contexts/crystal.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/contexts/crystal.m) · [Spin Dynamics Wiki: crystal.m](https://spindynamics.org/wiki/index.php?title=crystal.m)

## Contract

`answer=crystal(spin_system,pulse_sequence,parameters,assumptions)` assembles the Hamiltonian, relaxation operator, and kinetic operator for one single-crystal orientation, then calls the supplied pulse-sequence function handle. It does not sweep an orientation grid or perform powder averaging.

## Orientation, operators, and dimensions

`parameters.orientation` is a real, finite three-element Euler-angle vector in radians, specifying the system orientation relative to its input orientation. The context evaluates the Hamiltonian as the isotropic part plus the rotated anisotropic part, `H=I+orientation(Q,parameters.orientation)`, and Hermitian-symmetrises it. It builds relaxation at the same orientation and obtains kinetics from the spin system. The angular coordinates affect the anisotropic spin terms; the isotropic contribution is not rotated.

The spatial subspace has dimension one: the sequence receives `parameters.spc_dim=1`. The context sets `parameters.spn_dim=size(H,1)`, the Hamiltonian matrix dimension. These are context metadata, not a promise about the sequence's return value: the function returns whatever `pulse_sequence` returns. Its inputs are `spin_system`, the updated `parameters`, `H`, `R`, and `K`; it owns the observable and result shape.

## Inputs and conventions

- `pulse_sequence` must be a function handle. `assumptions` is a character string passed to `assume.m` before the Hamiltonian is built, so the operator basis/frame is the one requested by that assumption.
- `parameters.spins` lists channel species, for example `{'1H','13C'}`; `parameters.offset` gives the corresponding transmitter offsets in Hz. Missing offsets default to zero.
- `parameters.rframes` selects rotating-frame transformations by species and order. The source example uses second order for carbon-13 and third order for nitrogen-14; when using these transformations, the assumptions for those spins should be laboratory frame. Arbitrary order, including infinite order, is supported by the called rotating-frame routine.
- Put `'zeeman_op'` in `parameters.needs` to receive the laboratory-frame Zeeman Hamiltonian in `parameters.hzeeman`. Put `'aniso_eq'` there to recompute thermal equilibrium with the full anisotropic Hamiltonian at this orientation and receive it as `parameters.rho0`.
- Missing `decouple` and `rframes` fields default to empty; missing channel offsets default to zero. Other sequence-specific parameter fields pass through.

The source states the angle and offset units explicitly: radians for orientation, hertz for transmitter offsets. The Hamiltonian is assembled by Spinach's `hamiltonian` and `orientation` routines in the formalism selected through the spin-system assumptions; this context adds no separate unit conversion for the spin-system Hamiltonian.

## Source-supported examples

Use a single orientation such as `[0 0 0]` for the input orientation. The header's channel example is `{'1H','13C'}`; its rotating-frame example specifies order 2 for carbon-13 and order 3 for nitrogen-14. Those examples describe channel/frame choices, not simulated outputs.

## State-dependent chemistry boundary

This context rejects a function handle returned by `kinetics` with `Spinach:crystal:stateDependentKinetics`. Multi-reactant or callback-rate reaction records require a custom pulse sequence using `step`/`iserstep`, rather than static context assembly; see `examples/kinetics/nonlinear/bimolecular_closures.m` and `examples/microfluidics/reacting_flow_nmr.m`. Constant matrix kinetics remain supported.
