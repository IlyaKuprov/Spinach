# kernel/contexts/liquid.m

- Signature: `answer=liquid(spin_system,pulse_sequence,parameters,assumptions)`
- Source: [kernel/contexts/liquid.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/contexts/liquid.m)
- Wiki: [liquid.m](https://spindynamics.org/wiki/index.php?title=liquid.m)

## Contract

`liquid` constructs the liquid-phase spin Liouvillian components and calls `pulse_sequence(spin_system,parameters,H,R,K)`. The context returns whatever the pulse sequence returns. It applies the supplied `assumptions` before constructing the operators; the source documentation lists context-specific choices such as `nmr`, `epr`, and `labframe`, with the pulse-sequence documentation determining which are appropriate.

This context has no orientation or spatial grid: `parameters.spc_dim=1` and `parameters.spn_dim=size(H,1)`. The spin-space operator dimension follows the spin system's active basis and assumptions. The context does not define the shape of the pulse sequence's returned value.

## Channel and sequence inputs

`parameters.spins` is a nonempty cell array of isotope labels in channel order; the documented example is `{'1H','13C'}`. `parameters.offset` gives one transmitter offset per listed spin, in Hz, and defaults to zero offsets when omitted. `parameters.needs` requests optional sequence inputs:

- `zeeman_op` builds the laboratory-frame Zeeman operator and places it in `parameters.hzeeman`.
- `rho_eq` builds the thermal-equilibrium state with respect to the isotropic Hamiltonian and places it in `parameters.rho0`.
- `rdc` selects residual-dipolar-coupling handling. The context uses the order matrix through `residual(spin_system)`, then forms the coherent Hamiltonian from the isotropic part and the rank components at zero orientation. Relaxation and kinetics are also built for this mode.

The `parameters.rframes` cell array specifies rotating frames as {isotope, order} pairs. For example, `{{'13C',2},{'14N',3}}` requests a second-order transformation for the carbon-13 carrier and a third-order transformation for nitrogen-14; the source states that arbitrary orders, including infinite order, are supported. If omitted, no extra rotating-frame transformations are applied.

## Example from the source documentation

This context function assembles sequence inputs; it does not itself specify a complete pulse sequence.

## State-dependent chemistry boundary

This context rejects a function handle returned by `kinetics` with `Spinach:liquid:stateDependentKinetics`. Multi-reactant or callback-rate reaction records require a custom pulse sequence using `step`/`iserstep`, rather than static context assembly; see `examples/kinetics/nonlinear/bimolecular_closures.m` and `examples/microfluidics/reacting_flow_nmr.m`. Constant matrix kinetics remain supported.
