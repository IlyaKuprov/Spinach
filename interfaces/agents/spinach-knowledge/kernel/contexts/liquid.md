# kernel/contexts/liquid.m

- Signature: `answer=liquid(spin_system,pulse_sequence,parameters,assumptions)`

## Purpose

Builds the Liouvillian components for a liquid-phase simulation and passes them to a pulse-sequence function handle. With `'rdc'` requested, it processes residual anisotropic couplings using the supplied order matrix.

## Inputs

- `spin_system` — Spinach spin system.
- `pulse_sequence` — function handle called as `pulse_sequence(spin_system,parameters,H,R,K)`. See the experiments directory for pulse sequences supplied with Spinach.
- `parameters.spins` — non-empty cell array of isotope strings present in the spin system, in channel order (for example, `{'1H','13C'}`). Required.
- `parameters.offset` — numeric array of transmitter offsets, with one element per entry in `parameters.spins`. Defaults to zero offsets.
- `parameters.needs` — cell array containing any of `'rdc'`, `'zeeman_op'`, or `'rho_eq'`. Defaults to `{}`.
  - `'rdc'` requests residual anisotropic coupling processing.
  - `'zeeman_op'` places the laboratory-frame Zeeman Hamiltonian in `parameters.hzeeman`.
  - `'rho_eq'` places the thermal equilibrium state in `parameters.rho0`, computed for the isotropic Hamiltonian at the specified temperature.
- `parameters.rframes` — cell array of `{isotope,order}` pairs specifying rotating-frame transformations. The isotope must be present in the spin system. For example, `{{'13C',2},{'14N',3}}` requests second-order rotating-frame transformation for carbon-13 and third-order transformation for nitrogen-14. Defaults to `{}`. When used, assumptions for the respective spins should be laboratory frame. Arbitrary transformation order, including infinite order, is supported; see the header of `rotframe.m`.
- `parameters.decouple` — defaults to `{}` if absent.
- `parameters.*` — additional fields may be required by the pulse sequence; consult its documentation.
- `assumptions` — context-specific assumptions, such as `'nmr'`, `'epr'`, or `'labframe'`; see the pulse-sequence header. Must be a character string.

## Processing

The interface applies the assumptions and obtains relaxation `R` and kinetics `K`. In RDC mode, it computes these before liquid-crystal averaging, then obtains the coherent Hamiltonian as `H=I+orientation(Q,[0 0 0])`. Otherwise it obtains the isotropic Hamiltonian, relaxation, and kinetics directly. If requested, it also constructs `parameters.hzeeman` and/or `parameters.rho0`. It applies channel offsets and the requested rotating-frame transformations before calling the pulse sequence. The call receives `parameters.spc_dim=1` and `parameters.spn_dim=size(H,1)`.

## Output

- `answer` — whatever the pulse sequence returns.

[liquid.m documentation](https://spindynamics.org/wiki/index.php?title=liquid.m)