# kernel/contexts/crystal.m

- Signature: `answer=crystal(spin_system,pulse_sequence,parameters,assumptions)`

## Purpose

Builds the Hamiltonian and relaxation/kinetics inputs for a single-crystal orientation, then calls the supplied pulse-sequence function handle. The orientation is specified by three Euler angles.

## Parameters / inputs

- `pulse_sequence` — function handle for a pulse sequence in the experiments directory.
- `assumptions` — string passed to `assume.m` when the Hamiltonian is built.
- `parameters.spins` — cell array of spin species used by the pulse sequence, for example `{'1H','13C'}`.
- `parameters.offset` — transmitter offsets in Hz, one per spin in `parameters.spins`.
- `parameters.orientation` — row vector of three Euler angles in radians, specifying the system orientation relative to the input orientation.
- `parameters.rframes` — rotating-frame specification. The source gives `{{'13C',2},{'14N,3}}` as an example of second-order carbon-13 and third-order nitrogen-14 transformations. When used, the assumptions for those spins should be in the laboratory frame.
- `parameters.needs` — cell array of additional sequence requirements:
  - `'zeeman_op'` requests the laboratory-frame Zeeman Hamiltonian in `parameters.hzeeman`.
  - `'aniso_eq'` requests thermal equilibrium recomputed from the full anisotropic Hamiltonian at the current orientation and supplied as `parameters.rho0`.
- Other `parameters` subfields may be required by the pulse sequence; consult its documentation.

The wrapper also sets `parameters.spc_dim` to 1 and `parameters.spn_dim` to the spin-dynamics matrix dimension before calling the sequence.

## Output

Returns whatever the pulse sequence returns.

## Note

Arbitrary-order rotating-frame transformations, including infinite order, are supported. See the header of `rotframe.m` for details.