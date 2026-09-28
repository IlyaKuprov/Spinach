# kernel/contexts/device.m

- Signature: `answer=device(spin_system,pulse_sequence,parameters,assumptions)`

## Purpose

Builds the evolution generators for a spin-boson device at a fixed orientation of the spin subsystem, then passes them to the supplied pulse-sequence function handle.

## Parameters / inputs

- `pulse_sequence` — function handle for a pulse sequence in the experiments directory.
- `parameters.spins` — cell array of spin species used as channels, in channel order (for example, `{'E'}`). It may be omitted when no spin channels are needed.
- `parameters.offset` — transmitter offsets in Hz, one per listed spin species.
- `parameters.mode_offset` — detuning offsets in Hz, one per bosonic mode in declaration order. The transmitter sign convention applies: each offset contributes minus the offset times that mode's number operator, so a positive offset lowers the mode frequency.
- `parameters.decouple` — cell array of spin species to remove from the evolution generators and initial state (for example, `{'1H'}`); defaults to an empty cell array.
- `parameters.orientation` — Euler angles in radians using the active ZYZ convention, specifying the spin-subsystem orientation. Bosonic terms are unaffected; the default is `[0 0 0]`.
- `parameters.rframes` — numerical rotating-frame specification for spin species, for example `{{'E',2}}`; see the header of `rotframe.m`.
- `parameters.needs` — cell array of additional sequence requirements. `'rho_eq'` requests thermal equilibrium at the system temperature, including Bose-Einstein populations of the bosonic modes, placed in `parameters.rho0`.
- Other `parameters` subfields may be required by the pulse sequence; consult its documentation.
- `assumptions` — one of `'labframe'`, `'cavity'`, or `'spin-phonon'`; see the header of `assume.m`.

The wrapper sets `parameters.spc_dim` to 1 and `parameters.spn_dim` to the spin-dynamics matrix dimension before calling the sequence.

## Output

Returns whatever the pulse sequence returns.

## Notes

- The spin system must contain at least one bosonic mode. Pure-spin systems should use `liquid.m`, `crystal.m`, or `powder.m`.
- Dissipative bosonic modes require a Liouville-space formalism; coherent simulations may also use `zeeman-hilb`.