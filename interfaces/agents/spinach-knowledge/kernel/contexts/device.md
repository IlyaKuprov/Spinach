# kernel/contexts/device.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/contexts/device.m) · [Spin Dynamics Wiki: device.m](https://spindynamics.org/wiki/index.php?title=device.m)

## Contract

`answer=device(spin_system,pulse_sequence,parameters,assumptions)` constructs the Hamiltonian, relaxation operator, and kinetic operator for a spin-boson device, then passes them to the supplied pulse-sequence function handle. The system must contain at least one bosonic mode; the spin subsystem is evaluated at one fixed orientation, while bosonic terms are not rotated.

## Degrees of freedom and dimensions

Let `D=size(H,1)` for the Hamiltonian assembled in the system's composite spin-boson basis. The context sets `parameters.spn_dim=D` and `parameters.spc_dim=1`; there is no spatial-orientation grid in this context. It calls the sequence with `spin_system`, `parameters`, `H`, `R`, and `K`, and returns whatever the sequence returns, so the final output shape is sequence-specific.

The orientation is `parameters.orientation`, a real three-element vector of active ZYZ Euler angles in radians for the spin subsystem only. Its default is `[0 0 0]`. The context evaluates `H=I+orientation(Q,parameters.orientation)` and Hermitian-symmetrises the result; it also constructs relaxation at that orientation and obtains kinetics for the complete system.

## Modes, channels, and units

The system includes spin and bosonic particles. The context accepts spin-channel species only for spin particles; the source example is `parameters.spins={'E'}`, and this field may be omitted when no spin channels are needed. Spin transmitter offsets in `parameters.offset` are in Hz and follow the transmitter sign convention.

`parameters.mode_offset` supplies one detuning in Hz per bosonic mode, in declaration order. For mode `n`, the context subtracts `2*pi*mode_offset(n)*N_n` from `H`, where `N_n` is that mode's number operator; a positive detuning therefore lowers the mode frequency. The default mode offsets are zero. The context's orientation and offsets have these explicit units; other Hamiltonian terms retain the units and operator basis supplied by Spinach's system and Hamiltonian construction.

`assumptions` must select one of `'labframe'`, `'cavity'`, or `'spin-phonon'`; the context applies it through `assume.m` before constructing operators. `parameters.rframes` requests rotating frames for spin species (the source example is `{{'E',2}}`); it does not rotate bosonic terms. `parameters.decouple` selects spin decoupling and defaults to empty.

## Optional equilibrium state and formalism

Put `'rho_eq'` in `parameters.needs` to calculate the thermal equilibrium state at the system temperature, including Bose-Einstein populations of the bosonic modes; it is passed as `parameters.rho0`. Dissipative bosonic modes require a Liouville-space formalism. Coherent simulations may also use `zeeman-hilb`, as stated in the source header.

## Source-supported example

A spin-channel list may be `{'E'}`. For a system with multiple modes, provide one `mode_offset` value for each mode in its declaration order; use zero values for no mode detuning. These examples specify parameter structure, not a calculated device response.
