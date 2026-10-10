# experiments/rapidscan.m

Source: [experiments/rapidscan.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/rapidscan.m)
Spinach Wiki: [rapidscan.m](https://spindynamics.org/wiki/index.php?title=rapidscan.m)

- Signature: `[b_axis,spectrum]=rapidscan(spin_system,parameters)`

## Purpose

Simulates a time-domain rapid field-scan ESR experiment (Eaton-style). It constructs the microwave-frame generator and propagates the equilibrium density matrix while stepping through a prescribed magnetic-field sweep; it does not define or compile a pulse sequence. Call it directly, without a context function.

## Inputs

- `parameters.mw_pwr` is a non-negative scalar in rad/s used as the microwave-drive coefficient.
- `parameters.sweep` is a two-element vector of field offsets in Tesla, swept linearly around the centre field `spin_system.inter.magnet`.
- `parameters.nsteps` is a positive integer number of field points; `parameters.timestep` is a positive duration in seconds.

## Calculation and output

The routine starts from isotropic thermal equilibrium. It forms the electron microwave operator from `L+`, builds the laboratory-frame Zeeman and coupling Hamiltonians and relaxation superoperator, symmetrises the drift Hamiltonian, and rotates it to the electron microwave frame using the carrier. The microwave term is `H_mw=-mw_pwr*(Ep-Ep')/(2i)`; the generator adds this drive and `1i*R` to the rotating-frame drift.

The sweep offsets are `linspace(sweep(1),sweep(2),nsteps)`; adding `spin_system.inter.magnet` gives `b_axis` in Tesla. The Zeeman Hamiltonian is normalised by the centre field, then each point samples the current `L+` expectation and advances the state for one timestep under the offset-field generator using `step`. The outputs are a Tesla field-axis column and the corresponding complex `L+` amplitudes in `spectrum`.
