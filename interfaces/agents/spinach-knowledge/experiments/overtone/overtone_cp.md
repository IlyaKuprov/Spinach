# experiments/overtone/overtone_cp.m

Source: [canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/overtone/overtone_cp.m) · [Spinach wiki](https://spindynamics.org/wiki/index.php?title=overtone_cp.m)

Signature: `spectrum=overtone_cp(spin_system,parameters,H,R,K)`.

## Behaviour

This wrapper calculates the overtone reference frequency as `-2*spin(parameters.spins{1})*spin_system.inter.magnet/(2*pi)`, lifts `parameters.Nx` and `parameters.Hx` over the requested Fokker–Planck spatial dimension using an identity Kronecker factor, applies one spin-lock pulse step, then calls `overtone_a` for acquisition. The source labels `Nx` as the X Zeeman operator of the quadrupolar nucleus and `Hx` as the X Zeeman operator of a spin-1/2 nucleus. It weights those channels with the first and second entries of `rf_pwr`, respectively; it does not specify a source-to-destination magnetisation direction or a protein-specific transfer pathway. The existing background statement that quadrupolar nuclei have spin greater than 1/2 is general context; this wrapper only labels the channel as quadrupolar and does not validate spin or construct an electric-field-gradient tensor.

With `method='average'`, the code forms an average pulse operator using `omega=2*pi*(ovt_frq-rf_frq)` and applies its propagator to `rho0` for `rf_dur`. With `method='fplanck'`, it instead calls `shaped_pulse_af` using `H+rf_pwr(2)*Hx`, `Nx`, the initial state, the offset in Hz, the quadrupolar-channel power, and the pulse duration. These are distinct source branches; the wrapper exposes no separate MAS-rate or sample-orientation field. `spc_dim` is documented as a Fokker–Planck spatial dimension, not an angle.

## Inputs and limits

Required fields are `spins` (one-element cell array), `spc_dim` (positive integer), `method` (`'average'` or `'fplanck'`), `sweep` (two-element real frequency extent in Hz around the overtone frequency), `npoints` (positive integer), `rho0` (initial state), `coil` (detection state), `Nx`, `Hx`, `rf_frq` (spin-lock offset in Hz), `rf_pwr` (two powers in rad/s, quadrupolar then spin-1/2), and `rf_dur` (positive scalar in seconds). The consistency check requires `H`, `R`, and `K` to be numeric matrices; the source does not impose an equal-size check here.

The source comment warns that relaxation must be present for the matrix inversion in `overtone_a` to converge and that `R` must not be thermalised. This documents an input requirement, not a runtime result. The return value is the result from `overtone_a`; this wrapper does not reshape it. `sweep` and `npoints` specify the requested frequency interval and sampling, while the precise MATLAB output shape is defined by the acquisition helper.
