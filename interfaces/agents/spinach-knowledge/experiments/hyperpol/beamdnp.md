# experiments/hyperpol/beamdnp.m

- Signature: contact_curve=beamdnp(spin_system,parameters,H,R,K)

## Purpose and physical scope

Simulates the BEAM dynamic nuclear polarisation (DNP) pulse-block experiment described in Science Advances, DOI 10.1126/sciadv.abq0536. The reported curve is calculated from the supplied spin system and context matrices; it is not a measured contact curve. Electron-nuclear hyperfine coupling, relaxation, and exchange or other kinetics can contribute only when represented in the supplied H, R and K.

## Inputs and parameters

H, R and K are context-supplied Hamiltonian, relaxation, and kinetics matrices. The routine constructs their combined generator as H+1i*R+1i*K and explicitly constructs electron control operators for spin E.

- parameters.irr_powers: microwave amplitude, interpreted as the electron nutation frequency in Hz.
- parameters.rho0: initial state.
- parameters.coil: detection state.
- parameters.pulse_dur: two pulse durations in seconds, for the positive- and negative-x irradiation steps respectively.
- parameters.nloops: positive integer number of repeated BEAM blocks.

## Pulse sequence and returned signal

First, the routine applies an electron pi/2 flip about y under microwave irradiation. Its duration is computed from the supplied frequency as 1/(4*irr_powers) seconds. Each repeated block then applies irradiation along +x for pulse_dur(1), followed by irradiation along -x for pulse_dur(2), with the corresponding propagators composed according to the selected Spinach formalism. After the initial state and after every block, the routine evaluates hdot(coil,rho).

contact_curve is a vector of nloops+1 contact values: element 1 is the initial coil expectation, and each subsequent element follows one complete two-pulse block. Its natural horizontal coordinate is block number (initial point, then 1 through nloops); the function does not return a continuous-time or physical-time axis. The source supports zeeman-hilb, zeeman-liouv, and sphten-liouv propagation.

## Limits and interpretation

This implements the fixed BEAM pulse sequence and blockwise detection; it is not a general DNP optimiser or a complete model of polarisation transfer independent of the supplied spin Hamiltonian and relaxation/kinetics. The source specifies the DOI above but gives no numeric microwave frequency, pulse-duration pair, loop count, or validated experimental fit, so no such values are invented here.

Paper: https://doi.org/10.1126/sciadv.abq0536
Source: https://github.com/IlyaKuprov/Spinach/blob/main/experiments/hyperpol/beamdnp.m
